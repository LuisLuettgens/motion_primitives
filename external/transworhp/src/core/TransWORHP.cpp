#include "TransWORHP.h"

#include "TWproblem.h"
#include "butcher.h"
#include "TWfolder.h"
#include "twstatus.h"

#ifdef TRANSWORHP_GRAPHICS
#include "Viewer.h"
#else

using GLfloat = float;
void glColor4fv(const GLfloat */*v*/) {}

namespace tw {

int SDL_GetTicks() {
	return 0;
}
FakeThread thethread;
FakeThread *thethread0 = &thethread;

}

#endif

namespace tw {

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::stringstream;


// Ausgabe spaeter Initialisieren mit Zeiger auf diese Funktion.
void worhpoutput(int i, const char *message) {

	static_assert(WORHP_MAJOR > 1 || (WORHP_MAJOR == 1 &&
	                                  (WORHP_MINOR > 9 || (WORHP_MINOR == 9 && WORHP_PATCH != 0))),
	              "WORHP 1.9.1 or newer required!");

	Status flag = Status::NORMAL;

	if (i & WORHP_PRINT_WARNING) {
		flag = Status::WARN;
	}
	else if (i & WORHP_PRINT_ERROR) {
		flag = Status::ERR;
	}
	else if (i & WORHP_PRINT_GREEN) {
		flag = Status::GREEN;
	}
	else if (i & WORHP_PRINT_BLUE) {
		flag = Status::BLUE;
	}

	if (i & WORHP_PRINT_CONTINUE_NEXT) {
		flag = flag | Status::CONT_NEXT;
	}
	if (i & WORHP_PRINT_CONTINUE_PREV) {
		flag = flag | Status::CONT_PREV;
	}

	if (i & WORHP_PRINT_BOLD) {
		flag = flag | Status::BOLD;
	}

	MyStatus("WORHP", message, flag);
}

extern "C" {

#ifndef worhp_print_t
	using worhp_print_t = void (*) (int mode, const char s[]);
#endif // worhp_print_t

	void SetWorhpPrint(worhp_print_t f);
}


TransWorhp::TransWorhp(TransWorhpProblem *ph, TWparameter *twparam)
	: phase(ph),
	n_dis(ph->n_dis),
	n_ode(ph->n_ode), n_ctrl(ph->n_ctrl), n_param(ph->n_param), n_rand(ph->n_rand), n_neben(ph->n_neben), n_integral(ph->n_integral),
	Delta1(0), Delta2(0),
	n_var(0), n_con(0), n_zen(ph->n_zen),
	X(nullptr), G(nullptr), ZEN(nullptr),
	X_low(nullptr), X_upp(nullptr), G_low(nullptr), G_upp(nullptr),
	Lambda(nullptr), Mu(nullptr),
	viewer(nullptr),
	twparameter(twparam),
	twdiscretization(nullptr),
	butcher(n_ode, n_ctrl, twparameter->abserr, twparameter->relerr),
	transworhp_type(TransWORHP_type::unknown),
	DF_start(0), DF_nnz(0),
	DG_start(0), DG_nnz(0),
	HM_start(0), HM_start_diag(0), HM_nnz_ohne_diag(0),
	eps(twparameter->eps),
	USERDF(0), USERDG(0), USERHM(0),
	tmp_ode_1(n_ode), tmp_ode_2(n_ode),
	tmp_ctrl_1(n_ctrl), tmp_ctrl_2(n_ctrl),
	tmp_rand_1(n_rand), tmp_rand_2(n_rand),
	tmp_neben_1(n_neben), tmp_neben_2(n_neben),
	tmp_integral_1(n_integral), tmp_integral_2(n_integral), tmp_integral_12(n_integral)
{

	SetWorhpPrint(worhpoutput);

	T.reserve(n_dis);
	for (int i = 0; i < n_dis; i++) {
		T.push_back( i/(n_dis-1.) );
	}

	// Groesse besser anpassen (auch in twfolder::MeshRef() und tw::lagrange())//TODO
	lagrange_integral = vector<double>(n_integral*n_dis*2);
	lagrange_weight = vector<double>(n_integral, 1.0); //filled with 1.0
}

/**
* Destruktor
*/
TransWorhp::~TransWorhp() {

}


void TransWorhp::TimeAxis(double exponent) {
	for (int i = 0; i < n_dis; i++) {
		T[i] = std::pow(i/(n_dis-1.), exponent);
	}
}


void TransWorhp::LinearTimeAxis(double start, double end) {
	phase->t0 = start;
	phase->tf = end;
	for (int i = 0; i < n_dis; i++) {
		T[i] = start + i/(n_dis-1.) * (end-start);
	}
}


void TransWorhp::boxConToInitGuess() {
	// damit die Box-Schranken durch den Startwert nicht verletzt sind
	for (int i = 0; i < n_var; i++) {
		if (X_low[i] == X_upp[i]) {
			X[i] = X_low[i];
		}
	}
}


void TransWorhp::printStructure() const {

	if (DS_obj.useStructure() && twparameter->showOBJstructure) {
		stringstream s;
		s << "Structure of objective" << endl;
		s << DS_obj;
		MyStatus("TransWORHP", s.str(), Status::NORMAL);
	}

	if (DS_ode.useStructure() && twparameter->showODEstructure) {
		stringstream s;
		s << "Structure of ode system" << endl;
		s << DS_ode;
		MyStatus("TransWORHP", s.str(), Status::NORMAL);
	}

	if (DS_rand.useStructure() && twparameter->showRANDstructure) {
		stringstream s;
		s << "Structure of boundary constraints" << endl;
		s << DS_rand;
		MyStatus("TransWORHP", s.str(), Status::NORMAL);
	}

	if (DS_neben.useStructure() && twparameter->showNEBENstructure) {
		stringstream s;
		s << "Structure of state constraints" << endl;
		s << DS_neben;
		MyStatus("TransWORHP", s.str(), Status::NORMAL);
	}
}


void TransWorhp::Constraints() {
	 Constraints2(G);
}


int TransWorhp::HM_structure(int hessianstructure, WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix *HM, int offset) {

	int ind = HM_structure_ohne_Diag(hessianstructure, DF, DG, HM, offset);

	// Diagonale
	for (int i = 0; i < n_var; ++i) {
		if (HM) {
			HM->row[ind] = i + 1;
			HM->col[ind] = i + 1;
		}
		ind++;
	}

	return ind;
}


void TransWorhp::HM_calculate(int hessianvalues, WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM, double ScaleObj, double *Mu) {
	//cout << "HM_calc" << endl;
	if (hessianvalues == 0) {
		HM_calculate1(DF, DG, HM, ScaleObj, Mu);
	} else if (hessianvalues == 1) {
		HM_calculate2(DF, DG, HM, ScaleObj, Mu);
	} else if (hessianvalues == 2) {
		HM_calculate3(DF, DG, HM, ScaleObj, Mu);
	} else if (hessianvalues == 3) {
		HM_calculate4(DF, DG, HM, ScaleObj, Mu);
	} else if (hessianvalues == 4) {
		HM_calculate5(DF, DG, HM, ScaleObj, Mu);
	}
}


void TransWorhp::HM_calculate2(WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM, double ScaleObj, double *Mu) {
	HM_calculate1(DF,DG,HM,ScaleObj,Mu);
}

void TransWorhp::HM_calculate3(WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM, double ScaleObj, double *Mu) {
	HM_calculate1(DF,DG,HM,ScaleObj,Mu);
}

void TransWorhp::HM_calculate4(WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM, double ScaleObj, double *Mu) {
	HM_calculate1(DF,DG,HM,ScaleObj,Mu);
}

void TransWorhp::HM_calculate5(WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM, double ScaleObj, double *Mu) {
	HM_calculate1(DF,DG,HM,ScaleObj,Mu);
}


int TransWorhp::DoubleFrom(TransWorhp *ph) {

	cout << n_dis << " " << n_ode << endl;
	for (int i=0; i<n_dis; i++) {
		for (int j=0; j<n_ode; j++) {
			if (i%2) X[x_index(i,j)] = ph->X[ph->x_index(i/2,j)];
			else X[x_index(i,j)] = (ph->X[ph->x_index(i/2,j)] + ph->X[ph->x_index(i/2+1,j)] )*.5;
		}
		for (int j=0; j<n_ctrl; j++) {
			if (i%2) X[u_index(i,j)] = ph->X[ph->u_index(i/2,j)];
			else X[u_index(i,j)] = (ph->X[ph->u_index(i/2,j)] + ph->X[ph->u_index(i/2+1,j)] )*.5;
		}
	}
	for (int i=0; i<phase->n_param; i++) {
		X[p_index(i)] = ph->X[ph->p_index(i)];
	}
	return 0;
}

void TransWorhp::ToMATLAB(const std::string& filename) const {

	std::ofstream of(filename);

	of.setf(std::ios::scientific);
	of << std::setprecision(9);

	of << "% TransWORHP-Result" << std::endl;

	of << "% ";
	for (int i = 0; i < n_param; i++) {
		of << std::setw(20) << p(i);
	}
	of << std::endl;

	for (int i = 0; i < n_dis; i++) {

		if (phase->freetime) {
			of << std::setw(20) << (T[i] * p(0));
		} else {
			of << std::setw(20) << T[i];
		}

		for (int j = 0; j < n_ode; j++) {
			of << std::setw(20) << x(i, j);
		}
		for (int j = 0; j < n_ctrl; j++) {
			of << std::setw(20) << u(i, j);
		}
		of << std::endl;
	}
}


void TransWorhp::FromMATLAB(const string& filename) {

	std::ifstream stream(filename);

	if (!stream) {
		MyStatus("FromMATLAB", "There is no file: " + filename, Status::ERR);
		MyStatus("FromMATLAB", "No initial guess is set.", Status::ERR);
		return;
	} else {
		MyStatus("FromMATLAB", "Reading: " + filename, Status::NORMAL);
	}

	fromMATLAB_impl(stream);
}


void TransWorhp::DoubleToMATLAB(double *data, int n, const std::string& filename) {

	std::ofstream of(filename);

	of.setf(std::ios::scientific);
	of << std::setprecision(9);

	for (int j = 0; j < n; j++) {
		of << std::setw(20) << data[j];
		of << endl;
	}

}

void TransWorhp::MatrixToMATLAB(const WorhpMatrix &m, const std::string& filename) {

	std::ofstream of(filename);

	of.setf(std::ios::scientific);
	of << std::setprecision(9);


	if (m.kind==1) {
		for (int j=1;j<=m.nRow;j++) {
			for (int i=1;i<=m.nCol;i++) {

				int found=0;
				for (int z=0;z<m.nnz;z++) {
					if (j==m.row[z] && i==m.col[z]) {
						of << std::setw(20) << m.val[z];
						found=1;
						break;
					}
				}
				if (found==0) of << std::setw(20) << 0;
			}
			of << endl;
		}
	}


	if (m.kind==2) {
		for (int j=1;j<=m.nRow;j++) {
			for (int i=1;i<=m.nCol;i++) {

				int found=0;
				for (int z=0;z<m.nnz;z++) {
					if (j==m.row[z] && i==m.col[z]) {
						of << std::setw(20) << m.val[z];
						found=1;
						break;
					}
				}
				if (found==0) of << std::setw(20) << 0;
			}
			of << endl;
		}
	}

	if (m.kind==5) {
		for (int j=1;j<=m.nRow;j++) {
		//	for (int i=1;i<=m.nCol;i++) {

				int found=0;
				for (int z=0;z<m.nnz;z++) {
					if (j==m.row[z]) {
						of << std::setw(20) << m.val[z];
						found=1;
						break;
					}
				}
				if (found==0) of << std::setw(20) << 0;
		//	}
			of << endl;
		}
	}
}


void TransWorhp::PrintOCPstates(std::ostream *os) const {

	string s("OCP States");
	if (os) {
		*os << "# " << s << endl;
	} else {
		MyStatus("Viewer", s, Status::NORMAL);
	}

	s = "Phase [" + std::to_string(n_ode) + " " + std::to_string(n_ctrl) + " " + std::to_string(n_dis) + " " + std::to_string(n_param) + "] = " + std::to_string(n_var);

	if (os) {
		*os << s << endl;
	} else {
		MyStatus("OCP", s, Status::NORMAL);
	}


	if (X) {

		stringstream b;

		b << std::setw(13) << "T";
		for (int i = 0; i < n_ode; i++) {
			b << std::setw(13) << "X" << i;
		}
		for (int i = 0; i < n_ctrl; i++) {
			b << std::setw(13) << "U" << i;
		}

		if (os) {
			*os << b.str() << endl;
		} else {
			MyStatus("OCP", b.str(), Status::NORMAL);
		}

		for (int j = 0; j < n_dis; j++) {

			stringstream b;

			b <<  std::setw(13) << T[j];

			for (int i = 0; i < n_ode; i++) {
				b << std::setw(14) << x(j,i);
			}

			for (int i = 0; i < n_ctrl; i++) {
				b << std::setw(14) << u(j,i);
			}

			if (os) {
				*os << b.str() << endl;
			} else {
				MyStatus("OCP", b.str(), Status::NORMAL);
			}
		}

	}


	s = "Objective function = " + std::to_string(phase->obj());

	if (os) {
		*os << s << endl;
	} else {
		MyStatus("OCP", s, Status::NORMAL);
	}
}


void TransWorhp::PrintMultipliers(std::ostream *os) const {

	string s("Print Multipliers");
	if (os) {
		*os << "# " << s << endl;
	} else {
		MyStatus("Viewer", s, Status::NORMAL);
	}


	for (int i = 0; i < n_var; i++) {
		stringstream a;
		a << "Lambda("  << std::setw(3) << i << ")" ;

		stringstream b;
		b << std::setw(20) << Lambda[i];

		if (os) {
			*os << std::setw(20) << a.str() << b.str();
		} else {
			MyStatus(a.str(), b.str(), Status::NORMAL);
		}
	}

	for (int i = 0; i < n_con; i++) {
		stringstream a;
		a << "Mu("  << std::setw(3) << i << ")" ;

		stringstream b;
		b << std::setw(20) << Mu[i];
		if (os) {
			*os << std::setw(20) << a.str() << b.str();
		} else {
			MyStatus(a.str(), b.str(), Status::NORMAL);
		}
	}
}

void TransWorhp::Structure_Sizes(TWfolder *f, int hessianstructure, int &DF_nnz, int &DG_nnz, int &HM_nnz) {

	DF_nnz = DF_structure();
	DG_nnz = DG_structure(f);

	// Ablesen der Groesse fuer HM...
	WorhpMatrix DFtemp,DGtemp;
	ZeroWorhpMatrix(&DFtemp);
	DFtemp.nnz = DF_nnz;
	DFtemp.nRow = n_var;
	DFtemp.nCol = 1;
	DFtemp.kind = 5;

	InitWorhpMatrix(&DFtemp, (char*)"DFtemp", 0,0,0);
	DF_structure(&DFtemp);
	// PrintWorhpMatrix0(&DFtemp);

	ZeroWorhpMatrix(&DGtemp);
	DGtemp.nnz = DG_nnz;
	DGtemp.nRow = n_con;
	DGtemp.nCol = n_var;
	DGtemp.kind = 1;
	InitWorhpMatrix(&DGtemp, (char*)"DGtemp", 0,0,0);
	DG_structure(f,&DGtemp);
	// PrintWorhpMatrix0(&DGtemp);
	// hiermit:

	//MyStatus("Hessian", TWstrings::string_hessianstructure[twparameter.hessianstructure], Status::WARN);
	//MyStatus("Hessian", TWstrings::string_hessianvalues[twparameter.hessianvalues], Status::WARN);

	HM_nnz = HM_structure(hessianstructure, DFtemp,DGtemp);

	FreeWorhpMatrix(&DFtemp);
	FreeWorhpMatrix(&DGtemp);
}


vector<double> TransWorhp::getGrid() const {
	return T;
}


void TransWorhp::newGrid(vector<double> grid) {
	T = move(grid);
}

}
