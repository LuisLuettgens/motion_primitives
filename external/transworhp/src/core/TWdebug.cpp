#include "TWdebug.h"

#include "TransWORHP.h"
#include "ExplTransWORHP.h"
#include "TWproblem.h"
#include "TWfolder.h"

#include "worhp_info.h"

#include <sstream>
#include <algorithm>

using namespace std;

namespace tw {

void TWdebug::PrintDF(TWfolder *twfolder, ostream *os) {

	string s("Check DF");
	if (os) {
		*os << "# " << s << endl;
	} else {
		MyStatus("CHECK", s, Status::NORMAL);
	}

	// temporaere Matrix anlegen und Ableitung auf 2. Weise berechnen
	WorhpMatrix DFtemp;
	ZeroWorhpMatrix(&DFtemp);
	DFtemp.nnz = twfolder->worhp_w.DF.nnz;
	DFtemp.nRow = twfolder->worhp_w.DF.nRow;
	DFtemp.nCol = twfolder->worhp_w.DF.nCol;
	DFtemp.kind = twfolder->worhp_w.DF.kind;
	InitWorhpMatrix(&DFtemp,(char*)"DFtemp", 0,0,0);

	int offset = 0;

	for (TransWorhpProblem* ph : twfolder->phases) {
		offset = ph->solver->DF_structure(&DFtemp, offset);
	}


	for (TransWorhpProblem* ph : twfolder->phases) {
		ph->solver->DF_calculate(DFtemp, twfolder->worhp_w.ScaleObj);
	}


	for (int i = 0; i < twfolder->worhp_w.DF.nnz; i++) {

		stringstream a;

		a << "DF("  << setw(3) << twfolder->worhp_w.DF.row[i]-1;
		a << ")" ;

		stringstream b;
		Status code = matrixentry(b, twfolder->worhp_w.DF.val[i], DFtemp.val[i]);
		if (os) {
			*os << setw(20) << a.str() << b.str();
		} else {
			MyStatus(a.str(),b.str(),code);
		}
	}
}


void TWdebug::PrintDG(TWfolder *twfolder, ostream *os) {

	if (twfolder->phases.front()->solver->transworhp_type == TransWORHP_type::multipleShooting) {
		string s("Check DG - does not work for multiple Shooting!");
		if (os) {
			*os << "# " << s << endl;
		} else {
			MyStatus("CHECK", s, Status::ERR);
		}
		return;
	}

	string s("Check DG");
	if (os) {
		*os << "# " << s << endl;
	} else {
		MyStatus("CHECK", s, Status::NORMAL);
	}

	WorhpMatrix DGtemp;
	ZeroWorhpMatrix(&DGtemp);
	DGtemp.nnz = twfolder->worhp_w.DG.nnz;
	DGtemp.nRow = twfolder->worhp_o.m;
	DGtemp.nCol = twfolder->worhp_o.n;
	DGtemp.kind = 1;
	InitWorhpMatrix(&DGtemp, (char*)"DGtemp", 0,0,0);

	int offset = 0;

	for (TransWorhpProblem* ph : twfolder->phases) {
		offset = ph->solver->DG_structure(twfolder, &DGtemp, offset);
	}

	for (TransWorhpProblem* ph : twfolder->phases) {
		ph->solver->DG_calculate(twfolder, DGtemp);
	}

	for (int i = 0; i < twfolder->worhp_w.DG.nnz; i++) {

		stringstream a, b;

		// -1 da WORHP-Index in FORTRAN
		a << "DG("  << setw(3) << twfolder->worhp_w.DG.col[i]-1;
		a << "/" ;
		a << setw(3) << twfolder->worhp_w.DG.row[i]-1;
		a << ")" ;

		//genauere Angaben welcher Eintrag zu welchem Typ gehoert - ANFANG
		{
		const int row = twfolder->worhp_w.DG.row[i]-1;
		const int col = twfolder->worhp_w.DG.col[i]-1;

		// offset fuer Phase
		int colOff = 0;

		//Phase bestimmen
		size_t phase = 0;
		for (; phase < twfolder->phasenOffset.size()-1; phase++) {
			if (row >= twfolder->phasenOffset[phase] && row < twfolder->phasenOffset[phase+1]) {
				break;
			}
			colOff += twfolder->phases[phase]->solver->n_var;
		}

		if (phase != twfolder->phases.size()) { // keine folder-Beschraenkung (alle anderen Beschraenkungen)
			a << "-> phase" << setw(2) << phase << twfolder->phases.at(phase)->solver->type_G(row - twfolder->phasenOffset.at(phase));

			// Opt.Variable
			int var = col - colOff;

			a << "(";

			// kein Parameter
			if (var < (twfolder->phases[phase]->solver->n_var-twfolder->phases[phase]->n_param)) {
				var %= twfolder->phases[phase]->n_ode+twfolder->phases[phase]->n_ctrl;
				if (var < twfolder->phases[phase]->n_ode) {
					a << "dx";
				} else if (var < twfolder->phases[phase]->n_ode+twfolder->phases[phase]->n_ctrl) {
					var -= twfolder->phases[phase]->n_ode;
					a << "du";
				}
			} else { //Parameter
				var -= twfolder->phases[phase]->solver->n_var-twfolder->phases[phase]->n_param;
				var %= twfolder->phases[phase]->n_param;
				a << "dp";
			}

			a << setw(2) << var << ")";

		} else { // Beschraenkung aus folder
			a << "-> folder " << twfolder->type_G(row - twfolder->phasenOffset.at(phase));

			// erneutes Bestimmen der Phase: noetig, da vorher Phase ueber die Zeile bestimmt wurde
			size_t auxPhase = 0;
			int colOff2 = 0;
			for (; auxPhase < twfolder->phases.size(); auxPhase++) {
				if (col < twfolder->phases[auxPhase]->solver->n_var+colOff2) {
					break;
				}
				colOff2 += twfolder->phases[auxPhase]->solver->n_var;
			}

			int var = col - colOff2;

			a << "(";

			var %= twfolder->phases[auxPhase]->n_ode+twfolder->phases[auxPhase]->n_ctrl;
			if (var < twfolder->phases[auxPhase]->n_ode) {
				a << "dx";
			} else if (var < twfolder->phases[auxPhase]->n_ode+twfolder->phases[auxPhase]->n_ctrl) {
				var -= twfolder->phases[auxPhase]->n_ode;
				a << "du";
			}

			a << setw(2) << var << ")";
		}
		}
		// ENDE

		Status code = matrixentry(b, twfolder->worhp_w.DG.val[i], DGtemp.val[i]);
		if (os) {
			*os << setw(20) << a.str() << b.str();
		} else {
			MyStatus(a.str(),b.str(),code);
		}
	}
}



void TWdebug::PrintHM(TWfolder *twfolder, ostream *os) {

	string s("Check HM");
	if (os) {
		*os << "# " << s << endl;
	} else {
		MyStatus("CHECK", s, Status::NORMAL);
	}

	WorhpMatrix HMtemp;
	ZeroWorhpMatrix(&HMtemp);
	HMtemp.nnz = twfolder->worhp_w.HM.nnz;
	HMtemp.nRow = twfolder->worhp_o.m;
	HMtemp.nCol = twfolder->worhp_o.n;
	HMtemp.kind = 2;
	InitWorhpMatrix(&HMtemp, (char*)"HMtemp", 0,0,0);

	// weglassen!?
	twfolder->HM_structure(twfolder->worhp_w.DF, twfolder->worhp_w.DG, HMtemp);

	for (TransWorhpProblem* ph : twfolder->phases) {
		ph->solver->HM_calculate(twfolder->twparameter->hessianvalues, twfolder->worhp_w.DF, twfolder->worhp_w.DG, HMtemp, twfolder->worhp_w.ScaleObj, twfolder->worhp_o.Mu);
	}

	for (int i = 0; i < twfolder->worhp_w.HM.nnz; i++) {
		stringstream a;

		a << "HM("  << setw(3) << twfolder->worhp_w.HM.col[i]-1;
		a  << "/" ;
		a << setw(3) << twfolder->worhp_w.HM.row[i]-1;
		a << ")" ;

		stringstream b;
		Status code = matrixentry(b, twfolder->worhp_w.HM.val[i], HMtemp.val[i]);
		if (os) {
			*os << setw(20) << a.str() << b.str();
		} else {
			MyStatus(a.str(),b.str(),code);
		}
	}

	FreeWorhpMatrix(&HMtemp);
}


void TWdebug::PrintNLPconstraints(TWfolder *twfolder, ostream *os) {

	string s("NLP Constraints");
	if (os) {
		*os << "# " << s << endl;
	} else {
		MyStatus("CHECK", s, Status::NORMAL);
	}

	stringstream b;

	for (auto it = twfolder->phases.begin();it!=twfolder->phases.end();it++) {

		if (twfolder->phases.size()> 1) {
			b.str("");
			b << "Phase " << it-twfolder->phases.begin();
			if (os) {
				*os << setw(20) << " " << b.str() << endl;
			} else {
				MyStatus(" ", b.str(), Status::NORMAL);
			}
			PrintNLPconstraintsX(twfolder, *it, os, it-twfolder->phases.begin());
		} else {
			PrintNLPconstraintsX(twfolder, *it, os);
		}

	}

	for (auto it = twfolder->phases.begin();it!=twfolder->phases.end();it++) {

		if (twfolder->phases.size() > 1) {
			b.str("");
			b << "Phase " << it-twfolder->phases.begin();
			if (os) {
				*os << setw(20) << " " << b.str() << endl;
			} else {
				MyStatus(" ", b.str(), Status::NORMAL);
			}
			PrintNLPconstraintsG(twfolder, *it, os, it-twfolder->phases.begin());
		} else {
			PrintNLPconstraintsG(twfolder, *it, os);
		}
	}

	if (twfolder->phases.size() > 1) {
		if (os) {
			*os << setw(20) << "Folder" << endl;
		} else {
			MyStatus(" ", "Folder", Status::NORMAL);
		}
	}

	PrintNLPconstraintsG(twfolder, os);
}


void TWdebug::PrintNLPconstraintsX(TWfolder */*twfolder*/, TransWorhpProblem *tw, ostream *os, int phase) {

	if (tw->solver->transworhp_type == TransWORHP_type::fullDiscretization) { //normales TW

		const int nz = tw->solver->twdiscretization->stuetzstellen(tw->n_dis);

		stringstream a, b;

		for (int i=0; i < tw->solver->n_var; i++) {

			a.str("");
			b.str("");

			// Phase ausgeben
			if (phase > -1) {
				b << "Ph." << setw(2) << phase << " - ";
			}

			a << "X(" << setw(3) << i + tw->solver->Delta1 << ")" ;

			if (i >= nz*(tw->n_ctrl+tw->n_ode)) {
				b << "p(" << setw(3) << i-nz*(tw->n_ctrl+tw->n_ode) << ")" << "          ";
			} else {

				int x = i%(tw->n_ctrl+tw->n_ode);
				int y = i/(tw->n_ctrl+tw->n_ode);
				string add = "  ";
				if (tw->solver->twdiscretization->type==TWdiscretizationType::HermiteSimpson) {
					if (y%2==0) add ="  ";
					else add=".5";
					y/=(tw->solver->twdiscretization->innenpunkt+1);
				}
				else if (tw->solver->twdiscretization->type==TWdiscretizationType::Lobatto) {
					if (y%3==0) add ="  ";
					else if (y%3==1) add=".3";
					else add=".7";
					y/=(tw->solver->twdiscretization->innenpunkt+1);
				}

				if (x<tw->n_ode) {
					b << "x(" << setw(3) << y << add << "," << setw(3) << x << ")    ";
				} else {
					b << "u(" << setw(3) << y << add << "," << setw(3) << x-tw->n_ode << ")    ";
				}
			}

			bool ret = inequality(b, tw->solver->X_low[i], tw->solver->X[i], tw->solver->X_upp[i]);

			if (os) {
				*os << setw(20) << a.str() << b.str() << endl;
			} else {
				MyStatus(a.str(), b.str(), ret ? Status::GREEN : Status::RED);
			}
		}

	}

	else if (tw->solver->transworhp_type == TransWORHP_type::multipleShooting) { // multiple shooting

		stringstream a, b;

		// Zaehler fuer X
		int m = 0;

		ExplTransWorhp* solver = dynamic_cast<ExplTransWorhp*>(tw->solver.get());
		if (!solver) {
			return;
		}

		for (int i = 0; i < tw->n_dis; i++) {

			//Pruefen, ob Gitterpunkt Multiknoten ist
			if (solver->isMultinode(i)) {
				for (int k = 0; k < tw->n_ode; k++) {

					a.str("");
					b.str("");

					// Phase ausgeben
					if (phase > -1) {
						b << "Ph." << setw(2) << phase << " - ";
					}

					a << "X(" << setw(3) << m + solver->Delta1 << ")" ;
					b << "x(" << setw(3) << i <<  "," << setw(3) << k << ")    ";

					bool ret = inequality(b, solver->X_low[m], solver->X[m], solver->X_upp[m]);

					if (os) {
						*os << setw(20) << a.str() << b.str() << endl;
					} else {
						MyStatus(a.str(), b.str(), ret ? Status::GREEN : Status::RED);
					}
					m++;
				}
			}
			for (int k = 0; k < tw->n_ctrl; k++) {

				a.str("");
				b.str("");

				// Phase ausgeben
				if (phase > -1) {
					b << "Ph." << setw(2) << phase << " - ";
				}

				a << "X(" << setw(3) << m + solver->Delta1 << ")" ;
				b << "u(" << setw(3) << i <<  "," << setw(3) << k << ")    ";

				bool ret = inequality(b, solver->X_low[m], solver->X[m], solver->X_upp[m]);

				if (os) {
					*os << setw(20) << a.str() << b.str() << endl;
				} else {
					MyStatus(a.str(), b.str(), ret ? Status::GREEN : Status::RED);
				}
				m++;
			}

		}

		for (int k = 0; k < tw->n_param; k++) {

			a.str("");
			b.str("");

			// Phase ausgeben
			if (phase > -1) {
				b << "Ph." << setw(2) << phase << " - ";
			}

			a << "X(" << setw(3) << m + solver->Delta1 << ")" ;
			b << "p(" << setw(3) << k << ")        ";

			bool ret = inequality(b, solver->X_low[m], solver->X[m], solver->X_upp[m]);

			if (os) {
				*os << setw(20) << a.str() << b.str() << endl;
			} else {
				MyStatus(a.str(), b.str(), ret ? Status::GREEN : Status::RED);
			}
			m++;
		}
	}

	else if (tw->solver->transworhp_type == TransWORHP_type::pseudospectral) {

		const int nz = tw->n_dis;

		stringstream a, b;

		for (int i = 0; i < tw->solver->n_var; i++) {

			a.str("");
			b.str("");

			// Phase ausgeben
			if (phase > -1) {
				b << "Ph." << setw(2) << phase << " - ";
			}

			a << "X(" << setw(3) << i + tw->solver->Delta1 << ")" ;

			if (i >= nz*(tw->n_ctrl+tw->n_ode)) {
				b << "p(" << setw(3) << i-nz*(tw->n_ctrl+tw->n_ode) << ")" << "        ";
			} else {

				int x = i%(tw->n_ctrl+tw->n_ode);
				int y = i/(tw->n_ctrl+tw->n_ode);


				if (x<tw->n_ode) {
					b << "x(" << setw(3) << y << "," << setw(3) << x << ")    ";
				} else {
					b << "u(" << setw(3) << y << "," << setw(3) << x-tw->n_ode << ")    ";
				}
			}

			bool ret = inequality(b,tw->solver->X_low[i], tw->solver->X[i], tw->solver->X_upp[i]);

			if (os) {
				*os << setw(20) << a.str() << b.str() << endl;
			} else {
				MyStatus(a.str(), b.str(), ret ? Status::GREEN : Status::RED);
			}
		}
	}

	else if (tw->solver->transworhp_type == TransWORHP_type::pseudospectral_gauss) {

		// Zaehler fuer X
		int nr = 0;

		stringstream a, b;

		for (int k = 0; k < tw->n_ode; k++) {

			a.str("");
			b.str("");

			// Phase ausgeben
			if (phase > -1) {
				b << "Ph." << setw(2) << phase << " - ";
			}

			a << "X(" << setw(3) << nr + tw->solver->Delta1 << ")" ;
			b << "x(" << setw(3) << 0 <<  "," << setw(3) << k << ")    ";

			bool ret = inequality(b,tw->solver->X_low[nr], tw->solver->X[nr], tw->solver->X_upp[nr]);

			if (os) {
				*os << setw(20) << a.str() << b.str() << endl;
			} else {
				MyStatus(a.str(), b.str(), ret ? Status::GREEN : Status::RED);
			}
			nr++;
		}

		for (int i = 1; i < tw->n_dis-1; i++) {

			for (int k = 0; k < tw->n_ode; k++) {

				a.str("");
				b.str("");

				// Phase ausgeben
				if (phase > -1) {
					b << "Ph." << setw(2) << phase << " - ";
				}

				a << "X(" << setw(3) << nr + tw->solver->Delta1 << ")" ;
				b << "x(" << setw(3) << i <<  "," << setw(3) << k << ")    ";

				bool ret = inequality(b,tw->solver->X_low[nr], tw->solver->X[nr], tw->solver->X_upp[nr]);

				if (os) {
					*os << setw(20) << a.str() << b.str() << endl;
				} else {
					MyStatus(a.str(), b.str(), ret ? Status::GREEN : Status::RED);
				}
				nr++;
			}

			for (int k = 0; k < tw->n_ctrl; k++) {

				a.str("");
				b.str("");

				// Phase ausgeben
				if (phase > -1) {
					b << "Ph." << setw(2) << phase << " - ";
				}

				a << "X(" << setw(3) << nr + tw->solver->Delta1 << ")" ;
				b << "u(" << setw(3) << i <<  "," << setw(3) << k << ")    ";

				bool ret = inequality(b, tw->solver->X_low[nr], tw->solver->X[nr], tw->solver->X_upp[nr]);

				if (os) {
					*os << setw(20) << a.str() << b.str() << endl;
				} else {
					MyStatus(a.str(), b.str(), ret ? Status::GREEN : Status::RED);
				}
				nr++;
			}

		}

		for (int k = 0; k < tw->n_ode; k++) {

			a.str("");
			b.str("");

			// Phase ausgeben
			if (phase > -1) {
				b << "Ph." << setw(2) << phase << " - ";
			}

			a << "X(" << setw(3) << nr + tw->solver->Delta1 << ")" ;
			b << "x(" << setw(3) << tw->n_dis-1 <<  "," << setw(3) << k << ")    ";

			bool ret = inequality(b,tw->solver->X_low[nr], tw->solver->X[nr], tw->solver->X_upp[nr]);

			if (os) {
				*os << setw(20) << a.str() << b.str() << endl;
			} else {
				MyStatus(a.str(), b.str(), ret ? Status::GREEN : Status::RED);
			}
			nr++;
		}

		for (int k = 0; k < tw->n_param; k++) {

			a.str("");
			b.str("");

			// Phase ausgeben
			if (phase > -1) {
				b << "Ph." << setw(2) << phase << " - ";
			}

			a << "X(" << setw(3) << nr + tw->solver->Delta1 << ")" ;
			b << "p(" << setw(3) << k << ")        ";

			bool ret = inequality(b, tw->solver->X_low[nr], tw->solver->X[nr], tw->solver->X_upp[nr]);

			if (os) {
				*os << setw(20) << a.str() << b.str() << endl;
			} else {
				MyStatus(a.str(), b.str(), ret ? Status::GREEN : Status::RED);
			}
		}
	}
}


void TWdebug::PrintNLPconstraintsG(TWfolder */*twfolder*/, TransWorhpProblem *tw, ostream *os, int phase) {

	if (tw->solver->transworhp_type == TransWORHP_type::fullDiscretization) { // normales TW

		const int stufen = tw->solver->twdiscretization->stufen;

		for (int con = 0; con < tw->solver->n_con; con++) {
			stringstream a;
			a << "G(" << setw(3) << con + tw->solver->Delta2 << ")" ;

			stringstream b;

			// Phase ausgeben
			if (phase > -1) {
				b << "Ph." << setw(2) << phase << " - ";
			}

			if (con >= (tw->n_dis-1)*(tw->n_ode)*stufen + tw->n_rand) {

				const int con_neben = con - (tw->n_dis-1)*(tw->n_ode)*stufen - tw->n_rand;

				b << "neben(" << setw(3) << con_neben/tw->n_neben << "," << setw(3) << con_neben%tw->n_neben << ")" << "  ";
			}
			else if (con >= (tw->n_dis-1)*(tw->n_ode)*stufen) {

				const int con_rand = con - (tw->n_dis-1)*(tw->n_ode)*stufen;

				b << "rand(" << setw(3) << con_rand << ")" << "       ";
			}
			else {
				int x = con%(tw->n_ode);
				int y = con/(tw->n_ode);

				string add(" ");

				if (tw->solver->twdiscretization->type==TWdiscretizationType::HermiteSimpson) {

					x = con%(tw->n_ode*2);
					y = con/(tw->n_ode*2);

					if (x % 2 == 0) {
						add="H";
					} else {
						add="S";
					}
					x/=2;
				}
				if (tw->solver->twdiscretization->type==TWdiscretizationType::Lobatto) {

					x = con%(tw->n_ode*3);
					y = con/(tw->n_ode*3);

					if (x % 3 == 0 || x % 3 == 1) {
						add="H";
					} else {
						add="S";
					}
					x/=3;
				}
				//if (con%(n_ctrl+n_ode)<n_ode)
				b << "ode(" << setw(3) << y << "," << setw(3) << x << ")" << " " << add << "  ";
			//else
			//	b << "u(" << setw(3) << con/(n_ctrl+n_ode) << "," << setw(3) << con%(n_ctrl+n_ode)-n_ode << ")" << "   ";
			}

			const bool ret = inequality(b,tw->solver->G_low[con], tw->solver->G[con], tw->solver->G_upp[con]);

			if (os) {
				*os << setw(20) << a.str() << b.str() << endl;
			} else {
				MyStatus(a.str(), b.str(), ret ? Status::GREEN : Status::RED);
			}
		}

	}
	else if (tw->solver->transworhp_type == TransWORHP_type::multipleShooting) { //Mehrzielmethode

		ExplTransWorhp* solver = dynamic_cast<ExplTransWorhp*>(tw->solver.get());
		if (!solver) {
			return;
		}

		for (int i=0; i<solver->n_con; i++) {
			stringstream a;
			stringstream b;

			// Phase ausgeben
			if (phase > -1) {
				b << "Ph." << setw(2) << phase << " - ";
			}

			a << "G(" << setw(3) << i + solver->Delta2 << ")" ;


			// Stetigkeitsbedingungen
			if (i < tw->n_ode*(solver->n_multinodes-1)) {
				int x = i%(tw->n_ode);
				int y = i/(tw->n_ode);

				y = solver->multinodes[y+1];
				b << "cont(" << setw(3) << y << "," << setw(3) << x << ")" << " ";
			}
			// Randbedingungen
			else if (i < tw->n_rand + tw->n_ode*(solver->n_multinodes-1)) {
				b << "rand          ";
			}
			// Nebenbedingungen
			else if (i < tw->n_neben*tw->n_dis + tw->n_rand + tw->n_ode*(solver->n_multinodes-1)) {
				b << "neben         ";
			}
			//Box-NB
			else if (i < ((ExplTransWorhp*)tw)->n_boxNeben*(tw->n_dis-solver->n_multinodes) + tw->n_rand + tw->n_neben*tw->n_dis + tw->n_ode*(solver->n_multinodes-1)) {
				b << "boxNB         ";
			}

			bool ret = inequality(b, tw->solver->G_low[i], tw->solver->G[i], tw->solver->G_upp[i]);

			if (os) {
				*os << setw(20) << a.str() << b.str() << endl;
			} else {
				MyStatus(a.str(), b.str(), ret ? Status::GREEN : Status::RED);
			}
		}

	}
	else if (tw->solver->transworhp_type == TransWORHP_type::pseudospectral) {


		for (int i = 0; i < tw->solver->n_con; i++) {
			stringstream a;
			a << "G(" << setw(3) << i + tw->solver->Delta2 << ")" ;

			stringstream b;

			// Phase ausgeben
			if (phase > -1) {
				b << "Ph." << setw(2) << phase << " - ";
			}

			if (i >= tw->n_dis*tw->n_ode + tw->n_rand) { //Nebenbedingungen

				const int ii = i - (tw->n_dis*tw->n_ode + tw->n_rand);

				b << "neben(" << setw(3) << ii/tw->n_neben << "," << setw(3) << i%tw->n_neben << ")";
			} else if (i >= tw->n_dis*tw->n_ode) { //Randbedingung

				const int r = (i - tw->n_dis*tw->n_ode)%tw->n_rand;

				b << "rand(" << setw(3) << r << ")" << "     ";
			} else { //ODE-Bedingungen / Kollokationsbedingungen

				const int dis = i/tw->n_ode;
				const int ode = i%tw->n_ode;

				b << "col(" << setw(3) << dis << "," << setw(3) << ode << ")  ";
			}

			bool ret = inequality(b, tw->solver->G_low[i], tw->solver->G[i], tw->solver->G_upp[i]);

			if (os) {
				*os << setw(20) << a.str() << b.str() << endl;
			} else {
				MyStatus(a.str(), b.str(), ret ? Status::GREEN : Status::RED);
			}
		}

	}
	else if (tw->solver->transworhp_type == TransWORHP_type::pseudospectral_gauss) {


		for (int i = 0; i < tw->solver->n_con; i++) {

			stringstream a;
			a << "G(" << setw(3) << i + tw->solver->Delta2 << ")" ;

			stringstream b;

			// Phase ausgeben
			if (phase > -1) {
				b << "Ph." << setw(2) << phase << " - ";
			}

			if (i>=(tw->n_dis-2)*(tw->n_ode) + tw->n_ode + tw->n_rand) { //Nebenbedingungen

				const int ii = i - ((tw->n_dis-2)*tw->n_ode + tw->n_ode + tw->n_rand);

				b << "neben(" << setw(3) << ii/tw->n_neben+1 << "," << setw(3) << i%tw->n_neben << ")";
			} else if (i>=(tw->n_dis-2)*(tw->n_ode) + tw->n_ode) { //Randbedingung

				const int r = (i - ((tw->n_dis-2)*tw->n_ode + tw->n_ode))%tw->n_rand;

				b << "rand(" << setw(3) << r << ")" << "     ";

			} else if (i>=(tw->n_dis-2)*tw->n_ode) { // Endbedingung

				const int ode = (i-(tw->n_dis-2)*tw->n_ode)%tw->n_ode;

				b << "end(" << setw(3) << ode << ")      ";
			} else { //ODE-Bedingungen / Kollokationsbedingungen

				const int dis = i/tw->n_ode + 1;
				const int ode = i%tw->n_ode;

				b << "col(" << setw(3) << dis << "," << setw(3) << ode << ")  ";
			}

			bool ret = inequality(b, tw->solver->G_low[i], tw->solver->G[i], tw->solver->G_upp[i]);

			if (os) {
				*os << setw(20) << a.str() << b.str() << endl;
			} else {
				MyStatus(a.str(), b.str(), ret ? Status::GREEN : Status::RED);
			}
		}

	}

}




void TWdebug::PrintNLPconstraintsG(TWfolder *twfolder, ostream *os) {

	for (int i = 0; i < twfolder->n_con; i++) {
		stringstream a;
		a << "G(" << setw(3) << i + twfolder->G_offset << ")" ;

		stringstream b;
		b << "con(" << setw(3) << i << ")        ";

		const bool ret = inequality(b, twfolder->G_low[i], twfolder->G[i], twfolder->G_upp[i]);

		if (os) {
			*os << setw(20) << a.str() << b.str() << endl;
		} else {
			MyStatus(a.str(), b.str(), ret ? Status::GREEN : Status::RED);
		}
	}
}

//TODO welche Werte im WORHP G gehoehren zu welchen Beschraenkungen in TW
int TWdebug::CheckNLPconstraints(TWfolder* twfolder, ostream* os) {

	string s("Check NLP Constraints");
	if (os) {
		*os << "# " << s << endl;
	} else {
		MyStatus("CHECK", s, Status::NORMAL);
	}

	for (size_t i = 0; i < twfolder->phases.size(); i++) {

		TransWorhpProblem *tw = twfolder->phases[i];

		//nur fuer explTW
		ExplTransWorhp* solver = dynamic_cast<ExplTransWorhp*>(tw->solver.get());
		if (!solver) {
			return true;
		}

		for (int k = 0; k < solver->n_con; k++) {

			//cout << tw->G_low[k] << " " << tw->G_upp[k] << endl;

			if (solver->G_low[k] == -solver->Infty && solver->G_upp[k] == solver->Infty) {

				stringstream ss;



				if (k < tw->n_ode*(solver->n_multinodes-1)) {
					ss << "Stetigkeitsbedingung ";

				}
				else if (k >= tw->n_ode*(solver->n_multinodes-1) && k < tw->n_ode*(solver->n_multinodes-1)+tw->n_rand) {
					ss << "Randbedingung ";
				}
				else if (k >= tw->n_ode*(solver->n_multinodes-1)+tw->n_rand && k < tw->n_ode*(solver->n_multinodes-1)+tw->n_rand+tw->n_neben * solver->twdiscretization->stuetzstellen(tw->n_dis)) {
					ss << "Nebenbedingung ";
				}
				else {
					ss << "Box-Nebenbedingung ";
				}



				ss << k << " in Phase " << i << " ist unbeschraenkt!";



				if (os) {
					*os << "# " << ss.str() << endl;
				} else {
					MyStatus("CHECK", ss.str(), Status::ERR);
				}

				//return false;

			}
		}

	}

	return true;
}



bool TWdebug::inequality(ostream &b, double low, double val, double upp) {

	b << setw(10) << low << " <= " <<
		  setw(13) << val << " <= " <<
		  setw(10) << upp;

	if (low <= val+1e-6 && val-1e-6 <= upp) {
		return true;
	} else {

		double m = 0;
		if (val < low) m = val-low; // -
		if (val > upp) m = val-upp; // +

		b << "  <========= " << setw(10) << m;
		return false;
	}
}


Status TWdebug::matrixentry(ostream &b, double val, double ref) {

	b << setw(20) << val << setw(20) << ref;
	const double aux = fabs(val - ref);
	if (aux > 1e-5) {
		b << "  <========= " << setw(10) << aux;
		return Status::RED;
	}
	return Status::GREEN;
}

}
