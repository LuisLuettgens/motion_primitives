#include "TWfolder.h"

#include "TWproblem.h"
#include "TransWORHP.h"
#include "FullDisTransWORHP.h"
#include "TWstrings.h"
#include "TWdebug.h"
#include "twstatus.h"

#include "../base/defines.h"

#include "conversion.h"

#include "transworhp_version.h"

#include "worhp_info.h"

#include <cstring>
#include <algorithm>

#ifdef TW_WITH_BOOST
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;
#endif

#ifdef WIN32
#include <windows.h>
#endif

namespace tw {

int TWDEBUG = 0;

TWfolder::TWfolder(TWparameter *twparam, int _con)
	:	meshRefFlag(false), initFlag(false),
		twparameter(twparam), verbose(true), Interrupt(false), viewer(nullptr),
		G(nullptr), G_low(nullptr), G_upp(nullptr), Mu(nullptr),
		G_offset(0),
		n_con(_con),
		tmpG1(n_con), tmpG2(n_con),
		showOnly(false) {

	/* Check version of library and header files */
	//#ifndef WIN32
	//	  CHECK_WORHP_VERSION
	//#endif
}


TWfolder::~TWfolder() {

	if (worhp_o.initialised) {
		WorhpFree(&worhp_o, &worhp_w, &worhp_p, &worhp_c);
	}
}


void TWfolder::Add(TransWorhpProblem *tw) {

	if (tw->solver == nullptr) {
		MyStatus("Check", "TWfolder::Add: you have to set a solver. Call ph.setSolver(..)) befor "
		                  "twfolder.Add(..)!",
		         Status::ERR);
		exit(1);
	}

	phases.push_back(tw);
}


int TWfolder::Init(bool verbose_, Params *params) {

	if (phases.empty()) {
		MyStatus("Check", "Please add a phase to the TWfolder!", Status::ERR);
		return 1;
	}

	initFlag = true;
	verbose = verbose_;

	if (verbose) {
		MyStatus("TransWORHP", "-------------------------------------------------------------",
		         Status::NORMAL);
		MyStatus("TransWORHP", " This is TransWORHP " + std::string(TransWORHP_VERSION) +
		                           ", the OCP solver for WORHP, ",
		         Status::NORMAL);
		MyStatus("TransWORHP", "              developed by University of Bremen.", Status::NORMAL);
		MyStatus("TransWORHP", "-------------------------------------------------------------",
		         Status::NORMAL);

		///infoMagic();

		//twparameter->ParseXML(xml);

		MyStatus(
		    "TransWORHP",
		    TWstrings::string_discretization[static_cast<int>(twparameter->twdiscretization.type)],
		    Status::WARN);
	}

	/** Start PHASE */
	{
		// Anzahl Variablen und Nebenbed. zaehlen
		int nn = 0, mm = n_con;

		for (TransWorhpProblem* ph : phases) {

			ph->solver->twdiscretization = &twparameter->twdiscretization;
			ph->solver->init0(); // Struktur-Matrizen initialisieren

			nn += ph->solver->n_var;
			mm += ph->solver->n_con;
		}

		//cout << "init0   " << nn << " " << mm << endl;
		init0(nn, mm); // Constraint-Matrix initialisieren

	}
	/** End PHASE */

	int nn = 0, mm = 0, nzen = 0;
	try {
		 for (TransWorhpProblem* ph : phases) {

			ph->solver->GetSpace(nn, mm);
			nn += ph->solver->n_var;
			mm += ph->solver->n_con;
			nzen += ph->n_zen;
		}

		GetSpace(mm);
		mm += n_con;
	}
	catch (std::string &s) {
		MyStatus("TransWORHP", "exception: GetSpace() ", Status::ERR);

#ifdef TRANSWORHP_GRAPHICS
		exit(1);
#else
		/*worhp_o = 0;
		worhp_w = 0;
		worhp_p = 0;
		worhp_c = 0;*/
		return 100;
#endif
	}

	WorhpPreInit(&worhp_o, &worhp_w, &worhp_p, &worhp_c);

	/* Uncomment this to get more info on data structures */
	// WorhpDiag(&o, &w, &p, &c);


	// Read Worhp XML parameter file
	int status;

	char buf[100];
	strcpy(buf, twparameter->paramfile_worhp.c_str());
	//#ifdef USE_ZEN
	//::ReadParamsNoInit(&status, buf, &p);
	//#else

	if (params) {
		worhp_p = *params;
	} else {
		::ReadParams(&status, buf, &worhp_p);
	}
	//#endif


	//worhp_p.UseZen = true;

	if (!verbose) {
		worhp_p.NLPprint = 0;
	}

	status = OK; // ignore errors: WorhpReadParams falls back to standard values.
	worhp_p.MatrixCC = false;

	if (twparameter->USERDG == -1) worhp_p.UserDG = false;
	if (twparameter->USERDG == +1) worhp_p.UserDG = true;

	if (twparameter->USERDF == -1) worhp_p.UserDF = false;
	if (twparameter->USERDF == +1) worhp_p.UserDF = true;

	if (twparameter->USERHM == -1) worhp_p.UserHM = false;
	if (twparameter->USERHM == +1) worhp_p.UserHM = true;

	/* Specify number of variables and constraints. */
	worhp_o.n = nn;
	worhp_o.m = mm;

	worhp_o.k = nzen;

	/* Specify nonzeros of derivative matrices. */
	// w.DF.nnz = WorhpMatrix_Init_Dense; /* dense, structure init by solver */
	// w.DG.nnz = WorhpMatrix_Init_Dense; /* dense, structure init by solver */
	// w.HM.nnz = o.n;                    /* diagonal, structure init by user */

	try {
		worhp_w.DF.nnz = 0;
		worhp_w.DG.nnz = 0;
		worhp_w.HM.nnz = 0;

		if (verbose) {
			MyStatus("Hessian", TWstrings::string_hessianstructure[twparameter->hessianstructure],
			         Status::WARN);
			MyStatus("Hessian", TWstrings::string_hessianvalues[twparameter->hessianvalues],
			         Status::WARN);
		}

		for (TransWorhpProblem* ph : phases) {

			int DF_nnz, DG_nnz, HM_nnz;
			// Start PHASE: Indirekter Aufruf von TWfolder::DG_structure
			// (spaltenweise)
			ph->solver->Structure_Sizes(this, twparameter->hessianstructure, DF_nnz, DG_nnz, HM_nnz);

			worhp_w.DF.nnz += DF_nnz;
			worhp_w.DG.nnz += DG_nnz;
			worhp_w.HM.nnz += HM_nnz;

			//cout<< "LOGS " << DF_nnz << " " << DG_nnz << " " << HM_nnz << endl;
			//cout<< "LOGS " << w.DF.nnz << " " << w.DG.nnz << " " << w.HM.nnz << endl;
		}
	}
	catch (std::string &s) {
#ifdef WIN32
#ifdef _MSC_VER
		MessageBox(0, (LPCWSTR)s.c_str(), (LPCWSTR)L"TransWORHP", MB_ICONEXCLAMATION);
#else
		MyStatus("TransWORHP", "exception: " + s, Status::ERR);
#endif
#else
		MyStatus("TransWORHP", "exception: " + s, Status::ERR);
#endif

#ifdef TRANSWORHP_GRAPHICS
		exit(1);
#else
		return 102;
#endif
	}

	//cout << "WORHP INIT" << endl;
	/* Data structure initialisation. */
	WorhpInit(&worhp_o, &worhp_w, &worhp_p, &worhp_c);

	// MyStatus("?","WORHP_C.status " + ToString(worhp_c.status) + " " + ToString(FirstCall) , Status::WARN);
	if (worhp_c.status != FirstCall) {

#ifndef TRANSWORHP_GRAPHICS
		//MyStatus("TransWORHP", "Initialisation failed.", Status::ERR);
		return 101;
#else
#ifdef WIN32
#ifdef _MSC_VER
		MessageBox(0,(LPCWSTR)L"WORHP Initialisation failed.",(LPCWSTR)L"TransWORHP",MB_ICONEXCLAMATION);
#else
		MyStatus("TransWORHP", "Main: Initialisation failed.", Status::ERR);
#endif
#else
		MyStatus("TransWORHP", "Main: Initialisation failed.", Status::ERR);
#endif
		exit(1);
#endif
	}

	/* Uncomment this to have a look at the workspace slices */
	//PrintIWMT(&w);
	//PrintRWMT(&w);

	//cout << "CONNECT " << endl;
	 for (TransWorhpProblem* ph : phases) {
		//cout << "--CONNECT " << endl;
		ph->solver->Connect(worhp_o, worhp_p);
	}

	/** Start PHASE */
	Connect();
	/** End PHASE */

	for (TransWorhpProblem* ph : phases) {
		ph->solver->Boundary();
	}

	/** Start PHASE */
	Boundary();
	/** End PHASE */

	//TWdebug::PrintNLPconstraints(this);

	/*
	* Specify matrix structures in CC format.
	* Indexing conventions:
	* Use native 0-based indices to access the workspace slices obtained
	* from the IWS_PTR and RWS_PTR macros.
	* Use Fortran 1-based indices for specifying the matrix values or
	* structure in the respective val, row and col workspace slices.
	*/

	/* Only specify the HM structure when not using BFGS! */

	/*if (worhp_w.DF.NeedStructure)*/ {
		// cout << " DF Structure: " <<  DF_structure(&w.DF) << endl;
		int offset = 0;

		 for (TransWorhpProblem* ph : phases) {
			offset = ph->solver->DF_structure(&worhp_w.DF, offset);
			//cout << "OFF " << offset << endl;
		}
	}

	if (worhp_w.DG.NeedStructure) {
		// cout << " DG Structure: " <<  DG_structure(&w.DG) << endl;
		int offset = 0;
		int auxCon = 0;

		 for (TransWorhpProblem* ph : phases) {
			// Indirekter Aufruf von TWfolder::DG_structure
			offset = ph->solver->DG_structure(this, &worhp_w.DG, offset);

			phasenOffset.push_back(auxCon);
			auxCon += ph->solver->n_con;
		}
		phasenOffset.push_back(auxCon);
	}

	if (worhp_w.HM.NeedStructure) {
		HM_structure(worhp_w.DF, worhp_w.DG, worhp_w.HM);
	}

	return 0;
}


#ifdef TRANSWORHP_GRAPHICS
int TWfolder::Init(Viewer *view) {
#else
int TWfolder::Init(Viewer*) {
#endif

	// falls TW noch nicht mit Init() initialisiert wurde
	if (!initFlag) {
		Init();
	}

#ifdef TRANSWORHP_GRAPHICS

	bool old = false;

	if (viewer) {
		viewer->closeAll();
		old = true;
	}

	viewer = view;

	if (viewer) {
		viewer->init(this);
	}

	if (/*n_dis<100 && */ viewer) {
		if (twparameter->showDF) viewer->Matrix("DF", &worhp_w.DF);
		if (twparameter->showDG) viewer->Matrix("DG", &worhp_w.DG);
		if (twparameter->showHM) viewer->Matrix("HM", &worhp_w.HM);

		// Plots fuer MeshRef
		 for (TransWorhpProblem* ph : phases) {
			if (meshRefFlag) {
				FullDisTransWorhp* tw = dynamic_cast<FullDisTransWorhp*>(ph->solver.get());
				if (tw) {
				      viewer->disFehler(tw->SCHRITTWEITE, tw->FEHLER, tw->twdiscretization);
				}
			}
			viewer->restrictionPlot(ph->solver->T.data(),ph->solver->Lambda,ph->n_dis,ph->n_ctrl,ph->n_ode,ph->solver->twdiscretization);
		}
	}

	if (old) viewer->tilePlots(1);

#endif
	return 0;
}


std::string TWfolder::type_G(int row) const {
	std::stringstream a;
	a << " con  " << std::setw(2) << row%n_con << " ";
	return a.str();
}


int TWfolder::meshRef() {

	// verhindern, dass meshRef bei explTW und pmTW ausgefuehrt wird
	// dann Loop aufrufen
	auto ph = std::find_if(phases.begin(),phases.end(),[](TransWorhpProblem *tw){return tw->solver->transworhp_type!=TransWORHP_type::fullDiscretization;});
	if (ph != phases.end()) {
		MyStatus("Check", "using Loop insted of meshRef()", Status::ERR);
		Loop();
		return 1;
	}

	meshRefFlag = true;

	// neuer Zeitvektor in Phase
	std::vector<std::vector<double> > T_neu;
	T_neu.reserve(phases.size());

	// Anzahl der neunen Punkte in Phase
	std::vector<int> ndis_neu;
	ndis_neu.reserve(phases.size());

	// fuer Startschaetzung
	std::vector<exportTW> eTW;
	eTW.reserve(phases.size());

	// Maximum pro Phase
	std::vector<double> max;
	max.reserve(phases.size());

	// erstes Loop
	Loop(2, 1);

	for (size_t i = 0; i < phases.size(); i++) {

		max.push_back(1.0);

		FullDisTransWorhp* tw = dynamic_cast<FullDisTransWorhp*>(phases[i]->solver.get());

		// als Startschaetzung
		eTW.push_back(tw->outTW());

		T_neu.push_back(tw->refineAlg(max.at(i)));

		ndis_neu.push_back(T_neu.at(i).size());
	}

	// Genauigkeit (aus Param-XML)
	const double tol = twparameter->meshref_tol;

	// Maximum ueber alle Phasen
	double MAX = 1.0;

	// Zaehler in Schleife (Abbruch)
	int k = 0;

	while (MAX > tol) {
		//cout << "STARTE SCHLEIFE" << endl;

		// Abbruch bei zu vielen Gitter-Schritten
		if (k == twparameter->meshref_maxIter) {
			MyStatus("meshRef",
			         "too many iterations. maxIter=" + std::to_string(twparameter->meshref_maxIter),
			         Status::ERR);
			break;
		}

		//// Wechsel der Diskretisierung
		// if (k%2 == 0) {
		//	twparameter->twdiscretization = TWdiscretization(TWdiscretizationType::HermiteSimpson,2,1);
		//} else {
		//	twparameter->twdiscretization = TWdiscretization(TWdiscretizationType::Trapez,1,0);
		//}

		if (viewer) viewer->closeAll();

		for (size_t i = 0; i < phases.size(); i++) {

			FullDisTransWorhp* tw = dynamic_cast<FullDisTransWorhp*>(phases[i]->solver.get());

			tw->n_dis = ndis_neu.at(i);

			tw->newGrid(T_neu.at(i));

			tw->lagrange_integral.assign(tw->n_integral*tw->n_dis*2,0.0);
		}

		Reinit(1);

		if (/*twparameter->meshref_VER &&*/ viewer) {
			Init(viewer);

			for (size_t i = 0; i < phases.size(); i++) {
				FullDisTransWorhp* tw = dynamic_cast<FullDisTransWorhp*>(phases[i]->solver.get());
				tw->inTW(eTW.at(i));
			}

			Loop(2, 1); // hier darf nicht Loop(0,..) stehen, da sonst der viewer nicht nachkommt
		}
		else {
			for (size_t i = 0; i < phases.size(); i++) {
				FullDisTransWorhp* tw = dynamic_cast<FullDisTransWorhp*>(phases[i]->solver.get());
				tw->inTW(eTW.at(i));
			}
			Loop(0, 1);
		}

		for (size_t i = 0; i < phases.size(); i++) {
			FullDisTransWorhp* tw = dynamic_cast<FullDisTransWorhp*>(phases[i]->solver.get());
			eTW.at(i) = tw->outTW();
			T_neu.at(i) = tw->refineAlg(max.at(i));
			ndis_neu.at(i) = T_neu.at(i).size();
		}

		MAX = *std::max_element(max.begin(), max.end());

		//cout << "maximaler Fehler: " << MAX << endl;

		k++;
	}

	if (viewer) viewer->closeAll();

	for (size_t i = 0; i < phases.size(); i++) {

		TransWorhp *tw = phases[i]->solver.get();

		tw->n_dis = ndis_neu.at(i);

		tw->newGrid(T_neu.at(i));
	}

	Reinit(1);
	Init(viewer);

	for (size_t i = 0; i < phases.size(); i++) {
		FullDisTransWorhp* tw = dynamic_cast<FullDisTransWorhp*>(phases[i]->solver.get());
		tw->inTW(eTW.at(i));
	}

	Loop(0, 0);

	//cout << "nach letztem loop" << endl;

	return 0;
}


#ifndef TW_WITH_BOOST
bool TWfolder::verify(bool) const {
	
	MyStatus("verify", "please activate boost.", Status::NORMAL);
	return true;
}
#else
bool TWfolder::verify(bool verbose) const {
	
	auto ph = phases.front();

	std::vector<double> x(ph->n_ode);
	for (int i = 0; i < ph->n_ode; i++) {
		x[i] = ph->x(0,i);
	}

	size_t steps = integrate_const(make_controlled<runge_kutta_fehlberg78<std::vector<double>>>(twparameter->abserr, twparameter->relerr),
		[&ph]( const std::vector<double> &x , std::vector<double> &dxdt , double t ) {

			std::vector<double> tmp_u(ph->n_ctrl);
			ph->solver->GetControl(tmp_u.data(),t);

			ph->ode(dxdt.data(),t,x.data(),tmp_u.data(),&ph->solver->X[ph->p_index(0)]);
			
		}, x, ph->solver->T.front(), ph->solver->T.back(), twparameter->stepsize/*,
		[](const std::vector<double> &x , double t) {
			std::cout << t << std::endl;
		}*/);

	if (verbose) {
		MyStatus("verify", "integration steps: " + std::to_string(steps), Status::NORMAL);
		
		MyStatus("verify", "states: ", Status::NORMAL);
		
		for (int i = 0; i < ph->n_ode; i++) {
			MyStatus("|", "state_" + std::to_string(i) + ": " + std::to_string(x[i]), Status::NORMAL);
		}
		
		MyStatus("verify", "error: ", Status::NORMAL);
	}

	double max_err = 0.0;
	for (int i = 0; i < ph->n_ode; i++) {
		const double err = std::abs(x[i] - ph->x(ph->n_dis-1,i));
		if (err > max_err) {
			max_err = err;
		}
		if (verbose) {
			MyStatus("|",  "error state_" + std::to_string(i) + ": " + std::to_string(err), Status::NORMAL);
		}
	}

	if (max_err > twparameter->abserr) {
		MyStatus("verify", "max. error = " + std::to_string(max_err) + ", (" + std::to_string(twparameter->abserr) + ")", Status::ERR);
		return false;
	} else {
		MyStatus("verify", "ok - max. error = " + std::to_string(max_err), Status::NORMAL);
		return true;
	}
}
#endif


#ifdef TRANSWORHP_GRAPHICS
int TWfolder::Loop(int wait, int terminate) {
#else
int TWfolder::Loop(int,int) {
#endif

	//cout << SDL_GetTicks() << "*** TransWORHP::Loop..." << endl;

	if (!TWdebug::CheckNLPconstraints(this)) {
		return 1;
	}

	if (viewer) {
#ifdef TRANSWORHP_GRAPHICS
		// Alter Thread laeuft noch?
		if (!thethread0->Active()) {

			viewer->startTime = SDL_GetTicks();
			viewer->running = false;

			viewer->update();
			//cout << "LOOP" << wait << terminate << endl;
			viewer->Loop(wait, terminate);
		}
#endif
	}
	else WorhpLoop();

	//cout << "Loop ende " << endl;
	return 0;
}


void TWfolder::Show() {
	showOnly = true;

	MyStatus("Show-mode", "show loaded results", Status::WARN);
	MyStatus("Show-mode", "no optimization is performed", Status::WARN);
	Loop();
	MyStatus("Show-mode", "finished", Status::WARN);
}


double TWfolder::obj() {

	double ret = 0.0;

	for (TransWorhpProblem* ph : phases) {
		ret += ph->solver->Objective(worhp_w.ScaleObj);
	}

	return ret;
}


void TWfolder::g_boundary(double*, double*) {}


void TWfolder::con(double*) {}


bool TWfolder::con_structure(DiffStructure&) {
	return false;
}


bool TWfolder::con_diff(DiffStructure&, int) {
	return false;
}


int TWfolder::WorhpLoop() {

	if (showOnly) return 1;

	worhp_timing.Start(TIME_ALL);

	//DEBUG = 1;
	/* Reverse Communication loop */
	while (worhp_c.status < TerminateSuccess &&  worhp_c.status > TerminateError) {

		if (GetUserAction(&worhp_c, callWorhp)) {
			worhp_timing.Start(TIME_WORHP);
			Worhp(&worhp_o, &worhp_w, &worhp_p, &worhp_c);
			/* DO NOT call DoneUserAction here! */
			worhp_timing.Stop(TIME_WORHP);

			if (TWDEBUG) std::cout << "callW" << std::endl;
		}

		if (TWDEBUG)  {
			std::cout << "Loop" << std::endl;
			for (int i = 0; i < worhp_o.n; i++) {

				std::cout << std::setw(20) << worhp_o.Lambda[i];
				if ((i % 5) == 4) std::cout << std::endl;
			}
			std::cout << std::endl;
			for (int i = 0; i < worhp_o.m; i++) {
				std::cout << std::setw(20) << worhp_o.Mu[i];
				if ((i % 5) == 4) std::cout << std::endl;
			}
			std::cout << std::endl;
		}
		if (GetUserAction(&worhp_c, iterOutput)) {
			worhp_timing.Start(TIME_OUTPUT);
			IterationOutput(&worhp_o, &worhp_w, &worhp_p, &worhp_c);

#ifdef TRANSWORHP_GRAPHICS
			if (viewer) viewer->update();
#endif

			if (!step()) break;
			if (Interrupt) break;

			if (twparameter->twconsole.Active() && worhp_w.MajorIter < 2)
				twparameter->twconsole.Resize(/*twparameter->twconsole.fw,twparameter->twconsole.fh,twparameter->twconsole.bw,*/50);

			DoneUserAction(&worhp_c, iterOutput);
			worhp_timing.Stop(TIME_OUTPUT);

			if (TWDEBUG) std::cout << "iterout" << std::endl;
		}

		if (GetUserAction(&worhp_c, evalF)) {

			worhp_timing.Start(TIME_F);

			worhp_o.F = obj();

			DoneUserAction(&worhp_c, evalF);

			worhp_timing.Stop(TIME_F);

			if (TWDEBUG) {
				std::cout << "evalF" << std::endl;
				std::cout << std::setw(20) << worhp_o.F << std::endl;
			}
		}

		if (GetUserAction(&worhp_c, evalG)) {

			worhp_timing.Start(TIME_G);

			 for (TransWorhpProblem* ph : phases) {
				ph->solver->Constraints();
			}
			/** Start PHASE */
			Constraints();
			/** End PHASE */
			DoneUserAction(&worhp_c, evalG);

			worhp_timing.Stop(TIME_G);

			if (TWDEBUG) {
				std::cout << "evalG" << std::endl;
				for (int i = 0; i < worhp_o.m; i++) {
					std::cout << std::setw(20) << worhp_o.G[i];
					if ((i % 5) == 4) std::cout << std::endl;
				}
				std::cout << std::endl;
			}
		}

		if (GetUserAction(&worhp_c, evalDF)) {

			worhp_timing.Start(TIME_DF);

			 for (TransWorhpProblem* ph : phases) {
				ph->solver->DF_calculate(worhp_w.DF, worhp_w.ScaleObj);
			}
			DoneUserAction(&worhp_c, evalDF);
			worhp_timing.Stop(TIME_DF);

			if (TWDEBUG) std::cout << "evalDF" << std::endl;

			if (TWDEBUG)  {
				std::cout << "Loop" << std::endl;
				for (int i = 0; i < worhp_w.DF.nnz; i++) {
					std::cout << std::setw(20) << worhp_w.DF.val[i];
					if ((i % 5) == 4) std::cout << std::endl;
				}
				std::cout << std::endl;
			}
		}

		if (GetUserAction(&worhp_c, evalDG)) {

			worhp_timing.Start(TIME_DG);

			 for (TransWorhpProblem* ph : phases) {
				/** Start PHASE */
				// Indirekter Aufruf von TWfolder::DG_calculate
				ph->solver->DG_calculate(this, worhp_w.DG);
				/** End PHASE */
			}
			DoneUserAction(&worhp_c, evalDG);

			worhp_timing.Stop(TIME_DG);

			if (TWDEBUG)  {
				std::cout << "evalDG" << std::endl;

				for (int i = 0; i < worhp_w.DG.nnz; i++) {
					std::cout << std::setw(20) << worhp_w.DG.val[i];
					if ((i % 5) == 4) std::cout << std::endl;
				}
				std::cout << std::endl;
			}
		}

		if (GetUserAction(&worhp_c, evalHM)) {

			worhp_timing.Start(TIME_HM);

			HM_calculate();
			DoneUserAction(&worhp_c, evalHM);

			worhp_timing.Stop(TIME_HM);
		}

		if (GetUserAction(&worhp_c, fidif)) {

			worhp_timing.Start(TIME_FIDIF);

			bool do_dg = false, do_df = false;
			if (GetUserAction(&worhp_c, fidif_DG)) {
				do_dg = true;
			}
			if (GetUserAction(&worhp_c, fidif_DF)) {
				do_df = true;
			}

			WorhpFidif(&worhp_o, &worhp_w, &worhp_p, &worhp_c);

			// Testen der Fidif-Routine von WORHP... kann weg?
			if (TWDEBUG) std::cout << "evalFiDIF" << std::endl;

			if (do_dg && GetUserAction(&worhp_c, callWorhp)) {

			}

			if (do_df && GetUserAction(&worhp_c, callWorhp)) {

			}

			worhp_timing.Stop(TIME_FIDIF);
		}
	}

	worhp_timing.Stop(TIME_ALL);

	StatusMsg(&worhp_o, &worhp_w, &worhp_p, &worhp_c);

	if (verbose) {

		MyStatus("", "", Status::NORMAL);
		MyStatus("TransWORHP", "TransWORHP statistics", Status::NORMAL);
		MyStatus("TransWORHP", "      calls to WORHP", Status::NORMAL);

		std::stringstream ss;
		worhp_timing.PrintTimes(ss);
		MyStatus("TransWORHP", ss.str(), Status::NORMAL);

		MyStatus("TransWORHP", "      calls in OCP                                f     df/d(x,u)   df/dp", Status::NORMAL);

		 for (TransWorhpProblem* ph : phases) {

			 MyStatus("TransWORHP", "                     objective function  " +
			                            ph->solver->twcount_calls.obj.ToString(),
			          Status::NORMAL);
			 MyStatus("TransWORHP", "                     ode system          " +
			                            ph->solver->twcount_calls.ode.ToString(),
			          Status::NORMAL);
			 MyStatus("TransWORHP", "                     boundary constraints" +
			                            ph->solver->twcount_calls.rand.ToString(),
			          Status::NORMAL);
			 MyStatus("TransWORHP", "                     state constraints   " +
			                            ph->solver->twcount_calls.neben.ToString(),
			          Status::NORMAL);
			 MyStatus("TransWORHP", "                     integration         " +
			                            ph->solver->twcount_calls.integrate.ToString(),
			          Status::NORMAL);

			 ph->terminate();
		}
	}

	if (meshRefFlag) {

		std::stringstream a;

		FullDisTransWorhp* tw = dynamic_cast<FullDisTransWorhp*>(phases.at(0)->solver.get());

		a << "Refinment step " << tw->SCHRITTWEITE.size();
		MyStatus("MeshRef", a.str(), Status::NORMAL);

		for (size_t i = 0; i < phases.size(); i++) {

			if (phases.size() > 1) { //Ausgabe, welche Phase
				a.str("");
				a << "Phase: " << i;
				MyStatus("MeshRef", a.str(), Status::NORMAL);
			}

			FullDisTransWorhp* tw = dynamic_cast<FullDisTransWorhp*>(phases[i]->solver.get());

			tw->SCHRITTWEITE.push_back(tw->getGrid());

			tw->diskretisierungsfehler();
		}
	}

	//cout << "WorhpLoop() ende! " << endl;
	return 0;
}


void TWfolder::HM_structure(WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM) {

	// cout << " HM Structure: " <<  HM_structure(w.DF, w.DG, &w.HM) << endl;
	int offset = 0;

	for (TransWorhpProblem* ph : phases) {
		offset = ph->solver->HM_structure_ohne_Diag(twparameter->hessianstructure, DF, DG, &HM, offset);
		//cout << "OFFSET" << offset << endl;
	}

	int off1 = offset;

	for (TransWorhpProblem* ph : phases) {
		ph->solver->HM_start_diag = off1;
		off1 += ph->solver->n_var;
	}

	for (int i = 0; i < worhp_o.n; ++i) {
		HM.row[offset] = i + 1;
		HM.col[offset] = i + 1;
		offset++;
	}
}


int TWfolder::Reinit(bool verbose_) {

	if (worhp_o.initialised) {
		WorhpFree(&worhp_o, &worhp_w, &worhp_p, &worhp_c);
	}

	return Init(verbose_);
}


void TWfolder::HM_calculate() {

	for (TransWorhpProblem* ph : phases) {
		ph->solver->HM_calculate(twparameter->hessianvalues, worhp_w.DF, worhp_w.DG, worhp_w.HM, worhp_w.ScaleObj, worhp_o.Mu);
	}
}

/** Start PHASE */

void TWfolder::Connect() {

	G = &worhp_o.G[G_offset];
	G_low = &worhp_o.GL[G_offset];
	G_upp = &worhp_o.GU[G_offset];

	Mu = &worhp_o.Mu[G_offset];

	for (int i = 0; i < n_con; i++) {
		Mu[i] = 0.1;
	}
}


void TWfolder::Boundary() {

	for (int i = 0; i < n_con; i++) {
		G_low[i] = 0.0;
		G_upp[i] = 0.0;
	}
	g_boundary(G_low, G_upp);

}


// wird vom DG_structure der Phase aufgerufen
void TWfolder::DG_structure(int &ind, int colindex, WorhpMatrix *DG) const {

	// fuer alle Beschraenkungen
	for (int c = 0; c < n_con; c++) {

		if (DS_con.check(c, colindex - 1)) {
			if (DG) {
				DG->col[ind] = colindex;
				DG->row[ind] = c + G_offset + 1;
				DG->val[ind] = 0;
			}
			ind++;
		}
	}
}


void TWfolder::DG_calculate(int DG_start, int &ind, int colindex, WorhpMatrix &DG) {

	if (n_con > 0) {

		if (con_diff(DS_con, colindex - 1)) {

			// fuer alle Beschraenkungen
			for (int l = 0; l < n_con; l++) {

				if (DS_con.check(l, colindex - 1)) {

					DG.val[DG_start + ind] = DS_con.get(l, colindex - 1);
					ind++;
				}

			}
			return;
		}

		const double eps = 1e-6;

		//cout << "DG_calc colindex " << colindex-1 << endl;

		double v = worhp_o.X[colindex - 1];
		worhp_o.X[colindex - 1] += eps;
		Constraints2(tmpG1.data());
		worhp_o.X[colindex - 1] = v - eps;
		Constraints2(tmpG2.data());
		worhp_o.X[colindex - 1] = v;

		// fuer alle Beschraenkungen
		for (int c = 0; c < n_con; c++) {

			if (DS_con.check(c, colindex - 1)) {

				DG.val[DG_start + ind] = (tmpG1[c] - tmpG2[c]) / (2 * eps);
				ind++;
			}
		}
	}
}


void TWfolder::init0(int nn, int /*mm*/) {
	DS_con.Init(n_con, nn);
}


void TWfolder::GetSpace(int delta2) {

	G_offset = delta2;

	DS_con.useStructure(con_structure(DS_con));

	DS_con.finish();

	MyStatus("", "", Status::NORMAL);
}


void TWfolder::Constraints() {
	Constraints2(G);
}


void TWfolder::Constraints2(double *GG) {
	con(GG);
}

/** End PHASE */

void TWfolder::WriteSensitivity(std::ostream &os, char var, char pert) {

	int npert = ZenN(&worhp_o, &pert);
	int nvar = ZenN(&worhp_o, &var);

	std::cout << "D" << var << "/D" << pert << ": " << nvar << " x " << npert << std::endl;
        auto d = new double[nvar]; // Das sollte gross genug sein...

        for (int i = 1; i <= npert; i++) {

		ZenGetD(&worhp_o, &worhp_w, &worhp_p, &worhp_c,
			&var, &pert, &i, d);

		os.setf(std::ios::scientific);
		os << std::setprecision(15);
		for (int j = 0; j < nvar; j++) {
			os << " " << std::setw(24) << d[j];
		}

		os << std::endl;
	}

	delete[]d;
}


void TWfolder::WriteSensitivityBinary(std::ostream &os, char var, char pert) {

	int npert = ZenN(&worhp_o, &pert);
	int nvar = ZenN(&worhp_o, &var);

	std::cout << "D" << var << "/D" << pert << ": " << nvar << " x " << npert << std::endl;
        auto d = new double[nvar]; // Das sollte gross genug sein...

        for (int i = 1; i <= npert; i++) {

		ZenGetD(&worhp_o, &worhp_w, &worhp_p, &worhp_c, &var, &pert, &i, d);

		os.write((char*)d, sizeof(double) * nvar);
	}

	delete[]d;
}


void TWfolder::GetSensitivity(char var, char pert, int i, double *d) {

	ZenGetD(&worhp_o, &worhp_w, &worhp_p, &worhp_c, &var, &pert, &i, d);
}


bool TWfolder::step() {

	bool ret = true;

	for (TransWorhpProblem* ph : phases) {
		// fuer jede Phase einmal step() aufrufen
		if (!ph->step()) {
			ret = false;
		}
	}

	return ret;
}


void TWfolder::PrintTimes(std::ostream &os) {
	//     worhp_timing.PrintTimes(os);
	os << worhp_timing.timing[0].sumUser;
}

}
