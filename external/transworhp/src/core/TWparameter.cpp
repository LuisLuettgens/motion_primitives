#include "TWparameter.h"

#include "TransWORHP.h"
#include "TWconsole.h"
#include "butcher.h"
#include "twstatus.h"

#include "../base/defines.h"

#include "conversion.h"

#include <string>
#include <cstring>
#include <iostream>

#ifdef WIN32
#include <windows.h>
#endif

namespace tw {

TWparameter::TWparameter(std::string filename) :

	countParams(0), setParams(0),

	solver(TransWORHP_type::fullDiscretization),
	butchertableau(0), stepsize(1e-3),
	abserr(1e-6), relerr(1e-6),
	linInter(true), parallel(false),
	hessianvalues(0), hessianstructure(3),
	paramfile_transworhp(std::move(filename)),
	USERDF(-1),USERDG(-1),USERHM(-1),
	showDF(true), showDG(true), showHM(true),
	showGrid(false),
	showOBJstructure(true), showODEstructure(true), showRANDstructure(true), showNEBENstructure(true),
	eps(1e-6),
	NDIS(11), PLOT(1),

	pm_smoothmode(false), pm_displaypoints(64),
	pm_nodes_type(PMnodeType::chebyshev_lobatto),

	meshref_mod(0),
	meshref_err_mod(0),
	meshref_M1(5),
	meshref_R(0),
	meshref_M(-1),
	meshref_maxIter(100),
	meshref_kappa(0.1),
	meshref_tol(1e-8),

	meshref_VER(0),

	meshref_PLOT_SW(0),
	meshref_PLOT_ERR(0),
	meshref_PLOT_MP(0),
	meshref_PLOT_LAMBDA(0),

	xml(ReadParams(paramfile_transworhp))
{
	if (xml) {
		parseXML(xml);
	}
}

TWparameter::~TWparameter()
{

}

std::map<std::string,std::string> TWparameter::Arguments(int argv, char* argc[]) {

	std::map<std::string,std::string> ret;

	for (int i = 1; i < argv; i++) {
		if (strcmp(argc[i], "-n") == 0) {
			i++;
			NDIS = atoi(argc[i]);
		}
		else if (strcmp(argc[i], "+p") == 0) {
			PLOT = 1;
		}

		else if (strcmp(argc[i], "-p") == 0) {
			PLOT = 0;
		}

		else if (strcmp(argc[i], "-m") == 0) {
			i++;
			multinode = ToIntArray(argc[i]);
		}

		else if (strcmp(argc[i], "-fd") == 0) {
			solver = TransWORHP_type::fullDiscretization;
		}
		else if (strcmp(argc[i], "-expl") == 0) {
			solver = TransWORHP_type::multipleShooting;
		}
		else if (strcmp(argc[i], "-pm") == 0) {
			solver = TransWORHP_type::pseudospectral;
		}

		else if (strncmp(argc[i], "--", 2) == 0) {
			std::string key( argc[i] );
			key = key.substr(2);
			i++;
			ret.insert(std::pair<std::string,std::string>(std::move(key),argc[i]));
		}
	}
	return ret;
}


std::unique_ptr<XMLNode> TWparameter::ReadParams(const std::string &filename) {

	XMLParser parser;
	XMLNode *node = parser.Parse(filename);
	std::unique_ptr<XMLNode> x_clone;

	if (node) {
		XMLNode *xmlmain = parser.GetDocument();
		x_clone = std::unique_ptr<XMLNode>(xmlmain->Clone());
	} else {
		DebugStream d(standard_textoutputtype);
		d << beginfile;
		parser.GetError(d);
		d << endfile;
		std::cout << d.GetString();

#ifdef WIN32
		std::string a = "Parameter file \"" + filename + "\" wasn't found!";
		MessageBoxA(0, a.c_str(), "TransWORHP", MB_ICONEXCLAMATION);
#endif


#ifdef TRANSWORHP_GRAPHICS
		exit(1);
#else
		return 0;
#endif

	}

	return x_clone;
}



void TWparameter::parseXML(const std::unique_ptr<XMLNode> &xmlmain) {

	countParams += 4;
	XMLNode *n = xmlmain->GetFirstChild("WORHP");
	if (n) {
		std::string s( n->GetAttribute("param") );
		if (!s.empty()) {
			paramfile_worhp = std::move(s);
			setParams++;
		}

		XMLNode *xml = n->GetFirstChild("USERDF");
		if (xml) {
			USERDF = ToInt(xml->GetText());
			setParams++;
		}

		xml = n->GetFirstChild("USERDG");
		if (xml) {
			USERDG = ToInt(xml->GetText());
			setParams++;
		}

		xml = n->GetFirstChild("USERHM");
		if (xml) {
			USERHM = ToInt(xml->GetText());
			setParams++;
		}
	}

	countParams += 2;
	n = xmlmain->GetFirstChild("SOLVER");
	if (n) {
		std::string s( n->GetText() );

		if (s == "fd") {
			solver = TransWORHP_type::fullDiscretization;
			setParams++;
		}
		else if (s == "expl") {
			solver = TransWORHP_type::multipleShooting;
			setParams++;
		}
		else if (s == "pm") {
			solver = TransWORHP_type::pseudospectral;
			setParams++;
		}

		std::string s2( n->GetAttribute("ndis") );
		if (!s2.empty()) {
			NDIS = ToInt(s2);
			setParams++;
		}
	}

	countParams++;
	n = xmlmain->GetFirstChild("DISCRETIZATION");
	if (n) {
		std::string s( n->GetText() );

		int a = ToInt(s);

		if (s == "Euler" || a == 0) {
			twdiscretization = TWdiscretization(TWdiscretizationType::Euler,1,0);
		}
		else if (s == "Trapez" || a == 1) {
			twdiscretization = TWdiscretization(TWdiscretizationType::Trapez,1,0);
		}
		else if (s == "HermiteSimpson" || a == 2) {
			twdiscretization = TWdiscretization(TWdiscretizationType::HermiteSimpson,2,1);
		}
		else if (s == "Lobatto" || a == 3) {
			twdiscretization = TWdiscretization(TWdiscretizationType::Lobatto,3,2);

		}
		setParams++;
	}

	countParams += 6;
	n = xmlmain->GetFirstChild("BUTCHER");
	if (n) {
		std::string s( n->GetText() );
		butchertableau = ToInt(s);
		setParams++;

		s = n->GetAttribute("stepsize");
		if (!s.empty()) {
			stepsize = ToDouble(s);
			setParams++;
		}

		/* geht nicht, da Butcher::name noch leer
		for (int i = 0; i < Butcher::names.size(); i++) {
			if (s == Butcher::names[i]) {
				butchertableau = i;
				break;
			}
		}
		*/

		s = n->GetAttribute("abserr");
		if (!s.empty()) {
			abserr = ToDouble(s);
			setParams++;
		}
		s = n->GetAttribute("relerr");
		if (!s.empty()) {
			relerr = ToDouble(s);
			setParams++;
		}
		s = n->GetAttribute("linInter");
		if (!s.empty()) {
			linInter = ToBool(s);
			setParams++;
		}
		s = n->GetAttribute("parallel");
		if (!s.empty()) {
			parallel = ToBool(s);
			setParams++;
		}
	}

	countParams += 2;
	n = xmlmain->GetFirstChild("HESSIAN");
	if (n) {

		XMLNode *nn = n->GetFirstChild("STRUCTURE");

		if (nn) {
			std::string s( nn->GetText() );
			hessianstructure = ToInt(s);
			setParams++;

			if (s == "Diagonal") {
				hessianstructure = 0;
			} else if (s == "Full") {
				hessianstructure = 1;
			} else if (s == "Odeblocks2") {
				hessianstructure = 2;
			} else if (s == "Odeblocks") {
				hessianstructure = 3;
			}
		}

		nn = n->GetFirstChild("VALUES");

		if (nn) {
			std::string s( nn->GetText() );
			hessianvalues = ToInt(s);
			setParams++;

			if (s == "DiffDiff") {
				hessianvalues = 0;
			} else if (s == "DiffDG") {
				hessianvalues = 1;
			}
		}
	}

	countParams++;
	n = xmlmain->GetFirstChild("FINITEDIFF");
	if (n) {
		std::string s( n->GetText() );
		if (!s.empty()) {
			eps = ToDouble(s);
			setParams++;
		}
	}

	countParams += 2;
	n = xmlmain->GetFirstChild("GUI");
	if (n) {
		std::string s( n->GetText() );
		if (!s.empty()) {
			PLOT = ToInt(s);
			setParams++;
		}

		s = n->GetAttribute("grid");
		if (!s.empty()) {
			showGrid = ToBool(s);
			setParams++;
		}
	}

	countParams += 3;
	n = xmlmain->GetFirstChild("PLOT");
	if (n) {
		XMLNode *xml = n->GetFirstChild("SPARSITY");
		while (xml) {
			std::string s( xml->GetText() );
			if (s == "DF") {
				showDF = true;
				setParams++;
			} else if (s == "DG") {
				showDG = true;
				setParams++;
			} else if (s == "HM") {
				showHM = true;
				setParams++;
			}
			xml = n->GetNextChild("SPARSITY");
		}
	}

	countParams += 4;
	n = xmlmain->GetFirstChild("SHOWSTRUCTURE");
	if (n) {
		XMLNode *nn = n->GetFirstChild("OBJ");
		if (nn) {
			std::string s( nn->GetText() );
			showOBJstructure = ToBool(s);
			setParams++;
		}
		nn = n->GetFirstChild("ODE");
		if (nn) {
			std::string s( nn->GetText() );
			showODEstructure = ToBool(s);
			setParams++;
		}
		nn = n->GetFirstChild("RAND");
		if (nn) {
			std::string s( nn->GetText() );
			showRANDstructure = ToBool(s);
			setParams++;
		}
		nn = n->GetFirstChild("NEBEN");
		if (nn) {
			std::string s( nn->GetText() );
			showNEBENstructure = ToBool(s);
			setParams++;
		}
	}

#ifdef TRANSWORHP_GRAPHICS
	n = xmlmain->GetFirstChild("WINDOW");
	if (n) {
		twwindow.ParseXML(n, &countParams, &setParams);
	}
#endif

	n = xmlmain->GetFirstChild("CONSOLE");
	if (n) {
		twconsole.ParseXML(n, &countParams, &setParams);
	}

	countParams += 3;
	n = xmlmain->GetFirstChild("PSEUDOSPECTRAL");
	if (n) {
		XMLNode *nn = n->GetFirstChild("NODE_TYPE");
		if (nn) {
			std::string s( nn->GetText() );
			if (!s.empty()) {

				switch (ToInt(s)) {
					case 0:
						pm_nodes_type = PMnodeType::chebyshev_lobatto;
						setParams++;
						break;
					case 1:
						pm_nodes_type = PMnodeType::chebyshev_maxima;
						setParams++;
						break;
					case 2:
						pm_nodes_type = PMnodeType::legendre_lobatto;
						setParams++;
						break;
					case 3:
						pm_nodes_type = PMnodeType::legendre_gauss;
						setParams++;
						break;
				}
			}
		}
		nn = n->GetFirstChild("SMOOTHMODE");
		if (nn) {
			std::string s( nn->GetText() );
			if (!s.empty()) {
				pm_smoothmode = ToBool(s);
				setParams++;
			}
		}
		nn = n->GetFirstChild("DISPLAYPOINTS");
		if (nn) {
			std::string s( nn->GetText() );
			if (!s.empty()) {
				pm_displaypoints = ToInt(s);
				setParams++;
			}
		}
	}

	/* Gitteranpassung */
	countParams += 13;
	n = xmlmain->GetFirstChild("MESHREF");
	if (n) {

		std::string s( n->GetAttribute("modus") );
		if (!s.empty()) {
			meshref_mod = ToInt(s);
			setParams++;
		}

		XMLNode *nn = n->GetFirstChild("FEHLER_MOD");
		if (nn) {
			std::string s( nn->GetText() );
			if (!s.empty()) {
				meshref_err_mod = ToInt(s);
				setParams++;
			}
		}

		nn = n->GetFirstChild("M1");
		if (nn) {
			std::string s( nn->GetText() );
			if (!s.empty()) {
				meshref_M1 = std::abs(ToInt(s));
				setParams++;
			} else { // Standard: wenn dieser Parameter leer ist, wird -1 gesetzt
				meshref_M1 = -1;
				setParams++; //deshalb auch erhoehen
			}
		}

		nn = n->GetFirstChild("R");
		if (nn) {
			std::string s( nn->GetText() );
			if (!s.empty()) {
				meshref_R = std::abs(ToInt(s));
				setParams++;
			}
		}

		nn = n->GetFirstChild("KAPPA");
		if (nn) {
			std::string s( nn->GetText() );
			if (!s.empty()) {
				meshref_kappa = std::abs(ToDouble(s));
				setParams++;
			}
		}

		nn = n->GetFirstChild("M");
		if (nn) {
			std::string s( nn->GetText() );
			if (!s.empty()) {
				meshref_M = ToInt(s);
				setParams++;
			} else {
				// M leer => Default VALUE
				meshref_M = -1;
				setParams++;
			}
		}

		nn = n->GetFirstChild("maxIter");
		if (nn) {
			std::string s( nn->GetText() );
			if (!s.empty()) {
				meshref_maxIter = ToInt(s);
				setParams++;
			}
		}

		nn = n->GetFirstChild("VERLAUF");
		if (nn) {
			std::string s( nn->GetText() );
			if (!s.empty()) {
				meshref_VER = ToInt(s);
				setParams++;
			}
		}

		nn = n->GetFirstChild("TOL");
		if (nn) {
			std::string s( nn->GetText() );
			if (!s.empty()) {
				meshref_tol = ToDouble(s);
				setParams++;
			}
		}

		nn = n->GetFirstChild("PLOTS");
		if (nn) {
			XMLNode *nnn = nn->GetFirstChild("SCHRITTWEITE");
			if (nnn) {
				std::string s( nnn->GetText() );
				if (!s.empty()) {
					meshref_PLOT_SW = ToInt(s);
					setParams++;
				}
			}
			nnn = nn->GetFirstChild("FEHLER");
			if (nnn) {
				std::string s( nnn->GetText() );
				if (!s.empty()) {
					meshref_PLOT_ERR = ToInt(s);
					setParams++;
				}
			}
			nnn = nn->GetFirstChild("GITTERPUNKTE");
			if (nnn) {
				std::string s( nnn->GetText() );
				if (!s.empty()) {
					meshref_PLOT_MP = ToInt(s);
					setParams++;
				}
			}
			nnn = nn->GetFirstChild("LAMBDA");
			if (nnn) {
				std::string s( nnn->GetText() );
				if (!s.empty()) {
					meshref_PLOT_LAMBDA = ToInt(s);
					setParams++;
				}
			}
		}
	}

	MyStatus("TransWORHP", "Used TransWORHP parameter file: " + paramfile_transworhp , Status::NORMAL);
	MyStatus("TransWORHP", "Read " + std::to_string(setParams) + " of " +
	                           std::to_string(countParams) + " parameters",
	         (countParams == setParams) ? Status::NORMAL : Status::WARN);
}

void TWparameter::printParams() const {
	MyStatus("TransWORHP", "Print parameters:", Status::NORMAL);
	MyStatus("TransWORHP", "solver: " + std::to_string(static_cast<int>(solver)), Status::NORMAL);
	MyStatus("TransWORHP", "twdiscretization: " + std::to_string(static_cast<int>(twdiscretization.type)), Status::NORMAL);
	MyStatus("TransWORHP", "butchertableau: " + std::to_string(butchertableau), Status::NORMAL);
	MyStatus("TransWORHP", "stepsize: " + std::to_string(stepsize), Status::NORMAL);
	MyStatus("TransWORHP", "abserr: " + std::to_string(abserr), Status::NORMAL);
	MyStatus("TransWORHP", "relerr: " + std::to_string(relerr), Status::NORMAL);
	MyStatus("TransWORHP", "linInter: " + std::to_string(linInter), Status::NORMAL);
	MyStatus("TransWORHP", "parallel: " + std::to_string(parallel), Status::NORMAL);
	MyStatus("TransWORHP", "hessianvalues: " + std::to_string(hessianvalues), Status::NORMAL);
	MyStatus("TransWORHP", "hessianstructure: " + std::to_string(hessianstructure), Status::NORMAL);
	MyStatus("TransWORHP", "paramfile_worhp: " + paramfile_worhp, Status::NORMAL);
	MyStatus("TransWORHP", "paramfile_worhp: " + paramfile_transworhp, Status::NORMAL);
	MyStatus("TransWORHP", "USERDF: " + std::to_string(USERDF), Status::NORMAL);
	MyStatus("TransWORHP", "USERDG: " + std::to_string(USERDG), Status::NORMAL);
	MyStatus("TransWORHP", "USERHM: " + std::to_string(USERHM), Status::NORMAL);
	MyStatus("TransWORHP", "showDF: " + std::to_string(showDF), Status::NORMAL);
	MyStatus("TransWORHP", "showDG: " + std::to_string(showDG), Status::NORMAL);
	MyStatus("TransWORHP", "showHM: " + std::to_string(showHM), Status::NORMAL);
	MyStatus("TransWORHP", "showGrid: " + std::to_string(showGrid), Status::NORMAL);
	MyStatus("TransWORHP", "eps: " + std::to_string(eps), Status::NORMAL);
	MyStatus("TransWORHP", "NDIS: " + std::to_string(NDIS), Status::NORMAL);
	MyStatus("TransWORHP", "PLOT: " + std::to_string(PLOT), Status::NORMAL);

	MyStatus("TransWORHP", "Pseudospectral:", Status::NORMAL);
	MyStatus("TransWORHP", "pm_smoothmode: " + std::to_string(pm_smoothmode), Status::NORMAL);
	MyStatus("TransWORHP", "pm_displaypoints: " + std::to_string(pm_displaypoints), Status::NORMAL);
	MyStatus("TransWORHP", "pm_nodes_type: " + std::to_string(static_cast<int>(pm_nodes_type)), Status::NORMAL);

	MyStatus("TransWORHP", "Mesh Refinement:", Status::NORMAL);
	MyStatus("TransWORHP", "meshref_mod: " + std::to_string(meshref_mod), Status::NORMAL);
	MyStatus("TransWORHP", "meshref_err_mod: " + std::to_string(meshref_err_mod), Status::NORMAL);
	MyStatus("TransWORHP", "meshref_M1: " + std::to_string(meshref_M1), Status::NORMAL);
	MyStatus("TransWORHP", "meshref_R: " + std::to_string(meshref_R), Status::NORMAL);
	MyStatus("TransWORHP", "meshref_M: " + std::to_string(meshref_M), Status::NORMAL);
	MyStatus("TransWORHP", "meshref_maxIter: " + std::to_string(meshref_maxIter), Status::NORMAL);
	MyStatus("TransWORHP", "meshref_kappa: " + std::to_string(meshref_kappa), Status::NORMAL);
	MyStatus("TransWORHP", "meshref_tol: " + std::to_string(meshref_tol), Status::NORMAL);
	MyStatus("TransWORHP", "meshref_VER: " + std::to_string(meshref_VER), Status::NORMAL);
	MyStatus("TransWORHP", "meshref_PLOT_SW: " + std::to_string(meshref_PLOT_SW), Status::NORMAL);
	MyStatus("TransWORHP", "meshref_PLOT_ERR: " + std::to_string(meshref_PLOT_ERR), Status::NORMAL);
	MyStatus("TransWORHP", "meshref_PLOT_MP: " + std::to_string(meshref_PLOT_MP), Status::NORMAL);
	MyStatus("TransWORHP", "meshref_PLOT_LAMBDA: " + std::to_string(meshref_PLOT_LAMBDA), Status::NORMAL);
}

}
