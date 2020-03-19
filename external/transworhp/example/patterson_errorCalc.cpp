/*----------------------------------------------------------------
 *
 *  Example: Example 1 (Patterson, Hager, Rao, Benson, Huntington)
 *  unifedFrameworkAAS.pdf
 *
 * y*  = 4./(1+3.*exp(t))
 * u*  = y* /2
 * mu* = exp(2*ln(1+3*exp(t))-t)/(exp(-5)+6+9*exp(5))
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

#include <algorithm>
#include <numeric>

class calcError : public tw::TransWorhpProblem {
public:

	calcError(const tw::TWdimension &TWdata) : TransWorhpProblem(TWdata) {}

	void localinit() override {
		solver->LinearTimeAxis(0,5);
	}

	void selectWindows(tw::Viewer *viewer) override {
		viewer->AddStateView(0,"y");
		viewer->AddControlView(0,"u");
	}

	void x_init(double *x, int i, int dis) override {
		//x[0] = 4./(1+3.*exp(solver->T[i]));
		x[0] = -1.0/(dis-1)*i + 1.0;
	}

	void u_init(double *u, int i, int dis) override {
		//u[0] = (4./(1+3.*exp(solver->T[i])))/2.;
		u[0] = -0.5/(dis-1)*i + 0.5;
	}

	double obj() override {
		return -x(n_dis-1,0);
	}

	bool obj_structure(tw::DiffStructure &s) override {
		s(0,x_index(n_dis-1,0));
		return true;
	}

	bool obj_diff(tw::DiffStructure &s) override {
		s(0,x_index(n_dis-1,0)) = -1.;
		return true;
	}

	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = -x[0]+x[0]*u[0]-u[0]*u[0];
	}

	bool ode_structure(tw::DiffStructure &s) override {
		s(0,x_indexode(0));
		s(0,u_indexode(0));
		return true;
	}

	bool ode_diff(tw::DiffStructure& s, double t, const double* x, const double* u, const double* p) override {
		s(0,x_indexode(0)) = -1. + u[0];
		s(0,u_indexode(0)) = x[0] - 2.*u[0];
		return true;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		x_low[x_index(0,0)] = x_upp[x_index(0,0)] = 1.0;
	}
};


/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);

	tw::TWdimension TWdim;
	TWdim.ID = "Beispiel mit analytischer Lsg";
	TWdim.n_ode = 1;
	TWdim.n_ctrl = 1;

	const std::vector<tw::TransWORHP_type> solverID = {/*tw::TransWORHP_type::fullDiscretization,*/ tw::TransWORHP_type::pseudospectral};
	const std::vector<tw::PMnodeType> nodeType = {tw::PMnodeType::legendre_gauss, tw::PMnodeType::legendre_lobatto /*tw::PMnodeType::chebyshev_maxima*/};

	for (auto solver : solverID) {

	twparameter.solver = solver;

	for (auto node : nodeType) {

		twparameter.pm_nodes_type = node;

		std::ofstream fs("literatur_sol"+std::to_string(static_cast<int>(solver))+"no"+std::to_string(static_cast<int>(node))+"_loeschen.txt");
		fs << "dis xMaxErr xAveErr uMaxErr uAveErr iter" << std::endl;

		for (int dis = 4; dis < 31; dis+=1) {

			TWdim.n_dis = dis;
			TWdim.setMultinodes(2);

			tw::TWfolder folder(&twparameter,0);

			calcError ph(TWdim);
			ph.setSolver(&twparameter);
			folder.Add(&ph);

			std::unique_ptr<tw::Viewer> viewer;
			if (twparameter.PLOT) viewer = std::unique_ptr<tw::Viewer>(new tw::Viewer(&twparameter));

			folder.Init();
			folder.Init(viewer.get());

			// Werte aus paper
			folder.worhp_p.TolOpti = 1e-15;
			folder.worhp_p.TolFeas = 2e-15;

			folder.Loop();

			std::vector<double> errX(TWdim.n_dis-2);
			std::vector<double> errU(TWdim.n_dis-2);

			for (int i = 1; i < TWdim.n_dis-1; i++) {
				//std::cout << ph.solver->T[i] << std::endl;
				errX[i-1] = std::max(std::abs(ph.x(i,0) - 4./(1+3.*exp(ph.solver->T[i]))),1e-16);
				errU[i-1] = std::max(std::abs(ph.u(i,0) - (4./(1+3.*exp(ph.solver->T[i])))/2.),1e-16);
			}
			/*
			std::cout << "errX" << std::endl;
			for (auto ele : errX) {
				std::cout << ele << std::endl;
			}
			std::cout << "errU" << std::endl;
			for (auto ele : errU) {
				std::cout << ele << std::endl;
			}
			*/
			//std::cout << "dis=" << dis << std::endl;
			//std::cout << "max errX: ";

			if (folder.worhp_c.status == OptimalSolution) {
				fs << std::setw(2) << dis;
				std::cout << std::setw(15) << *std::max_element(errX.begin(),errX.end());
				fs << std::setw(15) << *std::max_element(errX.begin(),errX.end());
				std::cout << std::setw(15) << std::accumulate(errX.begin(), errX.end(), 0.0)/errX.size();
				fs << std::setw(15) << std::accumulate(errX.begin(), errX.end(), 0.0)/errX.size();
				std::cout << std::setw(15) << *std::max_element(errU.begin(),errU.end());
				fs << std::setw(15) << *std::max_element(errU.begin(),errU.end());
				std::cout << std::setw(15) << std::accumulate(errU.begin(), errU.end(), 0.0)/errU.size();
				fs << std::setw(15) << std::accumulate(errU.begin(), errU.end(), 0.0)/errU.size();
				std::cout << std::setw(5) << folder.worhp_w.MajorIter;
				fs << std::setw(5) << folder.worhp_w.MajorIter;
				std::cout << std::endl;
				fs << std::endl;
			}
			//std::cout << "\t";
			//std::cout << "max errU: ";
			//std::cout << std::setw(12) << *std::max_element(errU.begin(),errU.end()) << std::endl;
		}

	fs.close();

	if (solver != tw::TransWORHP_type::pseudospectral) break;

	}
	}

	return 0;
}
