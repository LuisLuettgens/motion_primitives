/*----------------------------------------------------------------
 *
 *  Example: Splineproblem mit Fehlerbestimmung
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

#include <algorithm>
#include <numeric>
#include <chrono>

class SplinePhase : public tw::TransWorhpProblem {
public:

	SplinePhase(const tw::TWdimension &TWdim) : TransWorhpProblem(TWdim) {}

	void selectWindows(tw::Viewer *viewer) override {

		viewer->AddStateView(0,"x_1");
		viewer->AddStateView(1,"x_2");
		viewer->AddStateView(2,"x_3");
		
		viewer->AddControlView(0,"u");
	}
	
	double obj() override {
		return .5 * x(n_dis-1,2);
	}

	bool obj_structure(tw::DiffStructure &s) override {
		s(0, x_index(n_dis-1,2));
		return true;
	}

	bool obj_diff(tw::DiffStructure &s) override {
		s(0,x_index(n_dis-1,2)) = .5;
		return true;
	}

	void x_init(double *x, int i, int dis) {
		//x[0] = -solver->T[i]*solver->T[i]+solver->T[i];
		//x[1] = -2*solver->T[i]+1;
	}

	void u_init(double *u, int i, int dis) {
		u[0] = -2.0;
	}

	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		dx[0] = x[1];
		dx[1] = u[0];
		dx[2] = u[0]*u[0];
	}

	bool ode_structure(tw::DiffStructure &s) override {
		s(0, x_indexode(1)); // dx[0] / dx[1]
		s(1, u_indexode(0)); // dx[1] / du[0]
		s(2, u_indexode(0)); // dx[2] / du[0]
		return true;
	}

	bool ode_diff(tw::DiffStructure &s, double t, const double *x, const double *u, const double *p) override {
		s(0, x_indexode(1))= 1.0;
		s(1, u_indexode(0))= 1.0;
		s(2, u_indexode(0))= 2*u[0];
		return true;
	}

	void var_boundary(double *x_low, double *x_upp) override {
		x_low[x_index( 0,0 ) ] = x_upp[x_index( 0,0 ) ] = 0.0;
		x_low[x_index( 0,1 ) ] = x_upp[x_index( 0,1 ) ] = 1.0;
		x_low[x_index( 0,2 ) ] = x_upp[x_index( 0,2 ) ] = 0.0;

		x_low[x_index( n_dis-1,0 ) ] = x_upp[x_index( n_dis-1,0 ) ] = 0.0;
		x_low[x_index( n_dis-1,1 ) ] = x_upp[x_index( n_dis-1,1 ) ] = -1.0;
	}
};

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	std::map<std::string,std::string> args = twparameter.Arguments(argv,argc);

	tw::TWdimension TWdim;
	TWdim.ID = "Spline";
	TWdim.n_ode = 3;
	TWdim.n_ctrl = 1;

	const std::vector<tw::TransWORHP_type> solverID = {tw::TransWORHP_type::fullDiscretization, tw::TransWORHP_type::multipleShooting, tw::TransWORHP_type::pseudospectral};
	const std::vector<tw::PMnodeType> nodeType = {tw::PMnodeType::legendre_gauss, tw::PMnodeType::legendre_lobatto/*, tw::PMnodeType::chebyshev_maxima*/};

	for (auto solver : solverID) {

	twparameter.solver = solver;
	std::cout << "solver=" << static_cast<int>(solver) << std::endl;

	for (auto node : nodeType) {

	twparameter.pm_nodes_type = node;
	if (solver == tw::TransWORHP_type::pseudospectral) std::cout << "nodeType=" << static_cast<int>(node) << std::endl;

	std::ofstream fs("spline_sol"+std::to_string(static_cast<int>(solver))+"no"+std::to_string(static_cast<int>(node))+"_fidiff.txt");
	fs << "dis x0MaxErr x0AveErr x1MaxErr x1AveErr uMaxErr uAveErr obj time" << std::endl;

	for (int dis = 2; dis < 41; dis++) {

		TWdim.n_dis = dis;
		//TWdim.setMultinodes(dis/2);
		TWdim.setMultinodes(2);

		tw::TWfolder folder(&twparameter,0);

		std::unique_ptr<SplinePhase> ph = std::unique_ptr<SplinePhase>(new SplinePhase(TWdim));
		ph->setSolver(&twparameter);
		folder.Add(ph.get());

		std::unique_ptr<tw::Viewer> viewer;
		if (twparameter.PLOT) viewer = std::unique_ptr<tw::Viewer>(new tw::Viewer(&twparameter));

		folder.Init();
		folder.Init(viewer.get());

		//WORHP Parameter anpassen
		folder.worhp_p.TolOpti = 1e-11;
		folder.worhp_p.TolFeas = 1e-11;

		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
		folder.Loop();
		std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();

		std::vector<double> errX0(TWdim.n_dis-2);
		std::vector<double> errX1(TWdim.n_dis-2);
		std::vector<double> errU(TWdim.n_dis-2);

		for (int i = 1; i < TWdim.n_dis-1; i++) {
			errX0[i-1] = std::abs(ph->x(i,0) - (-(ph->solver->T[i])*(ph->solver->T[i])+ph->solver->T[i]));
			errX1[i-1] = std::abs(ph->x(i,1) - (-2.*(ph->solver->T[i])+1.0));
			errU[i-1] = std::abs(ph->u(i,0) - (-2.));
		}
		//std::cerr << std::endl;

		if (folder.worhp_c.status == OptimalSolution) {
		if (!errX0.empty()) {
			std::cout << std::setw(2) << dis << " : ";
			fs << std::setw(2) << dis;
			std::cout << std::setw(15) << *std::max_element(errX0.begin(),errX0.end());
			fs << std::setw(15) << *std::max_element(errX0.begin(),errX0.end());
			std::cout << std::setw(15) << std::accumulate(errX0.begin(), errX0.end(), 0.0)/errX0.size();
			fs << std::setw(15) << std::accumulate(errX0.begin(), errX0.end(), 0.0)/errX0.size();
			std::cout << std::setw(15) << *std::max_element(errX1.begin(),errX1.end());
			fs << std::setw(15) << *std::max_element(errX1.begin(),errX1.end());
			std::cout << std::setw(15) << std::accumulate(errX1.begin(), errX1.end(), 0.0)/errX1.size();
			fs << std::setw(15) << std::accumulate(errX1.begin(), errX1.end(), 0.0)/errX1.size();
			std::cout << std::setw(15) << *std::max_element(errU.begin(),errU.end());
			fs << std::setw(15) << *std::max_element(errU.begin(),errU.end());
			std::cout << std::setw(15) << std::accumulate(errU.begin(), errU.end(), 0.0)/errU.size();
			fs << std::setw(15) << std::accumulate(errU.begin(), errU.end(), 0.0)/errU.size();
			std::cout << std::setw(15) << folder.worhp_o.F;
			fs << std::setw(15) << folder.worhp_o.F;
			std::cout << std::setw(15) << duration;
			fs << std::setw(15) << duration;
			std::cout << std::endl;
			fs << std::endl;
		}
		}
	}

	fs.close();

	if (solver != tw::TransWORHP_type::pseudospectral) break;
	}
	}

	return 0;
}
