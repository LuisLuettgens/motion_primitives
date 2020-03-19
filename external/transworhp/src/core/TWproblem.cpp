#include "TWproblem.h"

#include "TransWORHP.h"
#include "FullDisTransWORHP.h"
#include "ExplTransWORHP.h"
#include "LobattoPmTransWORHP.h"
#include "GaussPmTransWORHP.h"
#include "TWparameter.h"
#include "twstatus.h"

#ifdef TRANSWORHP_GRAPHICS
#include "Viewer.h"
#endif

#include <utility>

namespace tw {

TWdimension::TWdimension()
	: n_dis(2),
	n_ode(0), n_ctrl(0), n_param(0),
	n_rand(0), n_neben(0),
	n_integral(0), n_zen(0) {

}

TWdimension::TWdimension(std::string name, int dis, int ode, int ctrl, int param, int rand, int neben, int integral, int zen)
	: ID(std::move(name)),
	n_dis(dis),
	n_ode(ode), n_ctrl(ctrl), n_param(param),
	n_rand(rand), n_neben(neben),
	n_integral(integral), n_zen(zen) {

	multinode = {0, n_dis-1};
}

/** setzt Mehrzielknoten
  * @param multi Anzahl der Mehrzielknoten
  */
void TWdimension::setMultinodes(int multi) {

	// schlechte Eingaben abfangen
	if (n_dis < 2) {
		MyStatus("TransWORHP", "setMultinodes(..): set n_dis first!", Status::ERR);
		exit(1);
	}
	if (multi > n_dis) {
		MyStatus("TransWORHP", "too many multinodes! set to n_dis=" + std::to_string(n_dis),
		         Status::ERR);
		multi = n_dis;
	} else if (multi < 2) {
		MyStatus("TransWORHP", "set at least 2 multinodes!", Status::ERR);
		multi = 2;
	}

	// loeschen, falls vorher welche drin waren
	multinode.clear();

	multinode.reserve(multi);
	for (int i = 0; i < multi; i++) {
		const int m = (i*(n_dis-1))/(multi-1);
		multinode.push_back(m);
		//std::cout << m << " " ;
	}
}

////////////////////////////////////////////////////////////////////////////

TransWorhpProblem::TransWorhpProblem(const TWdimension &TWdata)
      : solver(nullptr),
	id(TWdata.ID), unique_id(0),
	multinodes(TWdata.multinode),
	box_neben(TWdata.BOXneben),
	n_dis(TWdata.n_dis),
	n_ode(TWdata.n_ode), n_ctrl(TWdata.n_ctrl), n_param(TWdata.n_param),
	n_rand(TWdata.n_rand), n_neben(TWdata.n_neben), n_integral(TWdata.n_integral), n_zen(TWdata.n_zen),
	freetime(false),
	t0(0.0),
	tf(1.0)
{
	if (multinodes.empty()) {
		multinodes = {0, n_dis-1};
	}
}


void TransWorhpProblem::setSolver(TWparameter *twparam) {

	switch (twparam->solver) {
		case  TransWORHP_type::fullDiscretization:
			solver = std::unique_ptr<FullDisTransWorhp>(
				new FullDisTransWorhp(this,twparam));
			break;
		case TransWORHP_type::multipleShooting :
			twparam->twdiscretization = TWdiscretization(TWdiscretizationType::MultipleShooting,0,0);
			twparam->USERHM = -1;
			solver = std::unique_ptr<ExplTransWorhp>(
				new ExplTransWorhp(this,twparam,multinodes,box_neben));
			break;
		case TransWORHP_type::pseudospectral : {
			if (twparam->pm_nodes_type != PMnodeType::legendre_gauss) {
				twparam->twdiscretization = TWdiscretization(TWdiscretizationType::Pseudospectral,0,0);
				twparam->USERHM = -1;
				solver = std::unique_ptr<LobattoPmTransWorhp>(
					new LobattoPmTransWorhp(this,twparam));
			} else {
				twparam->twdiscretization = TWdiscretization(TWdiscretizationType::Pseudospectral,0,0);
				twparam->USERHM = -1;
				solver = std::unique_ptr<GaussPmTransWorhp>(
					new GaussPmTransWorhp(this,twparam));
			}
			break;
			}
		default:
			exit(1);
	}
}

#ifdef TRANSWORHP_GRAPHICS

void TransWorhpProblem::OpenWindows(Viewer*) {}

std::string TransWorhpProblem::GetXTitle(int) {return std::string();}

std::string TransWorhpProblem::GetUTitle(int) {return std::string();}

void TransWorhpProblem::selectWindows(Viewer* viewer) {
	if (viewer) viewer->selectWindows();
}

#endif

void TransWorhpProblem::localinit() {}

bool TransWorhpProblem::obj_structure(DiffStructure&) {return false;}
bool TransWorhpProblem::obj_diff(DiffStructure&) {return false;}
void TransWorhpProblem::integral(double*, double, const double*, const double*, const double*) {}
bool TransWorhpProblem::integral_structure(DiffStructure&) {return false;}
bool TransWorhpProblem::integral_diff(DiffStructure&, double, const double*, const double*, const double*) {return false;}

bool TransWorhpProblem::ode_structure(DiffStructure&) {return false;}
bool TransWorhpProblem::ode_diff(DiffStructure&, double, const double*, const double*, const double*) {return false;}
bool TransWorhpProblem::ode_diff_p(DiffStructure&, double, const double*, const double*, const double*, int) {return false;}

void TransWorhpProblem::x_boundary(double*, double*) {}
void TransWorhpProblem::u_boundary(double*, double*) {}
void TransWorhpProblem::p_boundary(double*, double*) {}
void TransWorhpProblem::var_boundary(double*, double*) {}

void TransWorhpProblem::rand(double*) {}
void TransWorhpProblem::rand_boundary(double*, double*) {}
bool TransWorhpProblem::rand_structure(DiffStructure&) {return false;}
bool TransWorhpProblem::rand_diff(DiffStructure&) {return false;}

void TransWorhpProblem::neben(double*, double, const double*, const double*, const double*) {}
void TransWorhpProblem::neben_boundary(double*, double*) {}
bool TransWorhpProblem::neben_structure(DiffStructure&) {return false;}
bool TransWorhpProblem::neben_diff(DiffStructure&, double, const double*, const double*, const double*) {return false;}
bool TransWorhpProblem::neben_diff_p(DiffStructure&, double, const double*, const double*, const double*, int) {return false;}

void TransWorhpProblem::init() {}

void TransWorhpProblem::x_init(double*, int, int) {}
void TransWorhpProblem::u_init(double*, int, int) {}
void TransWorhpProblem::p_init(double*) {}
void TransWorhpProblem::zen_init(double*) {}

void TransWorhpProblem::terminate() {}
bool TransWorhpProblem::step() {return true;}


double TransWorhpProblem::x(int dis, int ode) const {
	return solver->x(dis,ode);
}
double TransWorhpProblem::u(int dis, int ctrl) const {
	return solver->u(dis,ctrl);
}
double TransWorhpProblem::p(int param) const {
	return solver->p(param);
}
double TransWorhpProblem::x__(int dis, int ode) const {
	return solver->x__(dis,ode);
}
double TransWorhpProblem::u__(int dis, int ctrl) const {
	return solver->u__(dis,ctrl);
}
int TransWorhpProblem::x_index(int dis, int ode) const {
	return solver->x_index(dis,ode);
}
int TransWorhpProblem::x_index__(int dis, int ode) const {
	return solver->x_index__(dis,ode);
}
int TransWorhpProblem::u_index(int dis, int ctrl) const {
	return solver->u_index(dis,ctrl);
}
int TransWorhpProblem::u_index__(int dis, int ctrl) const {
	return solver->u_index__(dis,ctrl);
}
int TransWorhpProblem::p_index(int param) const {
	return solver->p_index(param);
}
int TransWorhpProblem::x_indexode(int ode) const {
	return solver->x_indexode(ode);
}
int TransWorhpProblem::u_indexode(int ctrl) const {
	return solver->u_indexode(ctrl);
}
int TransWorhpProblem::p_indexode(int param) const {
	return solver->p_indexode(param);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ScaleTransWorhpProblem::ScaleTransWorhpProblem(const TWdimension &TWdata) : TransWorhpProblem(TWdata) {

}

void ScaleTransWorhpProblem::ode(double *dx, double t, const double *x, const double *u, const double *p) {

	double utmp[100];
	for (int i=0;i<n_ctrl;i++) {
		utmp[i] = u[i];
	}
	u_scale(utmp);

	SC_ode(dx,t,x,utmp,p);
}
bool ScaleTransWorhpProblem::ode_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) {
	double utmp[100];
	for (int i=0;i<n_ctrl;i++) {
		utmp[i] = u[i];
	}
	u_scale(utmp);

	return SC_ode_diff(s,t,x,utmp,p);
}
bool ScaleTransWorhpProblem::ode_diff_p(DiffStructure &s, double t, const double *x, const double *u, const double *p, int index) {

	double utmp[100];
	for (int i=0;i<n_ctrl;i++) {
		utmp[i] = u[i];
	}
	u_scale(utmp);

	return SC_ode_diff_p(s,t,x,utmp,p,index);
}
void ScaleTransWorhpProblem::u_boundary(double* u_low, double* u_upp) {

	u_scale(u_low);
	u_scale(u_upp);

	SC_u_boundary(u_low, u_upp);

	u_unscale(u_low);
	u_unscale(u_upp);
}

void ScaleTransWorhpProblem::u_init(double *u, int i, int dis) {

	u_scale(u);

	SC_u_init(u, i, dis);

	u_unscale(u);

}
void ScaleTransWorhpProblem::x_boundary(double* x_low, double* x_upp) {

	x_scale(x_low);
	x_scale(x_upp);

	SC_x_boundary(x_low, x_upp);

	x_unscale(x_low);
	x_unscale(x_upp);
}

void ScaleTransWorhpProblem::x_init(double *x, int i, int dis) {

	x_scale(x);

	SC_x_init(x, i, dis);

	x_unscale(x);
}

bool ScaleTransWorhpProblem::SC_ode_diff(DiffStructure&, double, const double*, const double*, const double*) {return false;}
bool ScaleTransWorhpProblem::SC_ode_diff_p(DiffStructure&, double, const double*, const double*, const double*, int) {return false;}

void ScaleTransWorhpProblem::SC_x_boundary(double*, double*) {}
void ScaleTransWorhpProblem::SC_u_boundary(double*, double*) {}

void ScaleTransWorhpProblem::SC_x_init(double*, int, int) {}
void ScaleTransWorhpProblem::SC_u_init(double*, int, int) {}


void ScaleTransWorhpProblem::u_scale(double*) {}
void ScaleTransWorhpProblem::u_unscale(double*) {}
void ScaleTransWorhpProblem::x_scale(double*) {}
void ScaleTransWorhpProblem::x_unscale(double*) {}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

EmptyTransWorhpProblem::EmptyTransWorhpProblem(const TWdimension &TWdata) : TransWorhpProblem(TWdata) {}

double EmptyTransWorhpProblem::obj() {
	return 0.0;
}

void EmptyTransWorhpProblem::ode(double*, double, const double*, const double*, const double*) {

}

}
