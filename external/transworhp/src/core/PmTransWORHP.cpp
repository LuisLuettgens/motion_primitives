#include "PmTransWORHP.h"

#include "TWproblem.h"
#include "TWfolder.h"
#include "butcher.h"

#ifdef TRANSWORHP_GRAPHICS
#include "Viewer.h"
#endif

#include <algorithm>

namespace tw {

PmTransWorhp::PmTransWorhp(TransWorhpProblem* ph, TWparameter *twparam)
	: TransWorhp(ph,twparam),
	smoothMode(twparameter->pm_smoothmode), displayPoints(twparameter->pm_displaypoints)
{
	
}


PmTransWorhp::~PmTransWorhp() {

}


void PmTransWorhp::LinearTimeAxis(double start, double end) {

	phase->t0 = start;
	phase->tf = end;
	
	for (int i = 0; i < n_dis; i++) {
		T[i] = timeToTzeroTfinal(time[i]);
	}

	for (int i = 0; i < displayPoints; i++) {
		displayTime[i] = i/(displayPoints-1.)*(phase->tf-phase->t0)+phase->t0;
	}
}


double PmTransWorhp::timeToOneOne(double time_) const {
	return 2.0*time_/(phase->tf-phase->t0) - (phase->tf+phase->t0)/(phase->tf-phase->t0);
}


double PmTransWorhp::timeToTzeroTfinal(double time_) const {
	return time_*(phase->tf-phase->t0)/2.0+(phase->tf+phase->t0)/2.0;
}


double PmTransWorhp::x__(int dis, int ode) const {
	return x(dis,ode);
}


double PmTransWorhp::u__(int dis, int ctrl) const {
	return u(dis,ctrl);
}


int PmTransWorhp::x_index__(int dis, int ode) const {
	return x_index(dis,ode);
}


int PmTransWorhp::u_index__(int dis, int ode) const {
	return u_index(dis,ode);
}


int PmTransWorhp::x_indexode(int ode) const {
	return ode;
}


int PmTransWorhp::u_indexode(int ctrl) const {
	return n_ode + ctrl;
}


int PmTransWorhp::p_indexode(int param) const {
	return n_ode + n_ctrl + param;
}


void PmTransWorhp::GetState(double *xx, double t) {

	//wirklich virtual?

	int i = 0;
	for (; i < n_dis; i++) {
		if (t < T[i]) break;
	}

	for (int j = 0; j < n_ode; j++) {
		if (i > 0 && i < n_dis) {
			xx[j] = x(i-1,j) +  (x(i,j) - x(i-1,j)) / (T[i]-T[i-1])*(t-T[i-1]);
		}
		else if (i >= n_dis) {
			xx[j] = x(n_dis-1,j);
		}
	}
}


void PmTransWorhp::GetControl(double */*uu*/, double /*t*/) {
	std::cout << "GetControl" << std::endl;
	exit(1);
}


#ifdef TRANSWORHP_GRAPHICS
void PmTransWorhp::setTemptimeForViewer(std::vector<double> &tmptime) {

	if (!smoothMode) {
		tmptime.reserve(n_dis);
		for (int i = 0; i < n_dis; i++) {
			// Trafo [-1,1]->[t_0,t_f] fuer Viewer
			tmptime.push_back(timeToTzeroTfinal(time[i]));
		}
	} else {
		tmptime = displayTime;
	}
}
#endif


void PmTransWorhp::GetBoundaryIndices(std::vector<int> &indices, int d) {

	if ( X_low[x_index(0,d)] == X_upp[x_index(0,d)] ) {
		indices.push_back(0);
	}
	if ( X_low[x_index(n_dis-1,d)] == X_upp[x_index(n_dis-1,d)] ) {
		if (smoothMode) {
			indices.push_back(displayPoints-1);
		} else {
			indices.push_back(n_dis-1);
		}
	}
}

void PmTransWorhp::ToMATLAB(const std::string& filename) const {

	std::ofstream of(filename);

	of.setf(std::ios::scientific);
	of << std::setprecision(9);

	of << "% TransWORHP-Result" << std::endl;

	of << "% ";
	for (int i = 0; i < n_param; i++) {
		of << std::setw(20) << p(i);
	}
	of << std::endl;

	if (smoothMode) {
		for (int i = 0; i < displayPoints; i++) {

			if (phase->freetime) {
				of << std::setw(20) << (displayTime[i] * p(0));
			} else {
				of << std::setw(20) << displayTime[i];
			}

			for (int j = 0; j < n_ode; j++) {
				of << std::setw(20) << displayXU[i*(n_ode+n_ctrl) + j];
			}
			for (int j = 0; j < n_ctrl; j++) {
				of << std::setw(20) << displayXU[i*(n_ode+n_ctrl) + n_ode + j];
			}
			of << std::endl;
		}
	} else {
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
}

}
