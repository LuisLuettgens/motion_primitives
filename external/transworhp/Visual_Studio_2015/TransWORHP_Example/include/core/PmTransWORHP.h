#pragma once

#include "TransWORHP.h"

namespace tw {

enum class PMnodeType {chebyshev_lobatto, chebyshev_maxima, legendre_lobatto, legendre_gauss};

class DllExport PmTransWorhp : public TransWorhp {

public:
	
	PmTransWorhp(TransWorhpProblem *ph, TWparameter *twparam);
	virtual ~PmTransWorhp();

	/** Varianten fuer Zugriff */
	double x__(int dis, int ode) const override;
	double u__(int dis, int ctrl) const override;

	/** Indexbestimmung */
	int x_index__(int dis, int ode) const override;
	int u_index__(int dis, int ode) const override;
	
	int x_indexode(int ode) const override;
	int u_indexode(int ctrl) const override;
	int p_indexode(int param) const override;
	
	void GetState(double *x, double t) override;
	void GetControl(double *uu, double t) override;
	
#ifdef TRANSWORHP_GRAPHICS
	void setTemptimeForViewer(std::vector<double> &tmptime) override;
#endif

	void GetBoundaryIndices(std::vector<int> &indices, int d) override;
	
	void LinearTimeAxis(double start, double end) override;

	void ToMATLAB(const std::string& filename) const override;

protected:

	/* L_i(t) */
	virtual double lagrangePoly(int i, double t) const = 0;
	/* d/dt_k L_i(t_k) */
	virtual double lagrangePolyDer(int i, int k) const = 0;
	
	virtual double interpolate(int index, double time) const = 0;
	
	/** time t \in [t0,tf] -> [-1,1] */
	double timeToOneOne(double time) const;
	/** time t \in [-1,1] -> [t0,tf] */
	double timeToTzeroTfinal(double time) const;
	
	/** Pseudospectral-Differantial-Matrix */
	std::vector<std::vector<double>> PDM_D;
	/** weights for barycentric form (Lagrange Poly) */
	std::vector<double> weights;
	
	/** time \in [-1,1] -> Legendre or Chebychev */
	std::vector<double> time;
	
	/** smooth time vector */
	std::vector<double> displayTime;
	/** X_0, U_0, X_1, U_1, ... */
	std::vector<double> displayXU;
	
	/** smoothMode for smooth plot */
	bool smoothMode;
	/** #points to plot */
	int displayPoints;
};

}
