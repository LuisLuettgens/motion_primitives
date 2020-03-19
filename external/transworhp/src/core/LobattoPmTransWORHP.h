#pragma once

#include "PmTransWORHP.h"

namespace tw {

class DllExport LobattoPmTransWorhp : public PmTransWorhp {

public:

	LobattoPmTransWorhp(TransWorhpProblem* ph, TWparameter *twparam);
	virtual ~LobattoPmTransWorhp();

	/**
	* gibt an, welche Beschraenkung an der Stelle i ist (ODE,RAND,NB,etc).
	* @param Kompomente des G Vektor
	* @return Typ, zB ODE 0, RAND 1, etc
	*/
	std::string type_G(int i) const override;

	/** Zugriffsfunktionen */
	double x(int dis, int ode) const override;
	double u(int dis, int ctrl) const override;
	double p(int param) const override;

	/** Indexbestimmung */
	int x_index(int dis, int ode) const override;
	int u_index(int dis, int ctrl) const override;
	int p_index(int param) const override;

#ifdef TRANSWORHP_GRAPHICS
	std::string setViewerTitle() override;
	void updateViewer(DataStorage* ds, std::vector<double> &tmptime) override;
#endif

	double Objective(double ScaleObj) override;
	void Lagrange() override;

	/** Speicher und Struktur */
	void GetSpace(int delta1, int delta2) override;
	/** Verbindung zu WORHP */
	void Connect(const OptVar &o, const Params &p) override;
	/** Groessen der Opt-Vektoren */
	void init0() override;

	/** Beschraenkungsschranken */
	void Boundary() override;
	/** Beschraenkungen mit anderem G Vektor */
	void Constraints2(double *GG, int DGflag=0) override;

	int DF_structure(WorhpMatrix *DF=nullptr, int offset=0) override;
	void DF_calculate(WorhpMatrix &DF, double ScaleObj) override;

	int DG_structure(const TWfolder *f, WorhpMatrix *DG=nullptr, int offset=0) override;
	void DG_calculate(TWfolder *f, WorhpMatrix &DG) override;

	int HM_structure_ohne_Diag(int, WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix *HM=nullptr, int offset=0) override;
	void HM_calculate1(WorhpMatrix &DF, WorhpMatrix &DG, WorhpMatrix &HM, double ScaleObj, double *Mu) override;

	int Integrate(int btableau) override;
private:
	int integrate(int index);
	int integrateRKF(int index, int &startflag);

	double& setx(int dis, int ode) override;

	void DG_diff_ode(double t, int k, int active_index);
	void DG_diff_ode_p(double t, int l, int dis, int active_index);

	void DG_diff_rand();

	void DG_diff_neben(double t, int dis);
	void DG_diff_neben_p(double t, int l, int dis);

	/** zeros Chebychev polynomial */
	void chebyshevTime();
	/** zeros Chebychev polynomial + {-1,+1} */
	void chebyshevLobattoTime();
	/** extrema Chebychev polynomial */
	void chebyshevExtremaTime();

	/** zeros derivative Legendre polynomial n-1 
	 * rootsOf( d/dt P_{N-1}(t) ) + {-1,+1}
	 */
	void legendreLobattoTime();

	/** zeros Legendre polynomial + {-1,+1} 
	 * rootsOf( P_{N-2}(t) ) + {-1,+1}
	 */
	void legendreLobattoTime2();

	/** fill Pseudospectral Matrix D with the (direct) derivatives of the Lagrange Poly */
	void fillLobattoPseudospectralMatrixDirect();
	/** fill Pseudospectral Matrix D using barycentric form of Lagrange Poly */
	void fillLobattoPseudospectralMatrixBarycentric();

	/* L_i(t) */
	double lagrangePoly(int i, double t) const override;
	/* d/dt_k L_i(t_k) */
	double lagrangePolyDer(int i, int k) const override;

	double interpolate(int index, double time) const override;
	
	void fromMATLAB_impl(std::ifstream &stream) override;

	std::vector<double> lobatto_weights;
};

}
