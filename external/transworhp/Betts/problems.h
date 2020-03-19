extern "C"
{
	void aquade_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void aquain_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void araode_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void araoin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void alprde_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void alprin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void brgrde_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void brgrin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void ashrde_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void ashrin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void bracde_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void bracin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void chmrde_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void chmrin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void crande_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void cranin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void ffrbde_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void ffrbin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void hangde_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void hangin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void heatde_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void heatin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void jshide_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void jshiin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void lnhtde_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void lnhtin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void lntsde_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void lntsin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void lowtde_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void lowtin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void kplrde_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void kplrin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void medide_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void mediin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void pndlde_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void pndlin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void rbrmde_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void rbrmin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void tb2sde_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void tb2sin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void trande_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void tranin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
	void zrmlde_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void zrmlin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	
/*	void bangde_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void bangin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
	*/
	void qlinde_(	int &iphase,double &time,const double *yvec,int &nyvec,const double *p,
			const int &np,double *frhs,const int &nfrhs,int &iferr);
	void qlinin_(	int &IPHASE,int &NPHS,int &METHOD,int &NSTG,int *NCF,
			int *NPF,int &NPV,int &NAV,int &NGRID,int &INIT,
			int &MAXMIN,const int &MXPARM,double *P0,double *PLB,double *PUB,
			char *plbl,const int &MXSTAT,double *Y0,double *Y1,double *YLB,
			double *YUB,double *STSKL,char *stlbl,const int &MXPCON,double *CLB,
			double *CUB,char *clbl,const int &MXTERM,double *COEF,int *ITERM,
			char *TITLE,int &IER);
}
