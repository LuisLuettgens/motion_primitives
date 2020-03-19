/*----------------------------------------------------------------
*
* IAV Modellierungsseminar 2015: Prinzipienfreie Black-Box Modellierung eines Zylindermassenmodells
*
*----------------------------------------------------------------*/

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"
using namespace std;

int anz = 200; 

class ZylinderPhase : public TransWorhp {
public:

	ZylinderPhase(int dis) : TransWorhp("Zylindermassenmodell", dis, anz, 0, 0, 0, 0) {}

	double obj() override {
		return 1.0;
	}


	void ode(double *dx, double t, const double *x, const double *u, const double *p) override {
		for (int i = 0; i < anz; i++) {
			dx[i] = 1.0;
		}
	}
};

/////////////////////////////////////////////////////////////////////////////

int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	map<string, string> args = twparameter.Arguments(argv, argc);
	TWfolder folder(&twparameter, 0);

	ZylinderPhase ph(twparameter.NDIS);
	folder.Add(&ph);

	Viewer *viewer = nullptr;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);
	folder.Loop();

	delete viewer;

	return 0;
}
