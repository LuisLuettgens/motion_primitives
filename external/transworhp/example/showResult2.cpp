/*----------------------------------------------------------------
 *
 *  Example: Anzeigen einer *.m Datei
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

#include "laufkatze_base.h"

double phasenPlot (int&, int &n_ode, int &n_ctrl, double*,double *x, double *, int &i, int &index) {
	if (index==0)		return x[i*(n_ode+n_ctrl)+5];
	else if (index==1)	return x[i*(n_ode+n_ctrl)];
	else return 0.0;
}

class MyEmptyTransWorhp : public tw::EmptyTransWorhpProblem {

public:
	std::unique_ptr<XMLNode> xmlmain;
	XMLNode* scenenode;

	MyEmptyTransWorhp(int dis) : EmptyTransWorhpProblem(tw::TWdimension("showResult",dis,9,2,1,0,0)),
		xmlmain(tw::TWparameter::ReadParams("laufkatze.xml")), scenenode(nullptr) {

		freetime = true;

		if (xmlmain) {
			scenenode = xmlmain->GetFirstChild("SCENE");
		}
	}

	void OpenWindows(tw::Viewer *viewer) override {
		viewer->PhasePlot("Phasen-Plot", phasenPlot,0,1);
		viewer->ThreeD("Hallenansicht", scenenode, laufkatze3d);
	}

	void selectWindows(tw::Viewer *viewer) override {

		viewer->AddStateView(0,"Position X");
		viewer->AddStateView(1,"Geschwindigkeit X");
		viewer->AddStateView(2,"Verschiebung Last");
		viewer->AddStateView(3,"Geschwindigkeit Last");
		viewer->AddStateView(4,"Beschleunigung X");
		viewer->AddStateView(5,"Position Y");
		viewer->AddStateView(6,"Geschwindigkeit Y");
		viewer->AddStateView(7,"Beschleunigung Y");
		viewer->AddStateView(8,"u^2");
		
		viewer->AddControlView(0,"Ruck X");
		viewer->AddControlView(1,"Ruck Y");
	}
};

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);

	twparameter.showDF = false;
	twparameter.showDG = false;
	twparameter.showHM = false;

	tw::TWfolder folder(&twparameter,0);

	MyEmptyTransWorhp tw(twparameter.NDIS);
	tw.setSolver(&twparameter);
	folder.Add(&tw);

	tw::Viewer* viewer = nullptr;
	if (twparameter.PLOT) viewer = new tw::Viewer(&twparameter);
	
	folder.Init(viewer);
	tw.solver->FromMATLAB("showResult.dat");

	folder.Show();
	
	delete viewer;

	return 0;
}
