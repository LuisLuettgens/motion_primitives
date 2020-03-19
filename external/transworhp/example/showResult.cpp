/*----------------------------------------------------------------
 *
 *  Example: Anzeigen einer *.m Datei
 *
 *----------------------------------------------------------------*/

#include "TWfolder.h"
#include "TWproblem.h"
#include "TransWORHP.h"

int main(int argv, char* argc[]) {

	tw::TWparameter twparameter("transworhp.xml");
	twparameter.Arguments(argv,argc);

	//twparameter.showDF = false;
	//twparameter.showDG = false;
	//twparameter.showHM = false;

	tw::TWfolder folder(&twparameter,0);
	tw::EmptyTransWorhpProblem tw(tw::TWdimension("showResult",twparameter.NDIS,9,2,1,0,0));
	tw.setSolver(&twparameter);
	tw.freetime = true;
	folder.Add(&tw);

	tw::Viewer* viewer = nullptr;
	if (twparameter.PLOT) viewer = new tw::Viewer(&twparameter);

	folder.Init(viewer);
	tw.solver->FromMATLAB("showResult.dat");

	folder.Show();

	delete viewer;

	return 0;
}
