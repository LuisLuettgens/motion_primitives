#pragma once

#include "twstatus.h"

#include <iostream>

namespace tw {

class TWfolder;
class TransWorhpProblem;

/** Ausgabe von Systemzustaenden und Debug-Informationen */
class TWdebug {
public:
	/** gibt die Eintraege der Ableitung der Zielfunktion aus. */
	static void PrintDF(TWfolder *twfolder, std::ostream *os=nullptr);
	/** gibt die Eintraege der Ableitung der (nicht trivialen) NB aus. */
	static void PrintDG(TWfolder *twfolder, std::ostream *os=nullptr);
	/** gibt die Eintraege der Hessematrix aus. */
	static void PrintHM(TWfolder *twfolder, std::ostream *os=nullptr);

	/** gibt alle Beschraenkungen aus. */
	static void PrintNLPconstraints(TWfolder *twfolder, std::ostream *os=nullptr);
	/** gibt die Beschraenkungen der Obimierungsvariablen aus. */
	static void PrintNLPconstraintsX(TWfolder *twfolder, TransWorhpProblem *tw, std::ostream *os=nullptr, int phase = -1);
	/** gibt die Beschraenkungen der NB aus. */
	static void PrintNLPconstraintsG(TWfolder *twfolder, TransWorhpProblem *tw, std::ostream *os=nullptr, int phase = -1);
	/** gibt die Beschraenkungen des TWfolder aus. */
	static void PrintNLPconstraintsG(TWfolder *twfolder, std::ostream *os=nullptr);

	/** Methode zum Ueberpruefen der NB. */
	static int CheckNLPconstraints(TWfolder *twfolder, std::ostream *os=nullptr);

private:
	/** Prueft, ob val zwischen low und upp liegt */
	static bool inequality(std::ostream &os, double low, double val, double upp);
	/** Prueft, ob val und ref (numerisch) gleich sind */
	static Status matrixentry(std::ostream& b, double val, double ref);
};

}
