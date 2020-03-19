/*----------------------------------------------------------------
 *
 * Header-Datei fuer Zustaende, Steuerungen und Parameter
 *
 *----------------------------------------------------------------*/

#include "Node.h"

using namespace std;

class Variable {
private:
	string name;			// Bezeichnung in der XML-Datei
	::Node *tree;			// Term
	int varNum;			// 0: Steuerung, 1: Zustand, 2: Beschraenkung
	double start;			// Anfangswert, nur falls Zustand
	double fin;			// Endwert, nur falls Zustand
	double lowBound;		// untere Begraenzung
	double uppBound;		// obere Begraenzung
	double initialGuess;		// Startschaetzung

public:
	// Konstruktor	
	Variable(int varNum);

	// Ausgabe-Funktion
	void print();

	// Funktionen, mit denen auf Variableninhalt zugegriffen wird
	string getName();
	Node* getTree();
	double getStart();
	double getFinal();
	double getUppBound();
	double getLowBound();
	double getGuess();

	// Funktionen, mit denen Variableninhalt gesetzt wird
	void setName(string _name);
	void setTree(Node *_tree);
	void setStart(double _start);
	void setFinal(double _final);
	void setUppBound(double _uppBound);
	void setLowBound(double _lowBound);
	void setGuess(double _guess);

	// Funktionen, die angeben ob Variableninhalt gesetzt ist
	bool hasName();
	bool hasTree();
	bool hasStart();
	bool hasFinal();
	bool hasUppBound();
	bool hasLowBound();
	bool hasGuess();
};
