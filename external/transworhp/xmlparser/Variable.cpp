/*----------------------------------------------------------------
*
* Klasse fuer Zustaende, Steuerungen und Parameter
*
*----------------------------------------------------------------*/

#include "Variable.h" 

// Konstruktor [0: Steuerung, 1: Zustand, 2: Beschraenkung]
Variable::Variable(int _varNum) {
	name = "";		
	varNum = _varNum;

	if (varNum != 0) {
		tree = new ::Node("");
	} else {
		tree = NULL;
	}
	
	start = nan("");
	fin = nan("");	
	lowBound = nan("");	
	uppBound = nan("");	
	initialGuess = nan("");	
}

// Ausgabe-Funktion
void Variable::print() {
	if (hasName()) { cout << "Name: " << name << endl; }
	if (hasTree()) { cout << "ODE: "; tree->print(); cout << endl; }
	if (hasStart()) { cout << "Startwert: " << start << endl; }
	if (hasFinal()) { cout << "Endwert: " << fin << endl; }
	if (hasUppBound() || hasLowBound()) { 
		cout << "Beschraenkungen von " << lowBound << " bis " << uppBound << endl; }
	if (hasGuess()) { cout << "Startschaetzung: " << initialGuess << endl; }
}

// Funktionen, mit denen auf Variableninhalt zugegriffen wird
string Variable::getName() {
	return name;
}
Node* Variable::getTree(){
	return (varNum != 0) ? tree : NULL;
}
double Variable::getStart(){
	return (varNum == 1) ? start : nan("");
}
double Variable::getFinal(){
	return (varNum == 1) ? fin : nan("");
}
double Variable::getUppBound(){
	return uppBound;
}
double Variable::getLowBound(){
	return lowBound;
}
double Variable::getGuess(){
	return initialGuess;
}

// Funktionen, mit denen Variableninhalt gesetzt wird
void Variable::setName(string _name){
	name = _name;
}
void Variable::setTree(Node *_tree){
	if (varNum != 0) {		
		tree = _tree;
	}
}
void Variable::setStart(double _start){
	if (varNum == 1) {		
		start = _start;
	}
}
void Variable::setFinal(double _final){
	if (varNum == 1) {	
		fin = _final;
	}
}
void Variable::setUppBound(double _uppBound){
	uppBound = _uppBound;
}
void Variable::setLowBound(double _lowBound){
	lowBound = _lowBound;
}
void Variable::setGuess(double _guess){
	initialGuess = _guess;
}

// Funktionen, die angeben ob Variableninhalt gesetzt ist
bool Variable::hasName(){
	return name.compare("");
}

bool Variable::hasTree(){
	return (tree != NULL) && (tree->getValue()).compare("");
}

bool Variable::hasStart(){
	return !isnan(start);
}
bool Variable::hasFinal(){
	return !isnan(fin);
}
bool Variable::hasUppBound(){
	return !isnan(uppBound);
}
bool Variable::hasLowBound(){
	return !isnan(lowBound);
}
bool Variable::hasGuess(){
	return !isnan(initialGuess);
}