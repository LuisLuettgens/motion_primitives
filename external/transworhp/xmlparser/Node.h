/*----------------------------------------------------------------
 *
 * Header-Datei fuer Node.cpp
 *
 *----------------------------------------------------------------*/

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>

using namespace std;

// Aufzaehlungen zum Speichern des interpretierten Knoteninhalts
enum ValueType {EMPTY, DOUBLE, NAME, OPERATION, FUNCTION};
enum NameType {STATE, CONTROL, PARAM};
enum NameTime {NONE, FIRST, LAST};
enum OperationType {ADD, SUBTRACT, MULTIPLY, DIVIDE, EXPONENTIATION};
enum FunctionType {SIN, COS, TAN, EXP};

class Node {

private:
	string value;
	ValueType valueType;

	double doubleValue;				// DOUBLE: 		Zahlenwert
	OperationType operationType;	// OPERATION: 	Operator
	FunctionType functionType;		// FUNCTION:	verwendete Funktion
	int nameIndex;					// NAME: 		Index-Nummer der Variable
	NameType nameType;				
	NameTime nameTime;				// NAME:		Zeitangabe der Variable

	Node *leftChild;
	Node *rightChild;

public:
	// Konstruktoren	
	Node(string _value);
	Node(string _value, Node *left, Node *right);

	// Ausgabe-Funktion
	void print();

	// Funktionen, mit denen auf Knoteninhalt zugegriffen wird
	string getValue();
	ValueType getValueType();
	double getDoubleValue();
	OperationType getOperationValue();
	FunctionType getFunctionValue();
	int getNameValueIndex();
        NameType getNameValueType();
	NameTime getNameTime();

	Node* getLeftChild();
	Node* getRightChild();

	// Funktionen, mit denen der Knoteninhalt gesetzt wird
	void setValue(string _value); // darf nicht mehr nach buildTree aufgerufen werden!
	void setLeftChild(Node *_leftChild);
	void setRightChild(Node *_rightChild);

	// Funktionen, die angeben, ob Knoteninhalt gesetzt ist 
	bool hasLeftChild();
	bool hasRightChild();
	bool isLeaf();

	// Funktionen, die Baum erstellen und interpretieren
	void buildTree(string* stateNames, string* controlNames, string* paramNames, int stateCount, int controlCount, int paramCout);
	size_t getNextNodePos(string text);
};