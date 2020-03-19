/*----------------------------------------------------------------
*
* Klasse fuer Baumknoten
*
*----------------------------------------------------------------*/

#include "Node.h" 

// Konstruktoren
Node::Node(string _value) {
	value = _value;
	valueType = EMPTY;
        leftChild = NULL;
        rightChild = NULL;
}

Node::Node(string _value,  Node *left, Node *right) {
	value = _value;
	valueType = EMPTY;
	leftChild = left;
	rightChild = right;
}

// Ausgabe-Funktion zum Testen
void Node::print() {
	if (hasLeftChild() && hasRightChild()) {
        	cout << "("; 
		leftChild->print(); 
		cout << ") " << value << " ("; 
		rightChild->print(); 
		cout << ")";
	} else {
        	cout << value;
	}
}

// Funktionen, mit denen auf Knoteninhalt zugegriffen wird
string Node::getValue() {
	return value;
}

ValueType Node::getValueType() {
	return valueType;
}

double Node::getDoubleValue() {
	return doubleValue;
}

OperationType Node::getOperationValue() {
	return operationType;
}

FunctionType Node::getFunctionValue() {
	return functionType;
}

int Node::getNameValueIndex() {
	return nameIndex;
}

NameType Node::getNameValueType() {
	return nameType;
}

NameTime Node::getNameTime() {
	return nameTime;
}

Node* Node::getLeftChild() {
	return leftChild;
}

Node* Node::getRightChild() {
	return rightChild;
}

// Funktionen, mit denen der Knoteninhalt gesetzt wird
void Node::setValue(string _value) {
	value = _value;
}

void Node::setLeftChild(Node *_leftChild) {
	leftChild = _leftChild;
}

void Node::setRightChild(Node *_rightChild) {
	rightChild = _rightChild;
}

// Funktionen, die angeben, ob Knoteninhalt gesetzt ist 
bool Node::hasLeftChild(){
	return leftChild != NULL;
}

bool Node::hasRightChild(){
	return rightChild != NULL;	
}

bool Node::isLeaf(){
	return (leftChild == NULL) && (rightChild == NULL);
}

// Funktion: parst Text in Baum
void Node::buildTree(string* stateNames, string* controlNames, string* paramNames, int stateCount, int controlCount, int paramCount){
	// falls Klammern aussen um Term rum -> Klammern rausloeschen
	if (value[0] == '(' && value[value.size()-1] == ')') {
		int brackets = 0;
		bool change = true;
		size_t i = value.find_first_of("()");

		while (i != string::npos) {
			if (value[i] == '(') { brackets = brackets+1; }
			if (value[i] == ')') { brackets = brackets-1; }

			if (brackets == 0 && i != value.size()-1) { change = false; }

			i = value.find_first_of("()",i+1);
		}
		
		if (change) {
			value = value.substr(1,value.size()-2);
		}
	}

	// bei einer Negation einer Variable wird der Operator als Minus behandelt und im linken Kind
	// ein leerer String erzeugt
	if (value == "") {
		value = "0";
	}

	// Hilfsfunktion, gibt Position eines Operators oder einer Funktion zurueck	
	size_t i = getNextNodePos(value);

	if (i != string::npos) {
		if (value[i] == 's' || value[i] == 'c' || value[i] == 't' || value[i] == 'e') {
			// valueType ist eine Funktion
			valueType = FUNCTION;

			// linkes Kind wird nie benoetigt, Funktionsargument ist im rechten Kind
			Node *rightCh = new Node(value.substr(i+3));
			value = value.substr(i,3);

			if (value.compare("sin") == 0) {
				functionType = SIN;
			} else if (value.compare("cos") == 0) {
				functionType = COS;
			} else if (value.compare("tan") == 0) {
				functionType = TAN;
			} else if (value.compare("exp") == 0) {
				functionType = EXP;
			} else {
				cout << "Error: " << value << " is no known function." << endl;
				exit(EXIT_FAILURE);
			}

			setRightChild(rightCh);
			rightCh->buildTree(stateNames, controlNames, paramNames, stateCount, controlCount, paramCount);

			return;
		} else {
			// valueType ist eine Operation
			valueType = OPERATION;

			Node *leftCh = new Node(value.substr(0,i));
			Node *rightCh = new Node(value.substr(i+1));
			value = value.substr(i,1);

			if (value[0] == '+') {
				operationType = ADD;
			} else if (value[0] == '-') {
				operationType = SUBTRACT;
			} else if (value[0] == '*') {
				operationType = MULTIPLY;
			} else if (value[0] == '/') {
				operationType = DIVIDE;
			} else if (value[0] == '^') {
				operationType = EXPONENTIATION;
			} else {
				cout << "Error: " << value[0] << " is no known operation." << endl;
				exit(EXIT_FAILURE);		
			}

			setLeftChild(leftCh);
			setRightChild(rightCh);
	
			leftCh->buildTree(stateNames, controlNames, paramNames, stateCount, controlCount, paramCount);
			rightCh->buildTree(stateNames, controlNames, paramNames, stateCount, controlCount, paramCount);

			return;
		}
	} else if (strtod(value.c_str(),NULL) != 0.0 || value == "0") {
		// valueType ist eine Zahl
		valueType = DOUBLE;
		doubleValue = atof(value.data());
	} else {
		// valueType ist eine Variable
		valueType = NAME;
		
		size_t i = value.find_first_of("(");
		if (i != string::npos) {
			size_t j = value.find_first_of(")");
			if (j == string::npos) {
				cout << "Error: " << value << " no ')' found" << endl;
				exit (EXIT_FAILURE);
			}

			string tmp = value.substr(i+1,j-i-1);	// Klammerinhalt
			string val = value.substr(0,i);		// Variable vor Klammer

			for (int j=0; j<stateCount; j++) {
				if (!val.compare(stateNames[j])) {
					if (tmp.compare("end")==0.0) {
						nameType = STATE;
						nameIndex = j;
						nameTime = LAST;						

						return;
					} else if (tmp.compare("0")==0.0) {
						nameType = STATE;
						nameIndex = j;
						nameTime = FIRST;						

						return;
					}
				}
			}
			for (int j=0; j<controlCount; j++) {
				if (!val.compare(controlNames[j])) {
					if (tmp.compare("end")==0.0) {
						nameType = CONTROL;
						nameIndex = j;
						nameTime = LAST;						

						return;
					} else if (tmp.compare("0")==0.0) {
						nameType = CONTROL;
						nameIndex = j;
						nameTime = FIRST;						

						return;
					}
				}
			}
		} else {
			for (int j=0; j<stateCount; j++) {
				if (!value.compare(stateNames[j])) {
					nameType = STATE;
					nameIndex = j;
					nameTime = NONE;						

					return;
				}
			}

			for (int j=0; j<controlCount; j++) {
				if (!value.compare(controlNames[j])) {
					nameType = CONTROL;
					nameIndex = j;
					nameTime = NONE;						

					return;
				}
			}

			for (int j=0; j<paramCount; j++) {
				if (!value.compare(paramNames[j])) {
					nameType = PARAM;
					nameIndex = j;
					nameTime = NONE;						

					return;
				}
			}
		}

		valueType = EMPTY;		
	}
}

// Funktion: Hilfsfunktion fuer buildTree, gibt Operatorposition zurueck
size_t Node::getNextNodePos(string text){
	int brackets = 0;

	// Suche nach '+'
	size_t i = text.find_first_of("+()");
	while (i != string::npos) {
		if (text[i] == '+' && brackets == 0) { return i; }

		if (text[i] == '(') { brackets = brackets+1; }
		if (text[i] == ')') { brackets = brackets-1; }

		i = text.find_first_of("+()",i+1);
	}

	// Suche nach '-'
	i = text.find_last_of("-()");
	while (i != string::npos) {
		if (text[i] == '-' && brackets == 0) { return i; }

		if (text[i] == '(') { brackets = brackets+1; }
		if (text[i] == ')') { brackets = brackets-1; }
		
		if ( i == 0) {
			break;
		} else {
			i = text.find_last_of("-()",i-1);
		}
	}
	
	// Suche nach '*'
	i = text.find_first_of("*()");
	while (i != string::npos) {
		if (text[i] == '*' && brackets == 0) { return i; }

		if (text[i] == '(') { brackets = brackets+1; }
		if (text[i] == ')') { brackets = brackets-1; }

		i = text.find_first_of("*()",i+1);
	}

	// Suche nach '/'
	i = text.find_last_of("/()");
	while (i != string::npos) {
		if (text[i] == '/' && brackets == 0) { return i; }

		if (text[i] == '(') { brackets = brackets+1; }
		if (text[i] == ')') { brackets = brackets-1; }

		if ( i == 0) {
			break;
		} else {
			i = text.find_last_of("/()",i-1);
		}
	}

	// Suche nach '^'
	i = text.find_last_of("^()");
	while (i != string::npos) {
		if (text[i] == '^' && brackets == 0) { return i; }

		if (text[i] == '(') { brackets = brackets+1; }
		if (text[i] == ')') { brackets = brackets-1; }

		if ( i == 0) {
			break;
		} else {
			i = text.find_last_of("^()",i-1);
		}
	}

	// Suche nach 'sin'
	i = text.find_first_of("s()");
	while (i != string::npos) {
		if ((text.substr(i,3)).compare("sin") == 0 && brackets == 0) { return i; }

		if (text[i] == '(') { brackets = brackets+1; }
		if (text[i] == ')') { brackets = brackets-1; }

		i = text.find_first_of("s()",i+1);
	}

	// Suche nach 'cos'
	i = text.find_first_of("c()");
	while (i != string::npos) {
		if ((text.substr(i,3)).compare("cos") == 0 && brackets == 0) { return i; }

		if (text[i] == '(') { brackets = brackets+1; }
		if (text[i] == ')') { brackets = brackets-1; }

		i = text.find_first_of("c()",i+1);
	}

	// Suche nach 'tan'
	i = text.find_first_of("t()");
	while (i != string::npos) {
		if ((text.substr(i,3)).compare("tan") == 0 && brackets == 0) { return i; }

		if (text[i] == '(') { brackets = brackets+1; }
		if (text[i] == ')') { brackets = brackets-1; }

		i = text.find_first_of("t()",i+1);
	}

	// Suche nach 'exp'
	i = text.find_first_of("e()");
	while (i != string::npos) {
		if ((text.substr(i,3)).compare("exp") == 0 && brackets == 0) { return i; }

		if (text[i] == '(') { brackets = brackets+1; }
		if (text[i] == ')') { brackets = brackets-1; }

		i = text.find_first_of("e()",i+1);
	}

	return string::npos;
}