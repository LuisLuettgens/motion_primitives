/*----------------------------------------------------------------
*
* Klasse fuer XMLTransWorhp:
* ruft mit eingelesenen Daten TransWORHP Funktionen auf
*
*----------------------------------------------------------------*/

#include "XMLTransWorhp.h" 

// Konstruktor
XMLTransWorhp::XMLTransWorhp(int dis, vector<Variable*> _states, vector<Variable*> _controls, vector<Variable*> _params, Node *_objAnchor, vector<Variable*> _randc, vector<Variable*> _nebenc) : TransWorhp("XMlTransWorhp",dis,_states.size(), _controls.size(),_params.size(),_randc.size(),_nebenc.size()) {
	states = _states; 
	controls = _controls; 
	params = _params; 

	randc = _randc;
	nebenc = _nebenc;
	objAnchor = _objAnchor;

	stateCount = states.size();
	controlCount = controls.size();
	paramCount = params.size();
	randCount = randc.size();
	nebenCount = nebenc.size();

	// Array fuer Rechentermstrukturen erstellen
	diffArray = new int*[stateCount+randCount+nebenCount+1];
	for (vector<Variable*>::size_type i = 0; i<stateCount+randCount+nebenCount+1; i++) {
		diffArray[i] = new int[stateCount+controlCount+paramCount];
		for (int j=0; j<stateCount+controlCount+paramCount; j++) {
			diffArray[i][j] = 0;
		}
	}

	// Array fuellen:
	for (vector<Variable*>::size_type i=0; i<stateCount+randCount+nebenCount+1; i++) {
		if (i < stateCount) {
			// ODE
			makeStructure(states[i]->getTree(), i);
		} else if (i < stateCount+randCount) {
			// Randbedingung
			makeStructure(randc[i-stateCount]->getTree(), i);
		} else if (i < stateCount+randCount+nebenCount) {
			// Nebenbedingung
			makeStructure(nebenc[i-stateCount-randCount]->getTree(), i);
		} else {
			// Zielfunktion
			makeStructure(objAnchor, i);
		}
	}
}

// Hilfsfunktion stellt diffArray auf
void XMLTransWorhp::makeStructure(Node *n, int num) {
	if (n->hasLeftChild()) {
		makeStructure(n->getLeftChild(), num);
	}

	if (n->hasRightChild()) {
		makeStructure(n->getRightChild(), num);
	}

	if (!n->hasLeftChild() && !n->hasRightChild()) {
		ValueType valType = n->getValueType();

		// falls keine Variable wird nichts gesetzt
		if (valType == NAME) {
			NameTime time = n->getNameTime();

			if (time == NONE) {
				switch (n->getNameValueType()) {
					case STATE: 	diffArray[num][n->getNameValueIndex()] = 1;	return;
					case CONTROL: 	diffArray[num][n->getNameValueIndex()+stateCount] = 1; return;
					case PARAM: 	diffArray[num][n->getNameValueIndex()+stateCount+controlCount] = 1;	return;
				}
			} else {
				if (time == FIRST) {
					switch (n->getNameValueType()) {
						case STATE: 	diffArray[num][n->getNameValueIndex()] = 2;	return;
						case CONTROL: 	diffArray[num][n->getNameValueIndex()+stateCount] = 2; return;
						case PARAM: 	diffArray[num][n->getNameValueIndex()+stateCount+controlCount] = 2;	return;
					}
				} else if (time == LAST) {
					switch (n->getNameValueType()) {
						case STATE: 	diffArray[num][n->getNameValueIndex()] = 3;	return;
						case CONTROL: 	diffArray[num][n->getNameValueIndex()+stateCount] = 3; return;
						case PARAM: 	diffArray[num][n->getNameValueIndex()+stateCount+controlCount] = 3;	return;
					}
				}
			}
		}
	}
}

// Zielfunktion
double XMLTransWorhp::obj() {
	return evaluateObj(objAnchor);
}

// Funktion: liest den Zielfunktionsbaum aus 
double XMLTransWorhp::evaluateObj(Node *n) {
	ValueType val = n->getValueType();
	
	if (val == OPERATION) {
		double valL = evaluateObj(n->getLeftChild());
	        double valR = evaluateObj(n->getRightChild());

		switch (n->getOperationValue()) {
			case ADD: 		return valL + valR;
			case SUBTRACT: 		return valL - valR;
			case MULTIPLY: 		return valL * valR;
			case DIVIDE: 		return valL / valR;
			case EXPONENTIATION: 	return pow(valL, valR);
		}
	} else if (val == FUNCTION) {
	        double valR = evaluateObj(n->getRightChild());

		switch (n->getFunctionValue()) {
			case SIN: return sin (valR);
			case COS: return cos (valR);
			case TAN: return tan (valR);
			case EXP: return exp (valR);
		}		
	} else if (val == NAME) {
		if (n->getNameValueType() == STATE) {
			if (n->getNameTime() == FIRST) {
				return x(0,n->getNameValueIndex());
			} else if (n->getNameTime() == LAST) {
				return x(n_dis-1, n->getNameValueIndex());
			}
		} else if (n->getNameValueType() == CONTROL) {
			if (n->getNameTime() == FIRST) {
				return u(0,n->getNameValueIndex());
			} else if (n->getNameTime() == LAST) {
				return u(n_dis-1, n->getNameValueIndex());
			}
		} else if (n->getNameValueType() == PARAM) {
			return p(n->getNameValueIndex());
		}	 
	} else if (val == DOUBLE) {
		return n->getDoubleValue();
	} else {
		cout << "Error: Objective or Rand: '" << val << "' is no parameter" << endl;
		exit (EXIT_FAILURE);
		return -1;
	}
}

double XMLTransWorhp::evaluateDiffObj(Node *n, int diff_param) {
	ValueType val = n->getValueType();

	if (val == OPERATION) {
		double diff_valL = evaluateDiffObj(n->getLeftChild(),diff_param);
		double diff_valR = evaluateDiffObj(n->getRightChild(),diff_param);

		double valL = evaluateObj(n->getLeftChild());
		double valR = evaluateObj(n->getRightChild());
	
		switch (n->getOperationValue()) {
			case ADD: 		return diff_valL + diff_valR;
			case SUBTRACT: 		return diff_valL - diff_valR;
			case MULTIPLY: 		return diff_valL*valR + diff_valR*valL;
			case DIVIDE: 		return (diff_valL*valR-diff_valR*valL)/(valR*valR);
			case EXPONENTIATION: 	return valR*pow(valL, valR-1)*diff_valL;
		}
	} else if (val == FUNCTION) {
		double diff_valR = evaluateDiffObj(n->getRightChild(),diff_param);
	        double valR = evaluateObj(n->getRightChild());

		switch (n->getFunctionValue()) {
			case SIN: return diff_valR * cos (valR);
			case COS: return - diff_valR * sin (valR);
			case TAN: return diff_valR / (cos (valR) * cos (valR));
			case EXP: return diff_valR * exp (valR);
		}		
	} else if (val == NAME) {
		if (diff_param < stateCount) {
			if (n->getNameValueType() == STATE && n->getNameValueIndex() == diff_param) {
				return 1;
			}
			return 0;
		} else if (diff_param>=stateCount && diff_param<stateCount+controlCount) {
			if (n->getNameValueType() == CONTROL && n->getNameValueIndex() == diff_param-stateCount) {
				return 1;
			}
			return 0;
		} else if (diff_param >= stateCount+controlCount) {
			if (n->getNameValueType() == PARAM && n->getNameValueIndex() == diff_param-stateCount-controlCount) {
				return 1;
			}
			return 0;
		}
	} else if (val == DOUBLE) {
		return 0;
	} else { 
		cout << "Error: Objective: '" << n->getValue() << "' is no parameter" << endl;
		exit (EXIT_FAILURE);
		return -1;
	}
}

bool XMLTransWorhp::obj_structure(DiffStructure &s) {
	for (int j=0; j<stateCount+controlCount+paramCount; j++) {
		if (diffArray[stateCount+randCount+nebenCount][j] == 1) {
			s(0,p_index(j-stateCount-controlCount));	
		} else if (diffArray[stateCount+randCount+nebenCount][j] == 3) {
			if (j<stateCount) {
				s(0,x_index(n_dis-1,j));
			} else if(j>=stateCount && j<stateCount+controlCount) {
				s(0,u_index(n_dis-1,j-stateCount));
			}
		}
	}

	return true;
}

bool XMLTransWorhp::obj_diff(DiffStructure &s) {
	for (int j=0; j<stateCount+controlCount+paramCount; j++) {
		if (diffArray[stateCount+randCount+nebenCount][j] == 1) {
			s(0,p_index(j-stateCount-controlCount)) = evaluateDiffObj(objAnchor,j);
		} else if (diffArray[stateCount+randCount+nebenCount][j] == 3) {
			if (j<stateCount) {
				s(0,x_index(n_dis-1,j)) = evaluateDiffObj(objAnchor,j);
			} else if(j>=stateCount && j<stateCount+controlCount) {
				s(0,u_index(n_dis-1,j-stateCount)) = evaluateDiffObj(objAnchor,j);
			}
		}
	}

	return true;
}

// ODE-System
void XMLTransWorhp::ode(double *dx, double t, const double *x, const double *u, const double *p) {
	for (vector<Variable*>::size_type i=0; i<stateCount; i++) {
		if (states[i]->hasTree()) {
			dx[i] = evaluateOde(states[i]->getTree(), x, u, p);
		}
	}
}

// Funktion: liest den ODE-Baum aus
double XMLTransWorhp::evaluateOde(Node *n, const double *x, const double *u, const double *p) {
	ValueType val = n->getValueType();
	
	if (val == OPERATION) {
		double valL = evaluateOde(n->getLeftChild(),x,u,p);
	        double valR = evaluateOde(n->getRightChild(),x,u,p);

		switch (n->getOperationValue()) {
			case ADD: 				return valL + valR;
			case SUBTRACT: 			return valL - valR;
			case MULTIPLY: 			return valL * valR;
			case DIVIDE: 			return valL / valR;
			case EXPONENTIATION: 	return pow(valL, valR);
		}
	} else if (val == FUNCTION) {
	        double valR = evaluateOde(n->getRightChild(),x,u,p);

		switch (n->getFunctionValue()) {
			case SIN: return sin (valR);
			case COS: return cos (valR);
			case TAN: return tan (valR);
			case EXP: return exp (valR);
		}		
	} else if (val == NAME) {
		switch (n->getNameValueType()) {
			case STATE: 	return x[n->getNameValueIndex()];
			case CONTROL: 	return u[n->getNameValueIndex()];
			case PARAM: 	return p[n->getNameValueIndex()];
		}	 
	} else if (val == DOUBLE) {
		return n->getDoubleValue();
	} else {
		cout << "Error: ODE or Neben: '" << n->getValue() << "' is no parameter" << endl;
		exit (EXIT_FAILURE);
		return -1;
	}
}

double XMLTransWorhp::evaluateDiffOde(Node *n, const double *x, const double *u, const double *p, int diff_param) {
	ValueType val = n->getValueType();

	if (val == OPERATION) {
		double diff_valL = evaluateDiffOde(n->getLeftChild(),x,u,p,diff_param);
		double diff_valR = evaluateDiffOde(n->getRightChild(),x,u,p,diff_param);

		double valL = evaluateOde(n->getLeftChild(),x,u,p);
		double valR = evaluateOde(n->getRightChild(),x,u,p);
	
		switch (n->getOperationValue()) {
			case ADD: 				return diff_valL + diff_valR;
			case SUBTRACT: 			return diff_valL - diff_valR;
			case MULTIPLY: 			return diff_valL*valR + diff_valR*valL;
			case DIVIDE: 			return (diff_valL*valR-diff_valR*valL)/(valR*valR);
			case EXPONENTIATION: 	return valR*pow(valL, valR-1)*diff_valL;
		}
	} else if (val == FUNCTION) {
		double diff_valR = evaluateDiffOde(n->getRightChild(),x,u,p,diff_param);
	        double valR = evaluateOde(n->getRightChild(),x,u,p);

		switch (n->getFunctionValue()) {
			case SIN: return diff_valR * cos (valR);
			case COS: return - diff_valR * sin (valR);
			case TAN: return diff_valR / (cos (valR) * cos (valR));
			case EXP: return diff_valR * exp (valR);
		}
	} else if (val == NAME) {
		if (diff_param < stateCount) {
			if (n->getNameValueType() == STATE && n->getNameValueIndex() == diff_param) {
				return 1;
			}
			return 0;
		} else if (diff_param>=stateCount && diff_param<stateCount+controlCount) {
			if (n->getNameValueType() == CONTROL && n->getNameValueIndex() == diff_param-stateCount) {
				return 1;
			}
			return 0;
		} else if (diff_param >= stateCount+controlCount) {
			if (n->getNameValueType() == PARAM && n->getNameValueIndex() == diff_param-stateCount-controlCount) {
				return 1;
			}
			return 0;
		}
	} else if (val == DOUBLE) {
		return 0;
	} else { 
		cout << "Error: ODE: '" << n->getValue() << "' is no parameter" << endl;
		exit (EXIT_FAILURE);
		return -1;
	}
}

bool XMLTransWorhp::ode_structure(DiffStructure &s) {
	for (int i=0; i<stateCount; i++) {
		for (int j=0; j<stateCount+controlCount+paramCount; j++) {
			if (diffArray[i][j] == 1) {
				if (j<stateCount) {
					s(i,x_index(0,j));
				} else if(j>=stateCount && j<stateCount+controlCount) {
					s(i,u_index(0,j-stateCount));
				} else if(j>=stateCount+controlCount) {
					s(i,p_indexode(j-stateCount-controlCount));
				}
			}
		}
	}

	return true;
}

bool XMLTransWorhp::ode_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) {
	for (int i=0; i<stateCount; i++) {
		for (int j=0; j<stateCount+controlCount; j++) {
			if (diffArray[i][j] == 1) {
				if (j < stateCount) {
					s(i,x_indexode(j)) = evaluateDiffOde(states[i]->getTree(),x,u,p,j);
				} else if ( j>=stateCount ) {
					s(i,u_indexode(j-stateCount)) = evaluateDiffOde(states[i]->getTree(),x,u,p,j);
				}
			}
		}
	}
	return true;
}

bool XMLTransWorhp::ode_diff_p(DiffStructure &s, double t, const double *x, const double *u, const double *p, int index) {
	for (int i=0; i<stateCount; i++) {
		for (int j=stateCount+controlCount; j<stateCount+controlCount+paramCount; j++) {
			if (diffArray[i][j] == 1) {				
				s(i,p_indexode(j-stateCount-controlCount)) = evaluateDiffOde(states[i]->getTree(),x,u,p,j);
			}
		}
	}
	return false;
}


// Box Beschraenkungen	
void XMLTransWorhp::x_boundary(double *x_low, double *x_upp) {
	for (vector<Variable*>::size_type i=0; i<stateCount; i++) {
		if (states[i]->hasLowBound()) { x_low[i] = states[i]->getLowBound();}
		if (states[i]->hasUppBound()) { x_upp[i] = states[i]->getUppBound();}
	}
}

void XMLTransWorhp::u_boundary(double *u_low, double *u_upp) {
	for (vector<Variable*>::size_type i=0; i<controlCount; i++) {
		if (controls[i]->hasLowBound()) { u_low[i] = controls[i]->getLowBound();}
		if (controls[i]->hasUppBound()) { u_upp[i] = controls[i]->getUppBound();}
	}
}

void XMLTransWorhp::p_boundary(double *p_low, double *p_upp) {
	for (vector<Variable*>::size_type i=0; i<paramCount; i++) {
		if (params[i]->hasLowBound()) { p_low[i] = params[i]->getLowBound();}
		if (params[i]->hasUppBound()) { p_upp[i] = params[i]->getUppBound();}
	}
}
	
// Anfangs- und Endwerte
void XMLTransWorhp::var_boundary(double *x_low, double *x_upp) {
	for (vector<Variable*>::size_type i=0; i<stateCount; i++) {			
		if (states[i]->hasStart()) { 
			x_low[x_index(0,i)] = x_upp[x_index(0,i)] = states[i]->getStart();
		}
		if (states[i]->hasFinal()) { 
			x_low[x_index(n_dis-1,i)] = x_upp[x_index(n_dis-1,i)] = states[i]->getFinal();
		}
	}
}

// Randwerte
void XMLTransWorhp::rand(double *r) {
	for (vector<Variable*>::size_type i=0; i<randc.size(); i++) {
		if (randc[i]->hasTree()) {
			r[i] = evaluateObj(randc[i]->getTree());
		}
	}
}

void XMLTransWorhp::rand_boundary(double *r_low, double *r_upp) {
	for (vector<Variable*>::size_type i=0; i<randc.size(); i++) {
		if (randc[i]->hasLowBound()) { r_low[i] = randc[i]->getLowBound();}
		if (randc[i]->hasUppBound()) { r_upp[i] = randc[i]->getUppBound();}
	}
}

bool XMLTransWorhp::rand_structure(DiffStructure &s) {	
	if (randCount == 0) {
		return false;
	}

	int b = stateCount;
	int c = controlCount;

	for (int i=0; i<randc.size(); i++) {
		for (int j=0; j<stateCount+controlCount+paramCount; j++) {
			if (diffArray[i+stateCount][j] == 1 && j>=b+c) {
				s(0,p_index(j-b-c));
			} else if (diffArray[i+stateCount][j] == 2 && j<b+c) {
				if (j<b) {
					s(0,x_index(0,j));
				} else {
					s(0,u_index(0,j-b));
				}
			} else if (diffArray[i+stateCount][j] == 3 && j<b+c) {
				if (j<b) {
					s(0,x_index(n_dis-1,j));
				} else {
					s(0,u_index(n_dis-1,j-b));
				}
			}
		}
	}

	return true;
}

bool XMLTransWorhp::rand_diff(DiffStructure &s){
	if (randCount == 0) {
		return false;
	}

	int b = stateCount;
	int c = controlCount;

	for (int i=0; i<randc.size(); i++) {
		for (int j=0; j<stateCount+controlCount+paramCount; j++) {
			if (diffArray[i+stateCount][j] == 1 && j>=b+c) {
				s(0,p_index(j-b-c)) = evaluateDiffObj(objAnchor,j);
			} else if (diffArray[i+stateCount][j] == 2 && j<b+c) {
				if (j<b) {
					s(0,x_index(0,j))  = evaluateDiffObj(objAnchor,j);
				} else {
					s(0,u_index(0,j-b)) = evaluateDiffObj(objAnchor,j);
				}
			} else if (diffArray[i+stateCount][j] == 3 && j<b+c) {
				if (j<b) {
					s(0,x_index(n_dis-1,j)) = evaluateDiffObj(objAnchor,j);
				} else {
					s(0,u_index(n_dis-1,j-b)) = evaluateDiffObj(objAnchor,j);
				}
			}
		}
	}

	return true;
}

// Nebenbedingungen
void XMLTransWorhp::neben(double *c, double t, const double *x, const double *u, const double *p) {
	for (vector<Variable*>::size_type i=0; i<nebenc.size(); i++) {
		if (nebenc[i]->hasTree()) {
			c[i] = evaluateOde(nebenc[i]->getTree(), x, u, p);
		}
	}
}

void XMLTransWorhp::neben_boundary(double *c_low, double *c_upp){
	for (vector<Variable*>::size_type i=0; i<nebenCount; i++) {
		if (nebenc[i]->hasLowBound()) { c_low[i] = nebenc[i]->getLowBound();}
		if (nebenc[i]->hasUppBound()) { c_upp[i] = nebenc[i]->getUppBound();}
	}
}

bool XMLTransWorhp::neben_structure(DiffStructure &s) {
	if (nebenCount == 0) {
		return false;
	}

	for (int i=0; i<nebenCount; i++) {
		for (int j=0; j<stateCount+controlCount+paramCount; j++) {
			if (diffArray[i+randCount+stateCount][j] == 1) {
				if (j<stateCount) {
					s(i,x_index(0,j));
				} else if(j>=stateCount && j<stateCount+controlCount) {
					s(i,u_index(0,j-stateCount));
				} else if(j>=stateCount+controlCount) {
					s(i,p_indexode(j-stateCount-controlCount));
				}
			}
		}
	}

	return true;
}

bool XMLTransWorhp::neben_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) {
	if (nebenCount == 0) {
		return false;
	}

	for (int i=0; i<nebenCount; i++) {
		for (int j=0; j<stateCount+controlCount+paramCount; j++) {
			if (diffArray[i+randCount+stateCount][j] == 1) {
				if (j<stateCount) {
					s(i,x_index(0,j)) = evaluateDiffOde(nebenc[i]->getTree(),x,u,p,j);
				} else if(j>=stateCount && j<stateCount+controlCount) {
					s(i,u_index(0,j-stateCount)) = evaluateDiffOde(nebenc[i]->getTree(),x,u,p,j);
				}
			}
		}
	}
	return true;
}

bool XMLTransWorhp::neben_diff_p(DiffStructure &s, double t, const double *x, const double *u, const double *p, int index) {
	if (nebenCount == 0) {
		return false;
	}

	for (int i=0; i<nebenCount; i++) {
		for (int j=stateCount+controlCount; j<stateCount+controlCount+paramCount; j++) {
			if (diffArray[i+randCount+stateCount][j] == 1) {
				s(i,p_indexode(j-stateCount-controlCount)) = evaluateDiffOde(nebenc[i]->getTree(),x,u,p,j);
			}
		}
	}

	return true;
}
	
// Startschaetzung
void XMLTransWorhp::x_init(double *x, int i, int dis) {	
	for (vector<Variable*>::size_type j=0; j<stateCount; j++) {
		if (states[j]->hasGuess()) { x[j] = states[j]->getGuess();}
	}
}

void XMLTransWorhp::u_init(double *u, int i, int dis) {
	for (vector<Variable*>::size_type j=0; j<controlCount; j++) {
		if (controls[j]->hasGuess()) { u[j] = controls[j]->getGuess();}
	}
}

void XMLTransWorhp::p_init(double *p) {
	for (vector<Variable*>::size_type j=0; j<paramCount; j++) {
		if (params[j]->hasGuess()) { p[j] = params[j]->getGuess(); }
	}
}
