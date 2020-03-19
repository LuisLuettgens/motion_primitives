/*----------------------------------------------------------------
*
*  
*  Calculation of a new trajectory by WorhpZen given a pertubation in the initial states or a shift in the final states.
*  The new trajectory will be stored in ph.solver
*
*----------------------------------------------------------------*/


#include "CourseCorrection.h"
#include <iostream>
#include "global.h"

using namespace std;

//windows users are not able to call WriteSensitivity(), therefor this function is used in nachkorrektur()
void WriteMySensitivity(std::ostream &os, char var, char pert) {

	int npert = ZenN(&(folder0->worhp_o), &pert);
	int nvar = ZenN(&(folder0->worhp_o), &var);

	std::cout << "D" << var << "/D" << pert << ": " << nvar << " x " << npert << std::endl;
	auto d = new double[nvar]; // Das sollte gross genug sein...

	for (int i = 1; i <= npert; i++) {
		ZenGetD(&(folder0->worhp_o), &(folder0->worhp_w), &(folder0->worhp_p), &(folder0->worhp_c),
			&var, &pert, &i, d);
	
		os.setf(std::ios::scientific);
		os << std::setprecision(15);
		for (int j = 0; j < nvar; j++) {
			os << " " << std::setw(24) << d[j];
		}

		os << std::endl;
	}
	delete[]d;
}


// calculates a new trajectory given a pertubation in the initial states using WorhpZen
void korrektur_startposition(vector<double> pertubation) {

	// constant dimensions for initializing arrays (Visual Studio problem only)
	const int n_dis   = ph0->n_dis;
	const int n_ode   = ph0->n_ode;
	const int n_ctrl  = ph0->n_ctrl;
	const int n_neben = ph0->n_neben;
	const int n_rand  = ph0->n_rand;

	int it = 0; //number of iterations for nachkorrektur

	vector<double> xperturbed(ph0->n_ode); 

	for (int i = 0; i < n_ode; i++) {
		xperturbed[i] = ph0->solver->X[ph0->x_index(0, i)] + pertubation[i];
	}

	char var[] = { 'X','G','F','L','M',0 };
	char pert[] = { 'P','R','Q','B',0 };
	ofstream os((string("ZenD") + 'X' + ".dat").c_str());

	WriteMySensitivity(os, var[0], pert[2]); //DXDQ
	WriteMySensitivity(os, var[0], pert[3]); //DXDB
    //Writes Zen files.

	const int n = (n_ode + n_ctrl) * n_dis + 1; //number of rows in ZenDXDB
	const int m = n_ode  * (n_dis - 1) + n_neben*n_dis; //number of rows in ZenDXDQ

	double DXDB[n][n];
	for (int xxx = 1;xxx < n + 1; xxx++) {
		ZenGetD(&folder0->worhp_o, &folder0->worhp_w, &folder0->worhp_p, &folder0->worhp_c,
			&var[0], &pert[3], &xxx, DXDB[xxx - 1]);
	} //makes DXDB file from ZEN.DX

	double DXDQ[m][n];
	for (int xxx = 1;xxx < m + 1; xxx++) {
		ZenGetD(&folder0->worhp_o, &folder0->worhp_w, &folder0->worhp_p, &folder0->worhp_c,
			&var[0], &pert[2], &xxx, DXDQ[xxx - 1]);
	} //makes DXDQ file from ZEN.DX 


	  //---------------------------Vorkorrektur---------------------------------------------------------------------------------------------	

	double X_0[n];
	int index = 0;
	for (int jjj = 0;jjj < n_dis; jjj++) {
		for (int idx = 0; idx < n_ode + n_ctrl; idx++) {
			if (idx < n_ode) {
				X_0[index] = ph0->solver->X[ph0->x_index(jjj, idx)];
			}
			else {
				X_0[index] = ph0->solver->X[ph0->u_index(jjj, idx - n_ode)];
			}
			index++;
		}
	}
	X_0[index] = ph0->solver->X[ph0->p_index(0)]; //Makes X_0 file from ph0.solver

	double TF; //This is the final time before correction.
	TF = ph0->solver->X[ph0->p_index(0)];

	double E[n_ode];
	for (int idx = 0; idx < ph0->n_ode; idx++) {
		E[idx] = xperturbed[idx] - ph0->solver->X[ph0->x_index(0, idx)];
	} // Calcukates the pertubation of the states

	  //To Update states and controls with sensitivities to states and ODE's (Echzeit2)
	double XX[n_dis][n_ode + n_ctrl];
	for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
			XX[jjj][idx] = X_0[idx + (ph0->n_ode + ph0->n_ctrl)*jjj];
		}
	}// makes matrix of the initial states and controls at each discretization point

	double B[n_ode + n_ctrl][n_ode][n_dis];
	for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
		for (int idx = 0;idx < ph0->n_ode;idx++) {
			for (int jjj = 0;jjj < ph0->n_dis;jjj++) {
				B[k][idx][jjj] = DXDB[idx][k + (ph0->n_ode + ph0->n_ctrl)*jjj];
			}
		}
	}//Makes a 3D matrix containing the sensitivity of the states and controls to the states, where k is the state or control to be updated

	double E_B[n_ode + n_ctrl][n_dis];
	for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
		for (int idx = 0; idx < ph0->n_dis; idx++) {
			E_B[k][idx] = 0.0;
			for (int jjj = 0; jjj < ph0->n_ode;jjj++) {
				E_B[k][idx] = E_B[k][idx] + E[jjj] * B[k][jjj][idx];
			}
		}
	}// Calculates correction vector for each state and control by multiplying E and B

	double X1[n_dis][n_ode + n_ctrl];
	for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
			X1[jjj][idx] = XX[jjj][idx] + E_B[idx][jjj];
		}
	}// First update to the states and controls due to sensitivities to the states with error as the initial pertubation

	double BT[n_ode];
	for (int idx = 0;idx < ph0->n_ode;idx++) {
		BT[idx] = DXDB[idx][n - 1];
	}//Vector of the sensitivity of the final time to the states.

	double TF1;
	double S = 0;
	for (int idx = 0; idx < ph0->n_ode; idx++) {
		S = S + E[idx] * BT[idx];
	} //Calculates correction for the time by multiplying E and BT

	TF1 = TF + S; //first time correction

				  //---------------------------Nachkorrektur---------------------------------------------------------------------------------------------	


	double X2[n_dis][n_ode + n_ctrl];
	for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
			X2[jjj][idx] = X1[jjj][idx];
		}
	}// copies the matrix X1 to X2 for Nachkorrektur

	double TF2 = TF1; //copies TF1 to TF2 for Nachkorrektur


	double dim = ph0->n_dis - 1;

	double x1[n_ode];
	double x2[n_ode];
	double u1[n_ctrl];
	double u2[n_ctrl];
	double dx1[n_ode];
	double dx2[n_ode];
	double EE[n_ode*(n_dis - 1)];
	double BB[n_ode + n_ctrl][n_ode*(n_dis - 1)][n_dis];
	double EE_BB[n_ode + n_ctrl][n_dis];
	double BBT[n_ode*(n_dis - 1)];

	for (int z = 0;z < it;z++) {

		double T2 = TF2 / dim; //equals one time step

		for (int jjj = 0;jjj < ph0->n_dis - 1;jjj++) {

			//Writing states and controls for discretization point j in x1 and u1
			for (int kkk = 0; kkk < ph0->n_ode; kkk++) {
				x1[kkk] = X2[jjj][kkk];
			}
			for (int kkk = 0; kkk < ph0->n_ctrl; kkk++) {
				u1[kkk] = X2[jjj][ph0->n_ode + kkk];
			}
			//Writing states and controls for discretization point j+1 in x2 and u2
			for (int kkk = 0; kkk < ph0->n_ode; kkk++) {
				x2[kkk] = X2[jjj + 1][kkk];
			}
			for (int kkk = 0; kkk < ph0->n_ctrl; kkk++) {
				u2[kkk] = X2[jjj + 1][ph0->n_ode + kkk];
			}

			const double p = 1.0; // To make ode not depend on the end time

			//Calculating the ode's given states and controls
			ph0->ode(dx1, 0.0, x1, u1, &p);
			ph0->ode(dx2, 0.0, x2, u2, &p);


			for (int idx = 0;idx < ph0->n_ode;idx++) {
				EE[idx + (ph0->n_ode)*jjj] = -(X2[jjj + 1][idx] - X2[jjj][idx] - (dx1[idx] + dx2[idx])*(T2 / 2.0));
// 				cout << "EE: " << X2[jjj + 1][idx] << "  " << X2[jjj][idx] << "  " << dx1[idx] << "  " << dx2[idx] << "  " << T2 << endl;
			}
		}// Calculates error between X2 and integrated X2



		for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
			for (int idx = 0;idx < ph0->n_ode*(ph0->n_dis - 1);idx++) {
				for (int jjj = 0;jjj < ph0->n_dis;jjj++) {
					//cout << DXDQ[idx][(k + (ph0->n_ode + ph0->n_ctrl)*jjj)] << endl;
					BB[k][idx][jjj] = DXDQ[idx][(k + (ph0->n_ode + ph0->n_ctrl)*jjj)];
				}
			}
		}//Makes a 3D matrix containing the sensitivity of the states and controls to the ODE's, where k is the state or control to be updated

		for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
			for (int idx = 0; idx < ph0->n_dis; idx++) {
				EE_BB[k][idx] = 0.0;
				for (int jjj = 0; jjj < ph0->n_ode*(ph0->n_dis - 1);jjj++) {
					EE_BB[k][idx] = EE_BB[k][idx] + EE[jjj] * BB[k][jjj][idx];
				}
			}
		}// Calculates correction vector for each state and control by multiplying EE and BB

		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
			for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
				X2[jjj][idx] = X2[jjj][idx] + EE_BB[idx][jjj];
			}
		}// Update to the states and controls due to sensitivities to the ODE's with error as the offset from the trapezoidal rule of integration


		for (int idx = 0;idx < ph0->n_ode*(ph0->n_dis - 1);idx++) {
			BBT[idx] = DXDQ[idx][n - 1];

		}//Vector of the sensitivity of the final time to the ODE's

		double S = 0.0;
		for (int idx = 0; idx < ph0->n_ode*(ph0->n_dis - 1); idx++) {
			S = S + EE[idx] * BBT[idx];
		} //Calculates correction for the time by multiplying EE and BBT
// 		cout << "S: " << S << endl;
		TF2 = TF2 + S; //Time correction

	} //End Nachkorrektur loop



	ofstream f("X2.txt");
	for (int jjj = 0; jjj < n_dis; jjj++) {
		for (int idx = 0; idx < n_ode + n_ctrl; idx++) {
			f << X2[jjj][idx] << endl; //creation of a .txt file with the new trajectory (for use in matlab)
			if (idx < n_ode) {
				ph0->solver->X[ph0->x_index(jjj, idx)] = X2[jjj][idx]; //writing the new states in ph0.solver
			}
			else {
				ph0->solver->X[ph0->u_index(jjj, idx - n_ode)] = X2[jjj][idx]; //writing the new controls in ph0.solver
			}
		}
	} 
	f.close();

}

//---------------------------------------------------------------------------------------------------------------------------------------------------------

// calculates a new trajectory given a shift in the final states using WorhpZen
void korrektur_endposition(vector<double> shift) {

	// constant dimensions for initializing arrays (Visual Studio problem only)
	const int n_dis = 23;
	const int n_ode = 7;
	const int n_ctrl = 2;

	int it = 100; //number of iterations for nachkorrektur

	vector<double> xshifted(ph0->n_ode);

	for (int i = 0; i < n_ode; i++) {
		xshifted[i] = ph0->solver->X[ph0->x_index(ph0->n_dis-1, i)] + shift[i]; //adds shift to the endposition
	}


	char var[] = { 'X','G','F','L','M',0 };
	char pert[] = { 'P','R','Q','B',0 };
	ofstream os((string("ZenD") + 'X' + ".dat").c_str());

	WriteMySensitivity(os, var[0], pert[2]); //DXDQ
	WriteMySensitivity(os, var[0], pert[3]); //DXDB
	//Writes Zen files.

	const int n = (n_ode + n_ctrl) * n_dis + 1; //number of rows in ZenDXDB
	const int m = n_ode  * (n_dis - 1) + n_dis; //number of rows in ZenDXDQ

	double DXDB[n][n];
	for (int xxx = 1;xxx < n + 1; xxx++) {
		ZenGetD(&folder0->worhp_o, &folder0->worhp_w, &folder0->worhp_p, &folder0->worhp_c,
			&var[0], &pert[3], &xxx, DXDB[xxx - 1]);
	} //makes DXDB file from ZEN.DX

	double DXDQ[m][n];
	for (int xxx = 1;xxx < m + 1; xxx++) {
		ZenGetD(&folder0->worhp_o, &folder0->worhp_w, &folder0->worhp_p, &folder0->worhp_c,
			&var[0], &pert[2], &xxx, DXDQ[xxx - 1]);
	} //makes DXDQ file from ZEN.DX 


	  //---------------------------Vorkorrektur---------------------------------------------------------------------------------------------	

	double X_0[n];
	int index = 0;
	for (int jjj = 0;jjj < n_dis; jjj++) {
		for (int idx = 0; idx < n_ode + n_ctrl; idx++) {
			if (idx < n_ode) {
				X_0[index] = ph0->solver->X[ph0->x_index(jjj, idx)];
			}
			else {
				X_0[index] = ph0->solver->X[ph0->u_index(jjj, idx - n_ode)];
			}
			index++;
		}
	}
	X_0[index] = ph0->solver->X[ph0->p_index(0)]; //Makes X_0 file from ph0.solver

	double TF; //This is the final time before correction.
	TF = ph0->solver->X[ph0->p_index(0)];

	double E[n_ode];
	for (int idx = 0; idx < ph0->n_ode; idx++) {
		E[idx] = xshifted[idx] - ph0->solver->X[ph0->x_index(ph0->n_dis-1, idx)];
	} // Calculates the pertubation of the final states

	  //To Update states and controls with sensitivities to states and ODE's (Echzeit2)
	double XX[n_dis][n_ode + n_ctrl];
	for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
			XX[jjj][idx] = X_0[idx + (ph0->n_ode + ph0->n_ctrl)*jjj];
		}
	}// makes matrix of the initial states and controls at each discretization point

	double B[n_ode + n_ctrl][n_ode][n_dis];
	for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
		for (int idx = n-(ph0->n_ode+ph0->n_ctrl+1);idx < n-(ph0->n_ctrl+1);idx++) {
			for (int jjj = 0;jjj < ph0->n_dis;jjj++) {
				B[k][idx-(n - (ph0->n_ode + ph0->n_ctrl + 1))][jjj] = DXDB[idx][k + (ph0->n_ode + ph0->n_ctrl)*jjj];
			}
		}
	}//Makes a 3D matrix containing the sensitivity of the states and controls to the states, where k is the state or control to be updated

	double E_B[n_ode + n_ctrl][n_dis];
	for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
		for (int idx = 0; idx < ph0->n_dis; idx++) {
			E_B[k][idx] = 0.0;
			for (int jjj = 0; jjj < ph0->n_ode;jjj++) {
				E_B[k][idx] = E_B[k][idx] + E[jjj] * B[k][jjj][idx];
			}
		}
	}// Calculates correction vector for each state and control by multiplying E and B

	double X1[n_dis][n_ode + n_ctrl];
	for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
			X1[jjj][idx] = XX[jjj][idx] + E_B[idx][jjj];
		}
	}// First update to the states and controls due to sensitivities to the states with error as the initial pertubation


	double BT[n_ode];
	for (int idx = 0;idx < ph0->n_ode;idx++) {
		BT[idx] = DXDB[idx][n - 1];
	}//Vector of the sensitivity of the final time to the states.

	double TF1;
	double S = 0;
	for (int idx = 0; idx < ph0->n_ode; idx++) {
		S = S + E[idx] * BT[idx];
	} //Calculates correction for the time by multiplying E and BT

	TF1 = TF + S; //first time correction

				  //---------------------------Nachkorrektur---------------------------------------------------------------------------------------------	


	double X2[n_dis][n_ode + n_ctrl];
	for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
			X2[jjj][idx] = X1[jjj][idx];
		}
	}// copies the matrix X1 to X2 for Nachkorrektur

	double TF2 = TF1; //copies TF1 to TF2 for Nachkorrektur


	double dim = ph0->n_dis - 1;

	double x1[n_ode];
	double x2[n_ode];
	double u1[n_ctrl];
	double u2[n_ctrl];
	double dx1[n_ode];
	double dx2[n_ode];
	double EE[n_ode*(n_dis - 1)];
	double BB[n_ode + n_ctrl][n_ode*(n_dis - 1)][n_dis];
	double EE_BB[n_ode + n_ctrl][n_dis];
	double BBT[n_ode*(n_dis - 1)];

	for (int z = 0;z < it;z++) {

		double T2 = TF2 / dim; //equals one time step

		for (int jjj = 0;jjj < ph0->n_dis - 1;jjj++) {

			//Writing states and controls for discretization point j in x1 and u1
			for (int kkk = 0; kkk < ph0->n_ode; kkk++) {
				x1[kkk] = X2[jjj][kkk];
			}
			for (int kkk = 0; kkk < ph0->n_ctrl; kkk++) {
				u1[kkk] = X2[jjj][ph0->n_ode + kkk];
			}
			//Writing states and controls for discretization point j+1 in x2 and u2
			for (int kkk = 0; kkk < ph0->n_ode; kkk++) {
				x2[kkk] = X2[jjj + 1][kkk];
			}
			for (int kkk = 0; kkk < ph0->n_ctrl; kkk++) {
				u2[kkk] = X2[jjj + 1][ph0->n_ode + kkk];
			}

			const double p = 1.0; // To make ode not depend on the end time

			//Calculating the ode's given states and controls
			ph0->ode(dx1, 0.0, x1, u1, &p);
			ph0->ode(dx2, 0.0, x2, u2, &p);

			for (int idx = 0;idx < ph0->n_ode;idx++) {
				EE[idx + (ph0->n_ode)*jjj] = -(X2[jjj + 1][idx] - X2[jjj][idx] - (dx1[idx] + dx2[idx])*(T2 / 2.0));
			}
		}// Calculates error between X2 and integrated X2



		for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
			for (int idx = 0;idx < ph0->n_ode*(ph0->n_dis - 1);idx++) {
				for (int jjj = 0;jjj < ph0->n_dis;jjj++) {
					//cout << DXDQ[idx][(k + (ph0->n_ode + ph0->n_ctrl)*jjj)] << endl;
					BB[k][idx][jjj] = DXDQ[idx][(k + (ph0->n_ode + ph0->n_ctrl)*jjj)];
				}
			}
		}//Makes a 3D matrix containing the sensitivity of the states and controls to the ODE's, where k is the state or control to be updated

		for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
			for (int idx = 0; idx < ph0->n_dis; idx++) {
				EE_BB[k][idx] = 0.0;
				for (int jjj = 0; jjj < ph0->n_ode*(ph0->n_dis - 1);jjj++) {
					EE_BB[k][idx] = EE_BB[k][idx] + EE[jjj] * BB[k][jjj][idx];
				}
			}
		}// Calculates correction vector for each state and control by multiplying EE and BB

		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
			for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
				X2[jjj][idx] = X2[jjj][idx] + EE_BB[idx][jjj];
			}
		}// Update to the states and controls due to sensitivities to the ODE's with error as the offset from the trapezoidal rule of integration


		for (int idx = 0;idx < ph0->n_ode*(ph0->n_dis - 1);idx++) {
			BBT[idx] = DXDQ[idx][n - 1];

		}//Vector of the sensitivity of the final time to the ODE's

		double S = 0.0;
		for (int idx = 0; idx < ph0->n_ode*(ph0->n_dis - 1); idx++) {
			S = S + EE[idx] * BBT[idx];
		} //Calculates correction for the time by multiplying EE and BBT

		TF2 = TF2 + S; //Time correction

	} //End Nachkorrektur loop



	ofstream f("X2.txt");
	for (int jjj = 0; jjj < n_dis; jjj++) {
		for (int idx = 0; idx < n_ode + n_ctrl; idx++) {
			f << X2[jjj][idx] << endl; //creation of a .txt file with the new trajectory (for use in matlab)
			if (idx < n_ode) {
				ph0->solver->X[ph0->x_index(jjj, idx)] = X2[jjj][idx]; //writing the new states in ph0.solver
			}
			else {
				ph0->solver->X[ph0->u_index(jjj, idx - n_ode)] = X2[jjj][idx]; //writing the new controls in ph0.solver
			}
		}
	}
	f.close();

}



//--------------------------------------------------------------------------------------------------------------------------

// both korrektur_startposition and korrektur_endposition in one function
void korrektur(vector<double> pertubation, vector<double> shift) {


	// constant dimensions for initializing arrays (Visual Studio problem only)
	const int n_dis = 23;
	const int n_ode = 7;
	const int n_ctrl = 2;

	int it = 100; //number of iterations for nachkorrektur
	
	vector<double> xperturbed(ph0->n_ode);

	for (int i = 0; i < n_ode; i++) {
		xperturbed[i] = ph0->solver->X[ph0->x_index(0, i)] + pertubation[i];
	}

	vector<double> xshifted(ph0->n_ode);

	for (int i = 0; i < n_ode; i++) {
		xshifted[i] = ph0->solver->X[ph0->x_index(ph0->n_dis - 1, i)] + shift[i]; //adds shift to the endposition
	}


	char var[] = { 'X','G','F','L','M',0 };
	char pert[] = { 'P','R','Q','B',0 };
	ofstream os((string("ZenD") + 'X' + ".dat").c_str());

	WriteMySensitivity(os, var[0], pert[2]); //DXDQ
	WriteMySensitivity(os, var[0], pert[3]); //DXDB
	//Writes Zen files.

	const int n = (n_ode + n_ctrl) * n_dis + 1; //number of rows in ZenDXDB
	const int m = n_ode  * (n_dis - 1) + n_dis; //number of rows in ZenDXDQ

	double DXDB[n][n];
	for (int xxx = 1;xxx < n + 1; xxx++) {
		ZenGetD(&folder0->worhp_o, &folder0->worhp_w, &folder0->worhp_p, &folder0->worhp_c,
			&var[0], &pert[3], &xxx, DXDB[xxx - 1]);
	} //makes DXDB file from ZEN.DX

	double DXDQ[m][n];
	for (int xxx = 1;xxx < m + 1; xxx++) {
		ZenGetD(&folder0->worhp_o, &folder0->worhp_w, &folder0->worhp_p, &folder0->worhp_c,
			&var[0], &pert[2], &xxx, DXDQ[xxx - 1]);
	} //makes DXDQ file from ZEN.DX 


	//---------------------------Vorkorrektur---------------------------------------------------------------------------------------------	

	double X_0[n];
	int index = 0;
	for (int jjj = 0;jjj < n_dis; jjj++) {
		for (int idx = 0; idx < n_ode + n_ctrl; idx++) {
			if (idx < n_ode) {
				X_0[index] = ph0->solver->X[ph0->x_index(jjj, idx)];
			}
			else {
				X_0[index] = ph0->solver->X[ph0->u_index(jjj, idx - n_ode)];
			}
			index++;
		}
	}
	X_0[index] = ph0->solver->X[ph0->p_index(0)]; //Makes X_0 file from ph0.solver

	double TF; //This is the final time before correction.
	TF = ph0->solver->X[ph0->p_index(0)];

	double EP[n_ode];
	for (int idx = 0; idx < ph0->n_ode; idx++) {
		EP[idx] = xperturbed[idx] - ph0->solver->X[ph0->x_index(0, idx)];
	} // Calcukates the pertubation of the states

	double ES[n_ode];
	for (int idx = 0; idx < ph0->n_ode; idx++) {
		ES[idx] = xshifted[idx] - ph0->solver->X[ph0->x_index(ph0->n_dis - 1, idx)];
	} // Calculates the pertubation of the final states

	  //To Update states and controls with sensitivities to states and ODE's (Echzeit2)
	double XX[n_dis][n_ode + n_ctrl];
	for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
			XX[jjj][idx] = X_0[idx + (ph0->n_ode + ph0->n_ctrl)*jjj];
		}
	}// makes matrix of the initial states and controls at each discretization point

	double BP[n_ode + n_ctrl][n_ode][n_dis];
	for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
		for (int idx = 0;idx < ph0->n_ode;idx++) {
			for (int jjj = 0;jjj < ph0->n_dis;jjj++) {
				BP[k][idx][jjj] = DXDB[idx][k + (ph0->n_ode + ph0->n_ctrl)*jjj];
			}
		}
	}//Makes a 3D matrix containing the sensitivity of the states and controls to the states, where k is the state or control to be updated

	double BS[n_ode + n_ctrl][n_ode][n_dis];
	for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
		for (int idx = n - (ph0->n_ode + ph0->n_ctrl + 1);idx < n - (ph0->n_ctrl + 1);idx++) {
			for (int jjj = 0;jjj < ph0->n_dis;jjj++) {
				BS[k][idx - (n - (ph0->n_ode + ph0->n_ctrl + 1))][jjj] = DXDB[idx][k + (ph0->n_ode + ph0->n_ctrl)*jjj];
			}
		}
	}//Makes a 3D matrix containing the sensitivity of the states and controls to the states, where k is the state or control to be updated

	double E_B[n_ode + n_ctrl][n_dis];
	for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
		for (int idx = 0; idx < ph0->n_dis; idx++) {
			E_B[k][idx] = 0.0;
			for (int jjj = 0; jjj < ph0->n_ode;jjj++) {
				E_B[k][idx] = E_B[k][idx] + ES[jjj] * BS[k][jjj][idx] + EP[jjj] * BP[k][jjj][idx];
			}
		}
	}// Calculates correction vector for each state and control by multiplying E and B

	double X1[n_dis][n_ode + n_ctrl];
	for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
			X1[jjj][idx] = XX[jjj][idx] + E_B[idx][jjj];
		}
	}// First update to the states and controls due to sensitivities to the states with error as the initial pertubation


	double BT[n_ode];
	for (int idx = 0;idx < ph0->n_ode;idx++) {
		BT[idx] = DXDB[idx][n - 1];
	}//Vector of the sensitivity of the final time to the states.

	double TF1;
	double S = 0;
	for (int idx = 0; idx < ph0->n_ode; idx++) {
		S = S + EP[idx] * BT[idx] + ES[idx] * BT[idx];
	} //Calculates correction for the time by multiplying E and BT

	TF1 = TF + S; //first time correction

	//---------------------------Nachkorrektur---------------------------------------------------------------------------------------------	


	double X2[n_dis][n_ode + n_ctrl];
	for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
			X2[jjj][idx] = X1[jjj][idx];
		}
	}// copies the matrix X1 to X2 for Nachkorrektur

	double TF2 = TF1; //copies TF1 to TF2 for Nachkorrektur


	double dim = ph0->n_dis - 1;

	double x1[n_ode];
	double x2[n_ode];
	double u1[n_ctrl];
	double u2[n_ctrl];
	double dx1[n_ode];
	double dx2[n_ode];
	double EE[n_ode*(n_dis - 1)];
	double BB[n_ode + n_ctrl][n_ode*(n_dis - 1)][n_dis];
	double EE_BB[n_ode + n_ctrl][n_dis];
	double BBT[n_ode*(n_dis - 1)];

	for (int z = 0;z < it;z++) {

		double T2 = TF2 / dim; //equals one time step

		for (int jjj = 0;jjj < ph0->n_dis - 1;jjj++) {

			//Writing states and controls for discretization point j in x1 and u1
			for (int kkk = 0; kkk < ph0->n_ode; kkk++) {
				x1[kkk] = X2[jjj][kkk];
			}
			for (int kkk = 0; kkk < ph0->n_ctrl; kkk++) {
				u1[kkk] = X2[jjj][ph0->n_ode + kkk];
			}
			//Writing states and controls for discretization point j+1 in x2 and u2
			for (int kkk = 0; kkk < ph0->n_ode; kkk++) {
				x2[kkk] = X2[jjj + 1][kkk];
			}
			for (int kkk = 0; kkk < ph0->n_ctrl; kkk++) {
				u2[kkk] = X2[jjj + 1][ph0->n_ode + kkk];
			}

			const double p = 1.0; // To make ode not depend on the end time

			//Calculating the ode's given states and controls (time independend, therefor it's t = 0.0)
			ph0->ode(dx1, 0.0, x1, u1, &p);
			ph0->ode(dx2, 0.0, x2, u2, &p);

			for (int idx = 0;idx < ph0->n_ode;idx++) {
				EE[idx + (ph0->n_ode)*jjj] = -(X2[jjj + 1][idx] - X2[jjj][idx] - (dx1[idx] + dx2[idx])*(T2 / 2.0));
			}
		}// Calculates error between X2 and integrated X2


		for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
			for (int idx = 0;idx < ph0->n_ode*(ph0->n_dis - 1);idx++) {
				for (int jjj = 0;jjj < ph0->n_dis;jjj++) {
					//cout << DXDQ[idx][(k + (ph0->n_ode + ph0->n_ctrl)*jjj)] << endl;
					BB[k][idx][jjj] = DXDQ[idx][(k + (ph0->n_ode + ph0->n_ctrl)*jjj)];
				}
			}
		}//Makes a 3D matrix containing the sensitivity of the states and controls to the ODE's, where k is the state or control to be updated

		for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
			for (int idx = 0; idx < ph0->n_dis; idx++) {
				EE_BB[k][idx] = 0.0;
				for (int jjj = 0; jjj < ph0->n_ode*(ph0->n_dis - 1);jjj++) {
					EE_BB[k][idx] = EE_BB[k][idx] + EE[jjj] * BB[k][jjj][idx];
				}
			}
		}// Calculates correction vector for each state and control by multiplying EE and BB

		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
			for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
				X2[jjj][idx] = X2[jjj][idx] + EE_BB[idx][jjj];
			}
		}// Update to the states and controls due to sensitivities to the ODE's with error as the offset from the trapezoidal rule of integration


		for (int idx = 0;idx < ph0->n_ode*(ph0->n_dis - 1);idx++) {
			BBT[idx] = DXDQ[idx][n - 1];

		}//Vector of the sensitivity of the final time to the ODE's

		double S = 0.0;
		for (int idx = 0; idx < ph0->n_ode*(ph0->n_dis - 1); idx++) {
			S = S + EE[idx] * BBT[idx];
		} //Calculates correction for the time by multiplying EE and BBT

		TF2 = TF2 + S; //Time correction

	} //End Nachkorrektur loop



	ofstream f("X2.txt");
	for (int jjj = 0; jjj < n_dis; jjj++) {
		for (int idx = 0; idx < n_ode + n_ctrl; idx++) {
			f << X2[jjj][idx] << endl; //creation of a .txt file with the new trajectory (for use in matlab)
			if (idx < n_ode) {
				ph0->solver->X[ph0->x_index(jjj, idx)] = X2[jjj][idx]; //writing the new states in ph0.solver
			}
			else {
				ph0->solver->X[ph0->u_index(jjj, idx - n_ode)] = X2[jjj][idx]; //writing the new controls in ph0.solver
			}
		}
	}
	f.close();


}







// 
// 
// // calculates a new trajectory given a pertubation in the initial states using WorhpZen
// void korrektur_startposition(vector<double> pertubation, void(*system)(double*, double*)) {
// 
// 	// constant dimensions for initializing arrays (Visual Studio problem only)
// 	const int n_dis  = ph0->n_dis;
// 	const int n_ode  = ph0->n_ode;
// 	const int n_ctrl = ph0->n_ctrl;
// 
// 	int it = 100; //number of iterations for nachkorrektur
// 
// 	vector<double> xperturbed(ph0->n_ode); 
// 
// 	for (int i = 0; i < n_ode; i++) {
// 		xperturbed[i] = ph0->solver->X[ph0->x_index(0, i)] + pertubation[i];
// 	}
// 
// 	char var[] = { 'X','G','F','L','M',0 };
// 	char pert[] = { 'P','R','Q','B',0 };
// 	int i = 0;
// 	while (var[i]) {
// 		ofstream os((string("ZenD") + var[i] + ".dat").c_str());
// 		int j = 0;
// 		while (pert[j]) {
// 			WriteMySensitivity(os, var[i], pert[j]);
// 			j++;
// 		}
// 		i++;
// 	} //Writes Zen files.
// 
// 
// 	const int n = (n_ode + n_ctrl) * n_dis + 1; //number of rows in ZenDXDB
// 	const int m = n_ode  * (n_dis - 1) + n_dis; //number of rows in ZenDXDQ
// 
// 	double DXDB[n][n];
// 	for (int xxx = 1;xxx < n + 1; xxx++) {
// 		ZenGetD(&(folder0->worhp_o), &(folder0->worhp_w), &(folder0->worhp_p), &(folder0->worhp_c),
// 			&var[0], &pert[3], &xxx, DXDB[xxx - 1]);
// 	} //makes DXDB file from ZEN.DX
// 	double DXDQ[m][n];
// 	for (int xxx = 1;xxx < m + 1; xxx++) {
// 		ZenGetD(&(folder0->worhp_o), &(folder0->worhp_w), &(folder0->worhp_p), &(folder0->worhp_c),
// 			&var[0], &pert[2], &xxx, DXDQ[xxx - 1]);
// 	} //makes DXDQ file from ZEN.DX 
// 
// 	  //---------------------------Vorkorrektur---------------------------------------------------------------------------------------------	
// 
// 	double X_0[n];
// 	int index = 0;
// 	for (int jjj = 0;jjj < n_dis; jjj++) {
// 		for (int idx = 0; idx < n_ode + n_ctrl; idx++) {
// 			if (idx < n_ode) {
// 				X_0[index] = ph0->solver->X[ph0->x_index(jjj, idx)];
// 			}
// 			else {
// 				X_0[index] = ph0->solver->X[ph0->u_index(jjj, idx - n_ode)];
// 			}
// 			index++;
// 		}
// 	}
// 	X_0[index] = ph0->solver->X[ph0->p_index(0)]; //Makes X_0 file from ph.solver
// 
// 	double TF; //This is the final time before correction.
// 	TF = ph0->solver->X[ph0->p_index(0)];
// 
// 	double X_start[n_ode];
// 	for (int idx = 0;idx < ph0->n_ode; idx++) {
// 		X_start[idx] = X_0[idx];
// 	} // makes vector of initial states
// 
// 	double E[n_ode];
// 	for (int idx = 0; idx < ph0->n_ode; idx++) {
// 		E[idx] = xperturbed[idx] - ph0->solver->X[ph0->x_index(0, idx)];
// 	} // Calcukates the pertubation of the states
// 
// 	  //To Update states and controls with sensitivities to states and ODE's (Echzeit2)
// 	double XX[n_dis][n_ode + n_ctrl];
// 	for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
// 		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
// 			XX[jjj][idx] = X_0[idx + (ph0->n_ode + ph0->n_ctrl)*jjj];
// 		}
// 	}// makes matrix of the initial states and controls at each discretization point
// 
// 	double B[n_ode + n_ctrl][n_ode][n_dis];
// 	for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
// 		for (int idx = 0;idx < ph0->n_ode;idx++) {
// 			for (int jjj = 0;jjj < ph0->n_dis;jjj++) {
// 				B[k][idx][jjj] = DXDB[idx][k + (ph0->n_ode + ph0->n_ctrl)*jjj];
// 			}
// 		}
// 	}//Makes a 3D matrix containing the sensitivity of the states and controls to the states, where k is the state or control to be updated
// 
// 	double E_B[n_ode + n_ctrl][n_dis];
// 	for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
// 		for (int idx = 0; idx < ph0->n_dis; idx++) {
// 			E_B[k][idx] = 0.0;
// 			for (int jjj = 0; jjj < ph0->n_ode;jjj++) {
// 				E_B[k][idx] = E_B[k][idx] + E[jjj] * B[k][jjj][idx];
// 			}
// 		}
// 	}// Calculates correction vector for each state and control by multiplying E and B
// 
// 	double X1[n_dis][n_ode + n_ctrl];
// 	for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
// 		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
// 			X1[jjj][idx] = XX[jjj][idx] + E_B[idx][jjj];
// 		}
// 	}// First update to the states and controls due to sensitivities to the states with error as the initial pertubation
// 
// 	double BT[n_ode];
// 	for (int idx = 0;idx < ph0->n_ode;idx++) {
// 		BT[idx] = DXDB[idx][n - 1];
// 	}//Vector of the sensitivity of the final time to the states.
// 
// 	double TF1;
// 	double S = 0;
// 	for (int idx = 0; idx < ph0->n_ode; idx++) {
// 		S = S + E[idx] * BT[idx];
// 	} //Calculates correction for the time by multiplying E and BT
// 
// 	TF1 = TF + S; //first time correction
// 
// 				  //---------------------------Nachkorrektur---------------------------------------------------------------------------------------------	
// 
// 
// 	double X2[n_dis][n_ode + n_ctrl];
// 	for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
// 		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
// 			X2[jjj][idx] = X1[jjj][idx];
// 		}
// 	}// copies the matrix X1 to X2 for Nachkorrektur
// 
// 	double TF2 = TF1; //copies TF1 to TF2 for Nachkorrektur
// 
// 
// 	double dim = ph0->n_dis - 1;
// 
// 	double x1[n_ode];
// 	double x2[n_ode];
// 	double u1[n_ctrl];
// 	double u2[n_ctrl];
// 	double dx1[n_ode];
// 	double dx2[n_ode];
// 	vector<double> EE(n_ode*(n_dis - 1));
// 	vector<double> BBT(n_ode*(n_dis - 1));
// 	vector<vector<double>> EE_BB(n_ode + n_ctrl, vector<double>(n_dis));
// 	vector<vector<vector<double>>> BB(n_ode + n_ctrl, vector<vector<double>>(n_ode*(n_dis - 1), vector<double>(n_dis)));
// 	
// 	
// 
// 	for (int z = 0;z < it;z++) {
// 
// 		double T2 = TF2 / dim; //equals one time step
// 		
// 		
// // 		for (int jjj = 0;jjj < ph0->n_dis - 1;jjj++) {
// // 
// // 			for (int kkk = 0; kkk < ph0->n_ode + ph0->n_ctrl; kkk++) {
// // 				x1[kkk] = X2[jjj][kkk];
// // 			}
// // 
// // 			for (int kkk = 0; kkk < ph0->n_ode + ph0->n_ctrl; kkk++) {
// // 				x2[kkk] = X2[jjj + 1][kkk];
// // 			}
// // 
// // 			void(*Funk)(double*, double*) = system;
// // 			Funk(x1, dx1);
// // 
// // 			Funk(x2, dx2);
// // 
// // 			for (int idx = 0;idx < ph0->n_ode;idx++) {
// // 				EE[idx + (ph0->n_ode)*jjj] = -(X2[jjj + 1][idx] - X2[jjj][idx] - (dx1[idx] + dx2[idx])*(T2 / 2.0));
// // 			}
// // 		}// Calculates error between X2 and integrated X2
// 
// 		for (int jjj = 0;jjj < ph0->n_dis - 1;jjj++) {
// 
// 			//Writing states and controls for discretization point j in x1 and u1
// 			for (int kkk = 0; kkk < ph0->n_ode; kkk++) {
// 				x1[kkk] = X2[jjj][kkk];
// 			}
// 			for (int kkk = 0; kkk < ph0->n_ctrl; kkk++) {
// 				u1[kkk] = X2[jjj][ph0->n_ode + kkk];
// 			}
// 			//Writing states and controls for discretization point j+1 in x2 and u2
// 			for (int kkk = 0; kkk < ph0->n_ode; kkk++) {
// 				x2[kkk] = X2[jjj + 1][kkk];
// 			}
// 			for (int kkk = 0; kkk < ph0->n_ctrl; kkk++) {
// 				u2[kkk] = X2[jjj + 1][ph0->n_ode + kkk];
// 			}
// 
// 			const double p = 1.0; // To make ode not depend on the end time
// 
// 			//Calculating the ode's given states and controls
// 			ph0->ode(dx1, 0.0, x1, u1, &p);
// 			ph0->ode(dx2, 0.0, x2, u2, &p);
// 
// 			for (int idx = 0;idx < ph0->n_ode;idx++) {
// 				EE[idx + (ph0->n_ode)*jjj] = -(X2[jjj + 1][idx] - X2[jjj][idx] - (dx1[idx] + dx2[idx])*(T2 / 2.0));
// 			}
// 		}// Calculates error between X2 and integrated X2
// 
// 
// 
// 		for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
// 			for (int idx = 0;idx < ph0->n_ode*(ph0->n_dis - 1);idx++) {
// 				for (int jjj = 0;jjj < ph0->n_dis;jjj++) {
// 					//cout << DXDQ[idx][(k + (ph0->n_ode + ph0->n_ctrl)*jjj)] << endl;
// 					BB[k][idx][jjj] = DXDQ[idx][(k + (ph0->n_ode + ph0->n_ctrl)*jjj)];
// 				}
// 			}
// 		}//Makes a 3D matrix containing the sensitivity of the states and controls to the ODE's, where k is the state or control to be updated
// 
// 		for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
// 			for (int idx = 0; idx < ph0->n_dis; idx++) {
// 				EE_BB[k][idx] = 0.0;
// 				for (int jjj = 0; jjj < ph0->n_ode*(ph0->n_dis - 1);jjj++) {
// 					EE_BB[k][idx] = EE_BB[k][idx] + EE[jjj] * BB[k][jjj][idx];
// 				}
// 			}
// 		}// Calculates correction vector for each state and control by multiplying EE and BB
// 
// 		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
// 			for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
// 				X2[jjj][idx] = X2[jjj][idx] + EE_BB[idx][jjj];
// 			}
// 		}// Update to the states and controls due to sensitivities to the ODE's with error as the offset from the trapezoidal rule of integration
// 
// 
// 		for (int idx = 0;idx < ph0->n_ode*(ph0->n_dis - 1);idx++) {
// 			BBT[idx] = DXDQ[idx][n - 1];
// 
// 		}//Vector of the sensitivity of the final time to the ODE's
// 
// 		double S = 0.0;
// 		for (int idx = 0; idx < ph0->n_ode*(ph0->n_dis - 1); idx++) {
// 			S = S + EE[idx] * BBT[idx];
// 		} //Calculates correction for the time by multiplying EE and BBT
// 
// 		TF2 = TF2 + S; //Time correction
// 
// // 		cout << "hello world!" << endl;
// 	} //End Nachkorrektur loop
// 
// 
// 
// 	ofstream f("X2.txt");
// 	for (int jjj = 0; jjj < n_dis; jjj++) {
// 		for (int idx = 0; idx < n_ode + n_ctrl; idx++) {
// 			f << X2[jjj][idx] << endl; //creation of a .txt file with the new trajectory (for use in matlab)
// 			if (idx < n_ode) {
// 				ph0->solver->X[ph0->x_index(jjj, idx)] = X2[jjj][idx]; //writing the new states in ph.solver
// 			}
// 			else {
// 				ph0->solver->X[ph0->u_index(jjj, idx - n_ode)] = X2[jjj][idx]; //writing the new controls in ph.solver
// 			}
// 		}
// 	} 
// 	f.close();
// 
// }
// 
// //---------------------------------------------------------------------------------------------------------------------------------------------------------
// 
// // calculates a new trajectory given a shift in the final states using WorhpZen
// void korrektur_endposition(vector<double> shift, void(*system)(double*, double*)) {
// 
// 	// constant dimensions for initializing arrays (Visual Studio problem only)
// 	const int n_dis  = ph0->n_dis;
// 	const int n_ode  = ph0->n_ode;
// 	const int n_ctrl = ph0->n_ctrl;
// 
// 	int it = 100; //number of iterations for nachkorrektur
// 
// 	vector<double> xshifted(ph0->n_ode);
// 
// 	for (int i = 0; i < n_ode; i++) {
// 		xshifted[i] = ph0->solver->X[ph0->x_index(ph0->n_dis-1, i)] + shift[i]; //adds shift to the endposition
// 	}
// 
// 	char var[] = { 'X','G','F','L','M',0 };
// 	char pert[] = { 'P','R','Q','B',0 };
// 	int i = 0;
// 	while (var[i]) {
// 		ofstream os((string("ZenD") + var[i] + ".dat").c_str());
// 		int j = 0;
// 		while (pert[j]) {
// 			WriteMySensitivity(os, var[i], pert[j]);
// 			j++;
// 		}
// 		i++;
// 	} //Writes Zen files.
// 
// 	const int n = (n_ode + n_ctrl) * n_dis + 1; //number of rows in ZenDXDB
// 	const int m = n_ode  * (n_dis - 1) + n_dis; //number of rows in ZenDXDQ
// 
// 	double DXDB[n][n];
// 	for (int xxx = 1;xxx < n + 1; xxx++) {
// 		ZenGetD(&(folder0->worhp_o), &(folder0->worhp_w), &(folder0->worhp_p), &(folder0->worhp_c),
// 			&var[0], &pert[3], &xxx, DXDB[xxx - 1]);
// 	} //makes DXDB file from ZEN.DX
// 
// 	double DXDQ[m][n];
// 	for (int xxx = 1;xxx < m + 1; xxx++) {
// 		ZenGetD(&(folder0->worhp_o), &(folder0->worhp_w), &(folder0->worhp_p), &(folder0->worhp_c),
// 			&var[0], &pert[2], &xxx, DXDQ[xxx - 1]);
// 	} //makes DXDQ file from ZEN.DX 
// 
// 
// 	  //---------------------------Vorkorrektur---------------------------------------------------------------------------------------------	
// 
// 	double X_0[n];
// 	int index = 0;
// 	for (int jjj = 0;jjj < n_dis; jjj++) {
// 		for (int idx = 0; idx < n_ode + n_ctrl; idx++) {
// 			if (idx < n_ode) {
// 				X_0[index] = ph0->solver->X[ph0->x_index(jjj, idx)];
// 			}
// 			else {
// 				X_0[index] = ph0->solver->X[ph0->u_index(jjj, idx - n_ode)];
// 			}
// 			index++;
// 		}
// 	}
// 	X_0[index] = ph0->solver->X[ph0->p_index(0)]; //Makes X_0 file from ph.solver
// 
// 	double TF; //This is the final time before correction.
// 	TF = ph0->solver->X[ph0->p_index(0)];
// 
// 	//double X_start[n_ode];
// 	//for (int idx = 0;idx < ph0->n_ode; idx++) {
// 	//	X_start[idx] = X_0[idx];
// 	//} // makes vector of initial states
// 
// 	double E[n_ode];
// 	for (int idx = 0; idx < ph0->n_ode; idx++) {
// 		E[idx] = xshifted[idx] - ph0->solver->X[ph0->x_index(ph0->n_dis-1, idx)];
// 	} // Calculates the pertubation of the final states
// 
// 	  //To Update states and controls with sensitivities to states and ODE's (Echzeit2)
// 	double XX[n_dis][n_ode + n_ctrl];
// 	for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
// 		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
// 			XX[jjj][idx] = X_0[idx + (ph0->n_ode + ph0->n_ctrl)*jjj];
// 		}
// 	}// makes matrix of the initial states and controls at each discretization point
// 
// 	double B[n_ode + n_ctrl][n_ode][n_dis];
// 	for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
// 		for (int idx = n-(ph0->n_ode+ph0->n_ctrl+1);idx < n-(ph0->n_ctrl+1);idx++) {
// 			for (int jjj = 0;jjj < ph0->n_dis;jjj++) {
// 				B[k][idx-(n - (ph0->n_ode + ph0->n_ctrl + 1))][jjj] = DXDB[idx][k + (ph0->n_ode + ph0->n_ctrl)*jjj];
// 			}
// 		}
// 	}//Makes a 3D matrix containing the sensitivity of the states and controls to the states, where k is the state or control to be updated
// 
// 	double E_B[n_ode + n_ctrl][n_dis];
// 	for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
// 		for (int idx = 0; idx < ph0->n_dis; idx++) {
// 			E_B[k][idx] = 0.0;
// 			for (int jjj = 0; jjj < ph0->n_ode;jjj++) {
// 				E_B[k][idx] = E_B[k][idx] + E[jjj] * B[k][jjj][idx];
// 			}
// 		}
// 	}// Calculates correction vector for each state and control by multiplying E and B
// 
// 	double X1[n_dis][n_ode + n_ctrl];
// 	for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
// 		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
// 			X1[jjj][idx] = XX[jjj][idx] + E_B[idx][jjj];
// 		}
// 	}// First update to the states and controls due to sensitivities to the states with error as the initial pertubation
// 
// 
// 	double BT[n_ode];
// 	for (int idx = 0;idx < ph0->n_ode;idx++) {
// 		BT[idx] = DXDB[idx][n - 1];
// 	}//Vector of the sensitivity of the final time to the states.
// 
// 	double TF1;
// 	double S = 0;
// 	for (int idx = 0; idx < ph0->n_ode; idx++) {
// 		S = S + E[idx] * BT[idx];
// 	} //Calculates correction for the time by multiplying E and BT
// 
// 	TF1 = TF + S; //first time correction
// 
// 				  //---------------------------Nachkorrektur---------------------------------------------------------------------------------------------	
// 
// 
// 	double X2[n_dis][n_ode + n_ctrl];
// 	for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
// 		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
// 			X2[jjj][idx] = X1[jjj][idx];
// 		}
// 	}// copies the matrix X1 to X2 for Nachkorrektur
// 
// 	double TF2 = TF1; //copies TF1 to TF2 for Nachkorrektur
// 
// 
// 	double dim = ph0->n_dis - 1;
// 
// 	double x1[n_ode + n_ctrl];
// 	double x2[n_ode + n_ctrl];
// 	double dx1[n_ode];
// 	double dx2[n_ode];
// 	double EE[n_ode*(n_dis - 1)];
// 	double BB[n_ode + n_ctrl][n_ode*(n_dis - 1)][n_dis];
// 	double EE_BB[n_ode + n_ctrl][n_dis];
// 	double BBT[n_ode*(n_dis - 1)];
// 
// 	for (int z = 0;z < it;z++) {
// 
// 		double T2 = TF2 / dim; //equals one time step
// 
// 		for (int jjj = 0;jjj < ph0->n_dis - 1;jjj++) {
// 
// 			for (int kkk = 0; kkk < ph0->n_ode + ph0->n_ctrl; kkk++) {
// 				x1[kkk] = X2[jjj][kkk];
// 			}
// 
// 			for (int kkk = 0; kkk < ph0->n_ode + ph0->n_ctrl; kkk++) {
// 				x2[kkk] = X2[jjj + 1][kkk];
// 			}
// 
// 			void(*Funk)(double*, double*) = system;
// 			Funk(x1, dx1);
// 
// 			Funk(x2, dx2);
// 
// 			for (int idx = 0;idx < ph0->n_ode;idx++) {
// 				EE[idx + (ph0->n_ode)*jjj] = -(X2[jjj + 1][idx] - X2[jjj][idx] - (dx1[idx] + dx2[idx])*(T2 / 2.0));
// 			}
// 		}// Calculates error between X2 and integrated X2
// 
// 
// 
// 		for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
// 			for (int idx = 0;idx < ph0->n_ode*(ph0->n_dis - 1);idx++) {
// 				for (int jjj = 0;jjj < ph0->n_dis;jjj++) {
// 					//cout << DXDQ[idx][(k + (ph0->n_ode + ph0->n_ctrl)*jjj)] << endl;
// 					BB[k][idx][jjj] = DXDQ[idx][(k + (ph0->n_ode + ph0->n_ctrl)*jjj)];
// 				}
// 			}
// 		}//Makes a 3D matrix containing the sensitivity of the states and controls to the ODE's, where k is the state or control to be updated
// 
// 		for (int k = 0;k < ph0->n_ode + ph0->n_ctrl;k++) {
// 			for (int idx = 0; idx < ph0->n_dis; idx++) {
// 				EE_BB[k][idx] = 0.0;
// 				for (int jjj = 0; jjj < ph0->n_ode*(ph0->n_dis - 1);jjj++) {
// 					EE_BB[k][idx] = EE_BB[k][idx] + EE[jjj] * BB[k][jjj][idx];
// 				}
// 			}
// 		}// Calculates correction vector for each state and control by multiplying EE and BB
// 
// 		for (int idx = 0; idx < ph0->n_ode + ph0->n_ctrl;idx++) {
// 			for (int jjj = 0; jjj < ph0->n_dis; jjj++) {
// 				X2[jjj][idx] = X2[jjj][idx] + EE_BB[idx][jjj];
// 			}
// 		}// Update to the states and controls due to sensitivities to the ODE's with error as the offset from the trapezoidal rule of integration
// 
// 
// 		for (int idx = 0;idx < ph0->n_ode*(ph0->n_dis - 1);idx++) {
// 			BBT[idx] = DXDQ[idx][n - 1];
// 
// 		}//Vector of the sensitivity of the final time to the ODE's
// 
// 		double S = 0.0;
// 		for (int idx = 0; idx < ph0->n_ode*(ph0->n_dis - 1); idx++) {
// 			S = S + EE[idx] * BBT[idx];
// 		} //Calculates correction for the time by multiplying EE and BBT
// 
// 		TF2 = TF2 + S; //Time correction
// 
// 	} //End Nachkorrektur loop
// 
// 
// 
// 	ofstream f("X2.txt");
// 	for (int jjj = 0; jjj < n_dis; jjj++) {
// 		for (int idx = 0; idx < n_ode + n_ctrl; idx++) {
// 			f << X2[jjj][idx] << endl; //creation of a .txt file with the new trajectory (for use in matlab)
// 			if (idx < n_ode) {
// 				ph0->solver->X[ph0->x_index(jjj, idx)] = X2[jjj][idx]; //writing the new states in ph.solver
// 			}
// 			else {
// 				ph0->solver->X[ph0->u_index(jjj, idx - n_ode)] = X2[jjj][idx]; //writing the new controls in ph.solver
// 			}
// 		}
// 	}
// 	f.close();
// 
// }



