/*----------------------------------------------------------------
 *
 * Flugzeug: Einfache Trajektorie
 *
 *----------------------------------------------------------------*/

#define _USE_MATH_DEFINES
#include <cmath>

#ifdef WIN32
#include "windows.h"
#endif

#include "TransWORHP.h"


using namespace std;

/** Modell von Bittner, Matthias <m.bittner@tum.de> */
typedef struct {
	const double MatrixConcatenate[3];
} ConstB_Model3Dof_T;

const ConstB_Model3Dof_T Model3Dof_ConstB = {
	{ 1.0, 0.0, 0.0 }
};

void Model3Dof(const double Model3Dof_U_States[8],
               const double Model3Dof_U_ControlsCMD[4],
               double Model3Dof_Y_StatesDot[8],
               double CON[1]) {
	double rtb_N_Z_B_UB;
	double rtb_Mue_K_UB;
	double rtb_TrigonometricFunction;
	double rtb_TrigonometricFunction2_g;
	double rtb_Transpose3;
	double rtb_Product6_b[3];
	double tmp[9];
	int i;
	double rtb_Sum_0;
	double rtb_MatrixConcatenate1_e_idx;
	double rtb_MatrixConcatenate2_c_idx;
	double rtb_Sum1_idx;

	/* Trigonometry: '<S8>/Trigonometric Function4' incorporates:
	 *  Inport: '<Root>/States'
	 */
	rtb_N_Z_B_UB = cos(Model3Dof_U_States[4]);

	/* Product: '<S8>/Product1' incorporates:
	 *  Inport: '<Root>/States'
	 *  Trigonometry: '<S8>/Trigonometric Function2'
	 */
	rtb_Mue_K_UB = Model3Dof_U_States[6] * cos(Model3Dof_U_States[3]) *
	               rtb_N_Z_B_UB;

	/* Product: '<S8>/Product2' incorporates:
	 *  Inport: '<Root>/States'
	 *  Trigonometry: '<S8>/Trigonometric Function1'
	 */
	rtb_N_Z_B_UB *= Model3Dof_U_States[6] * sin(Model3Dof_U_States[3]);

	/* Trigonometry: '<S17>/Trigonometric Function1' incorporates:
	 *  Inport: '<Root>/States'
	 */
	rtb_TrigonometricFunction = cos(Model3Dof_U_States[5]);

	/* Trigonometry: '<S17>/Trigonometric Function' incorporates:
	 *  Inport: '<Root>/States'
	 */
	rtb_TrigonometricFunction2_g = sin(Model3Dof_U_States[5]);

	/* Gain: '<S17>/Gain1' */
	rtb_MatrixConcatenate1_e_idx = -rtb_TrigonometricFunction2_g;

	/* Gain: '<S17>/Gain2' */
	rtb_MatrixConcatenate2_c_idx = rtb_TrigonometricFunction2_g;

	/* Product: '<S12>/Product4' incorporates:
	 *  Constant: '<S12>/C_L_Alpha5'
	 *  Inport: '<Root>/ControlsCMD'
	 */
	rtb_TrigonometricFunction2_g = Model3Dof_U_ControlsCMD[0] * 4.75;

	/* Math: '<S12>/Transpose3' */
	rtb_Transpose3 = rtb_TrigonometricFunction2_g * rtb_TrigonometricFunction2_g;

	/* Product: '<S12>/Product7' incorporates:
	 *  Constant: '<S12>/C_L_Alpha7'
	 *  Inport: '<Root>/ControlsCMD'
	 */
	rtb_TrigonometricFunction2_g = Model3Dof_U_ControlsCMD[1] * -0.589355;

	/* Gain: '<S5>/Gain1' incorporates:
	 *  Constant: '<S12>/C_L_Alpha1'
	 *  Constant: '<S12>/C_L_Alpha3'
	 *  Constant: '<S12>/C_L_Alpha4'
	 *  Math: '<S12>/Transpose1'
	 *  Product: '<S12>/Product1'
	 *  Product: '<S12>/Product3'
	 *  Sum: '<S12>/Sum2'
	 *  Sum: '<S12>/Sum3'
	 */
	rtb_Product6_b[0] = -((rtb_TrigonometricFunction2_g *
	                       rtb_TrigonometricFunction2_g * 1.6967701979282437 + rtb_Transpose3 * 0.05134)
	                      + 0.0761);

	/* Product: '<S14>/Product1' incorporates:
	 *  Constant: '<S14>/C_L_Alpha4'
	 *  Inport: '<Root>/ControlsCMD'
	 */
	rtb_Product6_b[1] = Model3Dof_U_ControlsCMD[1] * -0.589355;

	/* Gain: '<S5>/Gain4' incorporates:
	 *  Constant: '<S13>/C_L_Alpha10'
	 *  Constant: '<S13>/C_L_Alpha9'
	 *  Inport: '<Root>/ControlsCMD'
	 *  Product: '<S13>/Product5'
	 *  Sum: '<S13>/Sum'
	 */
	rtb_Product6_b[2] = -(Model3Dof_U_ControlsCMD[0] * 4.75 + 0.055);

	/* Product: '<S5>/Product1' incorporates:
	 *  Constant: '<S5>/Constant3'
	 *  Constant: '<S5>/Density_Const'
	 *  Inport: '<Root>/States'
	 *  Math: '<S5>/Transpose1'
	 */
	rtb_TrigonometricFunction2_g = Model3Dof_U_States[6] * Model3Dof_U_States[6] *
	                               0.5 * 1.225;

	/* Saturate: '<S5>/Saturation2' incorporates:
	 *  Inport: '<Root>/States'
	 */
	if (Model3Dof_U_States[6] >= 0.05) {
		rtb_Transpose3 = Model3Dof_U_States[6];
	} else {
		rtb_Transpose3 = 0.05;
	}

	/* Sum: '<S5>/Sum1' incorporates:
	 *  Constant: '<S5>/ACC_Gravity_1'
	 *  Constant: '<S5>/ACC_Gravity_2'
	 *  Constant: '<S5>/S'
	 *  Inport: '<Root>/ControlsCMD'
	 *  Product: '<S5>/Product3'
	 *  Product: '<S5>/Product4'
	 *  Product: '<S5>/Product6'
	 *  Saturate: '<S5>/Saturation2'
	 */
	rtb_Transpose3 = rtb_Product6_b[0] * 8.928 * rtb_TrigonometricFunction2_g /
	                 6796.0084499999994 + Model3Dof_U_ControlsCMD[3] * 30.0 * 0.8 /
	                 rtb_Transpose3;
	rtb_Sum1_idx = rtb_Product6_b[1] * 8.928 * rtb_TrigonometricFunction2_g /
	               6796.0084499999994;
	rtb_Sum_0 = rtb_Product6_b[2] * 8.928 * rtb_TrigonometricFunction2_g /
	            6796.0084499999994;

	/* Concatenate: '<S17>/Matrix Concatenate3' incorporates:
	 *  Product: '<S6>/Product6'
	 *  SignalConversion: '<S17>/ConcatBufferAtMatrix Concatenate1In2'
	 *  SignalConversion: '<S17>/ConcatBufferAtMatrix Concatenate2In3'
	 */
	tmp[0] = Model3Dof_ConstB.MatrixConcatenate[0];
	tmp[3] = Model3Dof_ConstB.MatrixConcatenate[1];
	tmp[6] = Model3Dof_ConstB.MatrixConcatenate[2];
	tmp[1] = 0.0;
	tmp[4] = rtb_TrigonometricFunction;
	tmp[7] = rtb_MatrixConcatenate1_e_idx;
	tmp[2] = 0.0;
	tmp[5] = rtb_MatrixConcatenate2_c_idx;
	tmp[8] = rtb_TrigonometricFunction;

	/* Product: '<S6>/Product6' incorporates:
	 *  SignalConversion: '<S6>/ConcatBufferAtMatrix ConcatenateIn1'
	 *  SignalConversion: '<S6>/ConcatBufferAtMatrix ConcatenateIn2'
	 */
	for (i = 0; i < 3; i++) {
		rtb_Product6_b[i] = 0.0;
		rtb_Product6_b[i] += tmp[i] * rtb_Transpose3;
		rtb_Product6_b[i] += tmp[i + 3] * rtb_Sum1_idx;
		rtb_Product6_b[i] += tmp[i + 6] * rtb_Sum_0;
	}

	/* Saturate: '<S11>/Saturation2' incorporates:
	 *  Inport: '<Root>/States'
	 */
	if (Model3Dof_U_States[6] >= 0.05) {
		rtb_TrigonometricFunction2_g = Model3Dof_U_States[6];
	} else {
		rtb_TrigonometricFunction2_g = 0.05;
	}

	/* End of Saturate: '<S11>/Saturation2' */

	/* Trigonometry: '<S11>/Trigonometric Function2' incorporates:
	 *  Inport: '<Root>/States'
	 */
	rtb_Transpose3 = cos(Model3Dof_U_States[4]);

	/* Saturate: '<S11>/Saturation1' */
	if (rtb_Transpose3 >= 1.0) {
		rtb_Transpose3 = 1.0;
	} else {
		if (rtb_Transpose3 <= 0.05) {
			rtb_Transpose3 = 0.05;
		}
	}

	/* Product: '<S11>/Divide' incorporates:
	 *  Product: '<S11>/Product2'
	 *  Saturate: '<S11>/Saturation1'
	 *  SignalConversion: '<S11>/ConcatBufferAtMatrix ConcatenateIn2'
	 */
	rtb_TrigonometricFunction = 1.0 / rtb_Transpose3 * (rtb_Product6_b[1] *
	                            9.80665) / rtb_TrigonometricFunction2_g;

	/* Product: '<S11>/Divide1' incorporates:
	 *  Gain: '<S11>/Gain6'
	 *  Inport: '<Root>/States'
	 *  Product: '<S11>/Product2'
	 *  SignalConversion: '<S11>/ConcatBufferAtMatrix ConcatenateIn3'
	 *  Sum: '<S11>/Sum'
	 *  Trigonometry: '<S11>/Trigonometric Function1'
	 */
	rtb_TrigonometricFunction2_g = -((cos(Model3Dof_U_States[4]) + rtb_Product6_b
	                                  [2]) * 9.80665) * (1.0 / rtb_TrigonometricFunction2_g);

	/* Outport: '<Root>/StatesDot' incorporates:
	 *  Constant: '<S4>/ACC_Gravity_3'
	 *  Gain: '<S11>/Gain2'
	 *  Inport: '<Root>/ControlsCMD'
	 *  Inport: '<Root>/States'
	 *  Product: '<S11>/Product2'
	 *  Product: '<S4>/Product2'
	 *  Product: '<S8>/Product3'
	 *  SignalConversion: '<S11>/ConcatBufferAtMatrix ConcatenateIn1'
	 *  Sum: '<S11>/Sum'
	 *  Sum: '<S4>/Add1'
	 *  Trigonometry: '<S11>/Trigonometric Function3'
	 *  Trigonometry: '<S8>/Trigonometric Function3'
	 */
	Model3Dof_Y_StatesDot[0] = rtb_Mue_K_UB;
	Model3Dof_Y_StatesDot[1] = rtb_N_Z_B_UB;
	Model3Dof_Y_StatesDot[2] = -(Model3Dof_U_States[6] * sin(Model3Dof_U_States[4]));
	Model3Dof_Y_StatesDot[3] = rtb_TrigonometricFunction;
	Model3Dof_Y_StatesDot[4] = rtb_TrigonometricFunction2_g;
	Model3Dof_Y_StatesDot[5] = Model3Dof_U_ControlsCMD[2];
	Model3Dof_Y_StatesDot[6] = (-sin(Model3Dof_U_States[4]) + rtb_Product6_b[0]) *
	                           9.80665;
	Model3Dof_Y_StatesDot[7] = (Model3Dof_U_ControlsCMD[3] - Model3Dof_U_States[7])
	                           / 0.5;

	/* Outport: '<Root>/Outputs' incorporates:
	 *  Constant: '<S7>/Mue_K_LB'
	 *  Constant: '<S7>/Mue_K_UB'
	 *  Constant: '<S7>/N_Z_B_LB'
	 *  Constant: '<S7>/N_Z_B_UB'
	 *  Constant: '<S7>/POS_Z_LB'
	 *  Constant: '<S7>/VEL_ABS_K_LB'
	 *  Constant: '<S7>/VEL_ABS_K_UB'
	 *  Gain: '<S7>/Gain1'
	 *  Gain: '<S7>/Gain6'
	 *  Inport: '<Root>/States'
	 *  Sum: '<S7>/Sum1'
	 *  Sum: '<S7>/Sum2'
	 *  Sum: '<S7>/Sum3'
	 *  Sum: '<S7>/Sum4'
	 *  Sum: '<S7>/Sum5'
	 *  Sum: '<S7>/Sum6'
	 *  Sum: '<S7>/Sum7'
	 */
	/* Model3Dof_Y_Outputs[0] = -Model3Dof_U_States[2] - -10.0;
	 Model3Dof_Y_Outputs[1] = Model3Dof_U_States[6] - 25.0;
	 Model3Dof_Y_Outputs[2] = 102.9 - Model3Dof_U_States[6];
	 Model3Dof_Y_Outputs[3] = Model3Dof_U_States[5] - -1.8325957145940461;
	 Model3Dof_Y_Outputs[4] = 1.8325957145940461 - Model3Dof_U_States[5];
	 Model3Dof_Y_Outputs[5] = -rtb_Sum_0 - -2.0;
	 Model3Dof_Y_Outputs[6] = 12.0 - (-rtb_Sum_0);
	 */

	CON[0] = rtb_Sum_0;
}


/** globaler Zeiger auf OCP für Grafik */
TransWorhp *pph=0;

double flugzeugplot (int &len, int &ndgl, int &nsteuer, double *t,double *x, double *u, int &i, int &index) {

	if (index==0) {
		int g_index = (pph->n_dis-1) * pph->n_ode + pph->n_rand + i*pph->n_neben;
		return pph->G[g_index + 0];
	}

	return 0;
}


void flugzeug3d(glObject *obj, double *x, double t) {

	double scale = 100;
	
	glLineWidth(2);
	glColor3f(.8,.1,.1);
	glDisable(GL_LIGHTING);
	glBegin(GL_LINE_STRIP);
	
	
	for (int i=0;i<pph->n_dis;i++)
		glVertex3f((pph->x(i,0)-400)/scale,pph->x(i,1)/scale,pph->x(i,2)/scale*30 );
	
	glEnd();
	
	
	glEnable(GL_LIGHTING);
	
	float offset = 4;

	int n = obj->countObj();

	Vektor<float> v(0,0,0);
	
	glColor3f(153/255.,185/255.,131/255.);
	obj->Draw(v,0,0,2);
	
	obj->Draw(v,0,0,1);
	
	
	glPushMatrix();
	glTranslatef((x[0]-400)/scale, x[1]/scale, x[2]/scale*30);
	
	glRotatef(x[3]*180/M_PI,0,0,1);
	glRotatef(10* x[4]*180/M_PI,0,1,0);
	glRotatef(-90,0,0,1);
	glColor3f(1,1,1);
	
	obj->Draw(v,0,0,4);
	glPopMatrix();
	/*

	Vektor<float> v(0,0,0);

	obj->Draw(v,0,0,1); // Ebene
	obj->Draw(v,0,0,2); // Regal1
	obj->Draw(v,0,0,3); // Regal2


	glPushMatrix();
	glTranslatef(x[0], 0, 0);
	obj->Draw(v,0,0,4);
	obj->Draw(v,0,0,5);

	obj->Draw(v,0,0,8);
	obj->Draw(v,0,0,9);
	obj->Draw(v,0,0,10);
	obj->Draw(v,0,0,11);

	glPopMatrix();

	glPushMatrix();
	glTranslatef(x[0]-x[2], 0, 9 - x[5]);
	obj->Draw(v,0,0,6);
	obj->Draw(v,0,0,7);
	glPopMatrix();


	// Seile
	glLineWidth(2);
	glColor3f(.3,.3,.9);
	glDisable(GL_LIGHTING);

	float xx=1.2, yy=1;

	glBegin(GL_LINES);
	glVertex3f(x[0] +xx, +yy, 9.2);
	glVertex3f(x[0]-x[2]+xx, +yy, 9-x[5]);

	glVertex3f(x[0] -xx, +yy, 9.2);
	glVertex3f(x[0]-x[2]-xx, +yy, 9-x[5]);

	glVertex3f(x[0] +xx, -yy, 9.2);
	glVertex3f(x[0]-x[2]+xx, -yy, 9-x[5]);

	glVertex3f(x[0] -xx, -yy, 9.2);
	glVertex3f(x[0]-x[2]-xx, -yy, 9-x[5]);
	glEnd();*/


}


class FlugzeugPhase : public TransWorhp {
public:

	XMLNode* scenenode;
	
	FlugzeugPhase(XMLNode *xmlmain, int dis) : TransWorhp("Flugzeug", dis,8,4,1,0,1,1) {
		/** freie Endzeit */
		freetime = 1;
		
		/** Gewichtung fuer Integral-Anteil der Zielfunktion */
		lagrange_weight[0] = .01;
		
		if (xmlmain)
			scenenode = xmlmain->GetFirstChild("SCENE");
	}

	void OpenWindows(Viewer *viewer) {

		viewer->ThreeD("Flugbahn", scenenode, flugzeug3d);
		viewer->Data("Lastfaktor", flugzeugplot,0);

	}
	void GetXTitle(int d, char *s) {
		if (d==0) strcpy(s, "Position X");
		if (d==1) strcpy(s, "Position Y");
		if (d==2) strcpy(s, "Position Z");
		if (d==3) strcpy(s, "Kurswinkel chi");
		if (d==4) strcpy(s, "Steigwinkel gamma");
		if (d==5) strcpy(s, "Haengewinkel mu");
		if (d==6) strcpy(s, "Geschwindigkeit V");
		if (d==7) strcpy(s, "Schubhebelstellung delta_t");
	}

	void GetUTitle(int d, char *s) {
		if (d==0) strcpy(s, "Anstellwinkel alpha_CMD");
		if (d==1) strcpy(s, "Schiebewinkel beta_CMD");
		if (d==2) strcpy(s, "Haengewinkelrate dot_mu_CMD");
		if (d==3) strcpy(s, "Schubhebelvorgabe delta_t_CMD");

	}

	void p_init(double *p) {
		/** Startschaetzung der Flugzeit */ 
		p[0] = 12;
	}

	void x_init(double *x, int i, int dis) {
		// x[5] = 5;
	}
	void u_init(double *u, int i, int dis) {
		// u[0] = 0;
		// u[1] = 0;
	}


	double obj() {
		/** Mayer-Anteil der Zielfunktion */
		return 1 * p(0);
	}


	bool obj_structure(DiffStructure &s) {
		s(0,p_index(0));
		return true;
	}

	bool obj_diff(DiffStructure &s) {
		s(0,p_index(0)) = 1;
		return true;
	}


	void integral(double *f, double t, const double *x, const double *u,
	              const double *p) {
		/** Lagrange-Anteil der Zielfunktion */
		f[0] = (u[0]*u[0] + u[1]*u[1] + u[2]*u[2] + u[3]*u[3]) * p[0];
	}

	bool integral_structure(DiffStructure &s) {

		s(0, u_indexode(0));
		s(0, u_indexode(1));
		s(0, u_indexode(2));
		s(0, u_indexode(3));
		s(0, p_indexode(0));
		return true;
	}

	bool integral_diff(DiffStructure &s, double t, const double *x, const double *u,
	                   const double *p) {

		s(0, u_indexode(0)) = 2*u[0]*p[0];
		s(0, u_indexode(1)) = 2*u[1]*p[0];
		s(0, u_indexode(2)) = 2*u[2]*p[0];
		s(0, u_indexode(3)) = 2*u[3]*p[0];
		s(0, p_indexode(0)) = u[0]*u[0] + u[1]*u[1]  + u[2]*u[2]  + u[3]*u[3];
		return true;
	}

	void ode(double *dx, double t, const double *x, const double *u, const double *p) {

		double con;
		Model3Dof(x, u, dx, &con);

		for (int i=0;i<8;i++) {
			dx[i] *= p[0];
		}
	}


	// Optional
	bool ode_structure(DiffStructure &s) {

		for (int i=0;i<8;i++) {

			s(i,p_indexode(0));
		}

		for (int i=1;i<8;i++) {
			if (i==7) continue;
			if (i==5) continue;
			if (i==2) continue;
			if (i==4) continue;
		
			s(i,x_indexode(4));
			s(i,x_indexode(5));
			s(i,x_indexode(6));
			s(i,x_indexode(7));
			s(i,u_indexode(0));
			s(i,u_indexode(1));
			s(i,u_indexode(2));
			s(i,u_indexode(3));
		}

		s(0,x_indexode(3));
		s(0,x_indexode(4));
		s(0,x_indexode(6));
		s(0,p_indexode(0));

		s(1,x_indexode(3));
		s(1,x_indexode(4));
		s(1,x_indexode(6));
		s(1,p_indexode(0));

		s(2,x_indexode(4));
		s(2,x_indexode(6));
		s(2,p_indexode(0));

		s(3,x_indexode(4));
		s(3,x_indexode(5));
		s(3,x_indexode(6));
		s(3,u_indexode(0));
		s(3,u_indexode(1));
		s(3,p_indexode(0));

		s(4,x_indexode(4));
		s(4,x_indexode(5));
		s(4,x_indexode(6));
		s(4,u_indexode(0));
		s(4,u_indexode(1));
		s(4,p_indexode(0));

		s(5,u_indexode(2));
		s(5,p_indexode(0));

		s(6,x_indexode(4));
		s(6,x_indexode(6));
		s(6,u_indexode(0));
		s(6,u_indexode(1));
		s(6,u_indexode(3));
		s(6,p_indexode(0));

		s(7,x_indexode(7));
		s(7,u_indexode(3));
		s(7,p_indexode(0));

		return true;
	}

	bool ode_diff(DiffStructure& s, double t, const double* x, const double* u, const double* p) {
	
		return false;
		
	}

	bool ode_diff_p(DiffStructure& s, double t, const double *x, const double *u, const double *p, int index) {
	
		double con;
		double dx[8];
		Model3Dof(x, u, dx, &con);

		for (int i=0;i<8;i++) {
			s(i,p_indexode(0)) = dx[i];
		}
		
		return true;
	}


	void u_boundary(double *u_low, double *u_upp) {

		/*u_low[0] = -1;
		u_upp[0] =  1;
		u_low[1] = -1;
		u_upp[1] =  1;*/

	}

	void x_boundary(double *x_low, double *x_upp) {

		//x_low[0] = 0;
		//x_upp[0] = 100;

		//x_low[1] = -3;
		//x_upp[1] = 3;

		x_low[2] = 0; // Additional!!!
		x_upp[2] = 10;

		//x_low[3] = -10;
		//x_upp[3] = 10;

		//x_low[4] = -4;
		//x_upp[4] = 4;

		x_low[5] = -1.8325957145940461;
		x_upp[5] = 1.8325957145940461;

		x_low[6] = 25;
		x_upp[6] = 102.9;

		//x_low[7] = -10;
		//x_upp[7] = 10;

		//x_low[8] = 0;

	}

	void p_boundary(double *p_low, double *p_upp) {
		p_low[0] = 1;
		p_upp[0] = 350;
	}

	void var_boundary(double *x_low, double *x_upp) {

		double start[] = {0,0,5,0,0,0,80,0};
		double ziel[] =  {1000,200,5,0,0,0,80,0};

		for (int i=0;i<8;i++) {
			x_low[x_index(0,i)] = start[i];
			x_upp[x_index(0,i)] = start[i];
		}
		for (int i=0;i<8;i++) {
			x_low[x_index(n_dis-1,i)] = ziel[i];
			x_upp[x_index(n_dis-1,i)] = ziel[i];
		}

	}


	void neben_boundary(double *c_low, double *c_upp) {

		c_low[0] = -12;
		c_upp[0] = +2;

	}

	void neben(double *c, double t, const double *x, const double *u, const double *p) {

		double dx[8];
		Model3Dof(x, u, dx, c);

	}

	bool neben_structure(DiffStructure &s) {

		s(0,x_index(0,8));
		s(0,x_index(0,6));
	//	s(0,p_indexode(0));

		return true;
	}
	/*
	bool neben_diff(DiffStructure &s, double t, const double *x, const double *u, const double *p) {

		s(0,x_index(0,4)) =1 ;
		s(0,x_index(0,5)) =-1 ;
		//s.use(0,p_indexode(0));

		return true;
	}

	bool neben_diff_p(DiffStructure &s, double t, const double *x, const double *u, const double *p, int index) {

		s(0,p_indexode(0)) = .01;
		return true;
	}*/

};


/////////////////////////////////////////////////////////////////////////////



int main(int argv, char* argc[]) {

	TWparameter twparameter("transworhp.xml");
	twparameter.NDIS = 41;
	twparameter.Arguments(argv,argc);
	TWfolder folder(&twparameter,0);

	XMLNode *xml_laufkatze = TWparameter::ReadParams("flugzeug.xml");

	FlugzeugPhase ph(xml_laufkatze, twparameter.NDIS);
	folder.Add(&ph);

	pph=&ph;

	Viewer *viewer = 0;
	if (twparameter.PLOT) viewer = new Viewer(&twparameter);

	folder.Init();
	folder.Init(viewer);

	// Optional: Startschaetzung
	double start[] = {0,0,5,0,0,0,80,0};
	ph.X[ph.x_index(0,2)]=5;
	ph.X[ph.x_index(0,6)]=80;

	ph.Integrate(twparameter.butchertableau);
	// Ende Optional

	folder.Loop(2);

	delete viewer;
	delete xml_laufkatze;

	return 0;
}

