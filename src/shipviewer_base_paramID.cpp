#include "shipviewer_base_paramID.h"

using namespace std;

extern TransWorhpProblem *ph0;

#define nbrShips 1


void shipviewer3d ( glObject *obj, double *x, double t ) {
	int nbrKM = 2;
	int NS[nbrShips+1] = { 0, 10 };
	double tEnd = 600;

	int maxi;

	for ( int k=0; k<nbrShips; k++ ) {
		glPointSize ( 10 );
		glDisable ( GL_LIGHTING );
		glColor3f ( 0.5,0.5,0.5 );


		maxi = ph0->n_dis*t/ ( tEnd /60 );

		if ( maxi>ph0->n_dis ) maxi = ph0->n_dis;

		for ( int i=0; i<maxi; i++ ) {
			glPushMatrix();
			glTranslatef ( ph0->x ( i,NS[k]+1 ) *0.001,ph0->x ( i,NS[k]+0 ) *0.001,0 );

			glBegin ( GL_POINTS );
			glVertex3f ( 0,0,0 );
			glEnd();

			glPopMatrix();
		}
	}



	glLineWidth ( 2.5 );
	glDisable ( GL_LIGHTING );
	glBegin ( GL_LINES );

	glColor3f ( 0.8, 0.2, 0 );
	glVertex3f ( -nbrKM,0,0 );
	glVertex3f ( nbrKM,0,0 );

	glColor3f ( 0, 84./243, 159./243 );
	glVertex3f ( 0,-nbrKM,0 );
	glVertex3f ( 0,nbrKM,0 );

	glEnd();

	glLineWidth ( 1.5 );
	glDisable ( GL_LIGHTING );

	glLineWidth ( 1 );
	glBegin ( GL_LINES );

	for ( int i=-nbrKM; i<nbrKM+1; i++ ) {
//     glLineWidth(1);
		glColor3f ( 0,0,0 );
		glVertex3f ( -nbrKM,i,0 );
		glVertex3f ( nbrKM,i,0 );
	}

	for ( int i=-nbrKM; i<nbrKM+1; i++ ) {
//     glLineWidth(1);
		glColor3f ( 0,0,0 );
		glVertex3f ( i,-nbrKM,0 );
		glVertex3f ( i,nbrKM,0 );
	}

	glEnd();


	glLineWidth ( 0.5 );
	glDisable ( GL_LIGHTING );

//   glLineWidth(1);
	glBegin ( GL_LINES );

	for ( int i=-10*nbrKM; i<10*nbrKM+1; i++ ) {
//     glLineWidth(1);
		glColor3f ( 0,0,0 );
		glVertex3f ( -nbrKM,0.1*i,0 );
		glVertex3f ( nbrKM,0.1*i,0 );
	}

	for ( int i=-10*nbrKM; i<10*nbrKM+1; i++ ) {
//     glLineWidth(1);
		glColor3f ( 0,0,0 );
		glVertex3f ( 0.1*i,-nbrKM,0 );
		glVertex3f ( 0.1*i,nbrKM,0 );
	}

	glEnd();


	for ( int k=0; k<nbrShips; k++ ) {
		glEnable ( GL_LIGHTING );

		Vektor<float> v1 ( 0,0,0 );
		glPushMatrix();

		glTranslatef ( x[NS[k]+1]*0.001,x[NS[k]+0]*0.001,0 );
		glRotatef ( x[NS[k]+2]*180/M_PI,0,0,-1 );

		obj->Draw ( v1,0,0,1 );
		glPopMatrix();
	}

}

