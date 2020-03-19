#include "shipviewer_base_Nships.h"
#include "global.h"


using namespace std;


void shipviewer3d ( tw::glObject *obj, double *x, double t ) {
	
	int nbrKM = 9;
	
	int    *maxi = new    int[ ph0->nbrCShips ];
	double *dt   = new double[ ph0->nbrCShips ];

	//Draw trajectory
	for ( int k=0; k<ph0->nbrCShips; k++ ) {
		glPointSize ( 8 );
		glDisable ( GL_LIGHTING );
		glColor3f ( 0,0,0 );

		if(ph0->ENDTIME == COMMON_AND_FIXED)       maxi[k] = ph0->n_dis*t/ ( ph0->tEnd[0]*ph0->scale_t[0] /ph0->viewer_t_scale );
		if(ph0->ENDTIME == OPEN)                   maxi[k] = ph0->n_dis*t/ ( ph0->p(k)*ph0->scale_t[k] /ph0->viewer_t_scale );
		if(ph0->ENDTIME == OPEN_WITH_FIXED_RATIOS) maxi[k] = ph0->n_dis*t/ ( ph0->ratioToFirst[k]*ph0->p (0)*ph0->scale_t[0] /ph0->viewer_t_scale );
		
		if ( maxi[k]>ph0->n_dis-1 ) maxi[k] = ph0->n_dis-1;

		for ( int i=0; i<=maxi[k]; i++ ) {
			glPushMatrix();
			glTranslatef ( ph0->x ( i,ph0->NS[k]+1 )*0.001*ph0->scale[k][1],ph0->x ( i,ph0->NS[k]+0 )*0.001*ph0->scale[k][0],0 );

			glBegin ( GL_POINTS );
			glVertex3f ( 0,0,0 );
			glEnd();

			glPopMatrix();
		}
	}



	//Draw grid
	
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
	
// 	double shift_y = -5350;
// 	double shift_x = -33400;
	double shift_y = 0;
	double shift_x = 0;
	
// 	Draw harbor
	if(ph0->geom_constr.n_harb_polygons > 0) {
		glLineWidth ( 5. );
		glDisable ( GL_LIGHTING );
		glBegin ( GL_LINES );
		
		std::vector<Vec2d>* harbor_segments=ph0->geom_constr.getHarborSegements(0);

		for (unsigned int i=0; i<harbor_segments->size()-1; i++) {
			glColor3f ( 0,0,0 );
			glVertex3f ( (shift_y+(*harbor_segments)[i].y)/1000., (shift_x+(*harbor_segments)[i].x)/1000.,0 );
			glVertex3f ( (shift_y+(*harbor_segments)[i+1].y)/1000., (shift_x+(*harbor_segments)[i+1].x)/1000.,0 );
		}
		glColor3f ( 0,0,0 );
		glVertex3f ( (shift_y+(*harbor_segments)[harbor_segments->size()-1].y)/1000., (shift_x+(*harbor_segments)[harbor_segments->size()-1].x)/1000.,0 );
		glVertex3f ( (shift_y+(*harbor_segments)[0].y)/1000., (shift_x+(*harbor_segments)[0].x)/1000.,0 );
		
		glEnd();
	}
	
	//Draw obstacles	
	for (unsigned int j=0; j<ph0->geom_constr.getNbrObstacles(); j++) {
		
		std::vector<Vec2d>* obstacle_segments=ph0->geom_constr.getObstaclesSegements(j);
		unsigned int n_segments = obstacle_segments->size();
	
		glLineWidth ( 5. );
		glDisable ( GL_LIGHTING );
		glBegin ( GL_LINES );
		
		for (unsigned int i=0; i<n_segments-1; i++) {
			glColor3f ( 0,0,0 );
			glVertex3f ( (*obstacle_segments)[i].y/1000.,(*obstacle_segments)[i].x/1000.,0 );
			glVertex3f ( (*obstacle_segments)[i+1].y/1000., (*obstacle_segments)[i+1].x/1000.,0 );
		}
		glColor3f ( 0,0,0 );
		glVertex3f ( (*obstacle_segments)[n_segments-1].y/1000.,
					 (*obstacle_segments)[n_segments-1].x/1000.,0 );
		glVertex3f ( (*obstacle_segments)[0].y/1000., (*obstacle_segments)[0].x/1000.,0 );
		
		glEnd();
	
	}


// // 	Draw ships
// 	for ( int k=0; k<ph0->nbrCShips; k++ ) {
// 		glEnable ( GL_LIGHTING );
// 
// 		Vektor<float> v1 ( 0,0,0 );
// 		glPushMatrix();
// 
// 		glTranslatef ( (shift_y+x[ph0->NS[k]+1])*0.001,(shift_x+x[ph0->NS[k]+0])*0.001,0 );
// 		glRotatef ( x[ph0->NS[k]+2]*180/M_PI,0,0,-1 );
// 
// 		obj->Draw ( v1,0,0,1 );
// 		glPopMatrix();
// 	}

// // 	Draw AIS ships
// 	for ( int k=0; k<ph0->nbrAShips; k++ ) {
// 		glEnable ( GL_LIGHTING );
// 
// 		Vektor<float> v1 ( 0,0,0 );
// 		glPushMatrix();
// 
// 		glTranslatef ( (shift_y+ph0->AISShips[k].y0)*0.001,(shift_x+ph0->AISShips[k].x0)*0.001,0 );
// 		glRotatef ( ph0->AISShips[k].psi*180/M_PI,0,0,-1 );
// 
// 		obj->Draw ( v1,0,0,1 );
// 		glPopMatrix();
// 	}
	
	
// // 	Draw ships
	double x0, y0, psi;
// 
	for ( int k=0; k<ph0->nbrAShips; k++ ) {
		glEnable ( GL_LIGHTING );
		
		Vektor<float> v1 ( 0,0,0 );
		glPushMatrix();
		
		double time = t*ph0->viewer_t_scale;
		
		x0 = ph0->AISShips[k].x0 + time*ph0->AISShips[k].sog*cos(ph0->AISShips[k].cog);
		y0 = ph0->AISShips[k].y0 + time*ph0->AISShips[k].sog*sin(ph0->AISShips[k].cog);
		psi = ph0->AISShips[k].psi;
		
		glTranslatef ( y0*0.001,x0*0.001,0 );
		glRotatef ( psi*180./M_PI,0,0,-1 );

		obj->Draw ( v1,0,0,1 );
		glPopMatrix();
	}
	
	for ( int k=0; k<ph0->nbrCShips; k++ ) {
		glEnable ( GL_LIGHTING );
		
		Vektor<float> v1 ( 0,0,0 );
		glPushMatrix();
		
		if(ph0->ENDTIME == COMMON_AND_FIXED)       maxi[k] = ph0->n_dis*t/ ( ph0->tEnd[0]*ph0->scale_t[0] /ph0->viewer_t_scale );
		if(ph0->ENDTIME == OPEN)                   maxi[k] = ph0->n_dis*t/ ( ph0->p(k)*ph0->scale_t[k] /ph0->viewer_t_scale );
		if(ph0->ENDTIME == OPEN_WITH_FIXED_RATIOS) maxi[k] = ph0->n_dis*t/ ( ph0->ratioToFirst[k]*ph0->p(0)*ph0->scale_t[0] /ph0->viewer_t_scale );
		if ( maxi[k]>ph0->n_dis-2 ) maxi[k] = ph0->n_dis-2;
		if(ph0->ENDTIME == COMMON_AND_FIXED)       dt[k] = ph0->n_dis*t/ ( ph0->tEnd[0]*ph0->scale_t[0] /ph0->viewer_t_scale ) - maxi[k];
		if(ph0->ENDTIME == OPEN)                   dt[k] = ph0->n_dis*t/ ( ph0->p ( k )*ph0->scale_t[k] /ph0->viewer_t_scale ) - maxi[k];
		if(ph0->ENDTIME == OPEN_WITH_FIXED_RATIOS) dt[k] = ph0->n_dis*t/ ( ph0->ratioToFirst[k]*ph0->p(0)*ph0->scale_t[0] /ph0->viewer_t_scale ) - maxi[k];
		
		if(dt[k]<=1) {
			if(ph0->ENDTIME == COMMON_AND_FIXED) {
				cubicSpline(maxi[k], k, dt[k], ph0->tEnd[0]*ph0->scale_t[0]/(ph0->n_dis-1), 0, &x0, &y0, &psi);
			}
			if(ph0->ENDTIME == OPEN) {
				cubicSpline(maxi[k], k, dt[k], ph0->p(k)*ph0->scale_t[k]/(ph0->n_dis-1), 0, &x0, &y0, &psi);
			}
			if(ph0->ENDTIME == OPEN_WITH_FIXED_RATIOS) {
				cubicSpline(maxi[k], k, dt[k], ph0->ratioToFirst[k]*ph0->p(0)*ph0->scale_t[0]/ph0->n_dis, 0, &x0, &y0, &psi);
			}
		}
		if(dt[k]>1) {
			maxi[k] = ph0->n_dis-1;
			if(ph0->ENDTIME == COMMON_AND_FIXED) {
				cubicSpline(maxi[k], k, dt[k]-1, ph0->tEnd[0]*ph0->scale_t[0]/ph0->n_dis, 1, &x0, &y0, &psi);
			}
			if(ph0->ENDTIME == OPEN) {
				cubicSpline(maxi[k], k, dt[k]-1, ph0->p(k)*ph0->scale_t[k]/ph0->n_dis, 1, &x0, &y0, &psi);
			}
			if(ph0->ENDTIME == OPEN_WITH_FIXED_RATIOS) {
				cubicSpline(maxi[k], k, dt[k]-1, ph0->ratioToFirst[k]*ph0->p(0)*ph0->scale_t[0]/ph0->n_dis, 1, &x0, &y0, &psi);
			}
		}

		glTranslatef ( y0*0.001*ph0->scale[k][1],x0*0.001*ph0->scale[k][0],0 );
		glRotatef ( psi*180/M_PI*ph0->scale[k][2],0,0,-1 );

		obj->Draw ( v1,0,0,1 );
		glPopMatrix();
	}
	
		delete[] maxi;
		delete[] dt;

}

