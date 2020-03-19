/*    Optimal Trajectory Calculations     */

#ifdef WIN32
#include "windows.h"
#endif
#include "TransWORHP.h"
#include "shipviewer_base_paramID.h"
#include "Models/Model1v4.h"
#include <vector>
#include "TWfolder.h"
#include "TWproblem.h"



using namespace std;

TransWorhpProblem *ph0=0;

int N = 12;
string params ( "ParamXML/s1v1_m1v4_p2.xml" );
vector<const char*> slog = { 
// 									"SLOG/Maneuvers/Maneuvers_P02.slog",
// 									"SLOG/Maneuvers/Maneuvers_P03.slog",
// 									"SLOG/Maneuvers/Maneuvers_P04.slog",
// 									"SLOG/Maneuvers/Maneuvers_P05.slog",
// 									"SLOG/Maneuvers/Maneuvers_P06.slog",
	
							"SLOG/Maneuvers/Maneuvers_P07.slog",
									"SLOG/Maneuvers/Maneuvers_P08.slog",
							"SLOG/Maneuvers/Maneuvers_P09.slog",
									"SLOG/Maneuvers/Maneuvers_P10.slog",
							"SLOG/Maneuvers/Maneuvers_P11.slog",
							
							"SLOG/Maneuvers/Maneuvers_P12.slog",
							"SLOG/Maneuvers/Maneuvers_P13.slog",
									"SLOG/Maneuvers/Maneuvers_P14.slog",
							"SLOG/Maneuvers/Maneuvers_P15.slog",
							"SLOG/Maneuvers/Maneuvers_P16.slog",
							
// 							"SLOG/Maneuvers/Maneuvers_P17.slog",
// 									"SLOG/Maneuvers/Maneuvers_P18.slog",
// 							"SLOG/Maneuvers/Maneuvers_P19.slog",
// 							"SLOG/Maneuvers/Maneuvers_P20.slog",
							"SLOG/Maneuvers/Maneuvers_P21.slog",
							
// 							"SLOG/Maneuvers/Maneuvers_P24.slog",
// 							"SLOG/Maneuvers/Maneuvers_P25.slog",
// 									"SLOG/Maneuvers/Maneuvers_P26.slog",
							"SLOG/Maneuvers/Maneuvers_P27.slog",
// 									"SLOG/Maneuvers/Maneuvers_P31.slog",


// 									"SLOG/Maneuvers/Maneuvers_P22.slog",
// 									"SLOG/Maneuvers/Maneuvers_P23.slog",
// 									"SLOG/Maneuvers/Maneuvers_P30.slog",
// 									"SLOG/Maneuvers/Maneuvers_P33.slog",
// 									"SLOG/Maneuvers/Maneuvers_P35.slog",
	
// 							"SLOG/Maneuvers/Maneuvers_P28.slog", //with thruster
// 							"SLOG/Maneuvers/Maneuvers_P29.slog", //with thruster
// 							"SLOG/Maneuvers/Maneuvers_P32.slog", //error in SIMOPT when generating this data
// 							"SLOG/Maneuvers/Maneuvers_P34.slog", //error in SIMOPT when generating this data


};

// double tEnd = 601;

//CopyFile is a simple function that copies a file from arg1 to arg2
int CopyFile(string initialFilePath, string outputFilePath)
{
	
	ifstream initialFile(initialFilePath.c_str(), ios::in|ios::binary);
	ofstream outputFile(outputFilePath.c_str(), ios::out|ios::binary);
	//defines the size of the buffer
	initialFile.seekg(0, ios::end);
	long fileSize = initialFile.tellg();
	//Requests the buffer of the predefined size


	//As long as both the input and output files are open...
	if(initialFile.is_open() && outputFile.is_open())
	{
		short * buffer = new short[fileSize/2];
		//Determine the file's size
		//Then starts from the beginning
		initialFile.seekg(0, ios::beg);
		//Then read enough of the file to fill the buffer
		initialFile.read((char*)buffer, fileSize);
		//And then write out all that was read
		outputFile.write((char*)buffer, fileSize);
		delete[] buffer;
	}
	
	//If there were any problems with the copying process, let the user know
	else if(!outputFile.is_open())
	{
		cout<<"I couldn't open "<<outputFilePath<<" for copying!\n";
		return 0;
	}
	else if(!initialFile.is_open())
	{
		cout<<"I couldn't open "<<initialFilePath<<" for copying!\n";
		return 0;
	}
		
	initialFile.close();
	outputFile.close();

	return 1;
}

void saveParams ( TransWorhpProblem* ph, string xmlparams ) {
	XMLParser p;
	XMLNode *n_ship     = p.Parse ( xmlparams )->GetFirstChild ( "SHIP" );
	XMLNode *n_model    = p.Parse ( xmlparams )->GetFirstChild ( "MODEL" );
	XMLNode *n_states   = p.Parse ( xmlparams )->GetFirstChild ( "STATES" );
	XMLNode *n_controls = p.Parse ( xmlparams )->GetFirstChild ( "CONTROLS" );
	XMLNode *n = p.Parse ( xmlparams )->GetFirstChild ( "PARAMS" );
	XMLNode *xml = n->GetFirstChild ( "DOUBLE" );
	
	int i=0;
	while ( xml ) {
		string idfy ( xml->GetAttribute ( "idfy" ) );
		if( idfy == "TRUE" ) {
			xml->SetText(ToString(ph->solver->X[ph->p_index ( i )]));
			i++;
		}
		xml = n->GetNextChild ( "DOUBLE" );
	}
	
	TextOutputType_e tot = ASCII;


		DebugStream d(tot);
		d << beginfile;
		
		if (n_ship) {
			d << *n_ship;
		} else {
			p.GetError(d);
		}
		if (n_model) {
			d << *n_model;
		} else {
			p.GetError(d);
		}
		if (n_states) {
			d << *n_states;
		} else {
			p.GetError(d);
		}
		if (n_controls) {
			d << *n_controls;
		} else {
			p.GetError(d);
		}
		if (n) {
			d << *n;
		} else {
			p.GetError(d);
		}
		
		
		d << endfile;

		ofstream file(xmlparams);
		file << d.GetString();
		
}


class SimOpt {
	public:
		
		vector<double> v;
		int nbrsec;
		
		vector<string> varName = {
			"Pos_x",
			"Pos_y",
			"Heading",
			"VelocLin_tw_u",
			"AccelLin_tw_u",
			"VelocLin_tw_v",
			"VelocRot_tw_r",
			"PU_0__RPMAct",
			"PU_0__PitchAct",
			"PU_0__EngineTorque",
			"EOT1",
			"RS_0__RE_Delta",
			"RS_0__RE_Cmd",
			"All_Force_X",
			"All_Force_Y",
			"All_Moment_N",
			"Hull_Force_X",
			"Hull_Force_Y",
			"Hull_Moment_N",
			"PU_0__Force_X",
			"PU_0__Force_Y",
			"PU_0__Moment_N",
			"PU_1__Force_X",
			"PU_1__Force_Y",
			"PU_1__Moment_N",
			"RS_0__Force_X",
			"RS_0__Force_Y",
			"RS_0__Moment_N",
			"RS_1__Force_X",
			"RS_1__Force_Y",
			"RS_1__Moment_N",
			"Wind_Force_X",
			"Wind_Force_Y",
			"Wind_Moment_N",
		};
		
		vector<double> Pos_x;
		vector<double> Pos_y;
		vector<double> Heading;
		vector<double> VelocLin_tw_u;
		vector<double> AccelLin_tw_u;
		vector<double> VelocLin_tw_v;
		vector<double> VelocRot_tw_r;
		vector<double> PU_0__RPMAct;
		vector<double> PU_0__PitchAct;
		vector<double> PU_0__EngineTorque;
		vector<double> EOT1;
		vector<double> RS_0__RE_Delta;
		vector<double> RS_0__RE_Cmd;
		vector<double> All_Force_X;
		vector<double> All_Force_Y;
		vector<double> All_Moment_N;
		vector<double> Hull_Force_X;
		vector<double> Hull_Force_Y;
		vector<double> Hull_Moment_N;
		vector<double> PU_0__Force_X;
		vector<double> PU_0__Force_Y;
		vector<double> PU_0__Moment_N;
		vector<double> PU_1__Force_X;
		vector<double> PU_1__Force_Y;
		vector<double> PU_1__Moment_N;
		vector<double> RS_0__Force_X;
		vector<double> RS_0__Force_Y;
		vector<double> RS_0__Moment_N;
		vector<double> RS_1__Force_X;
		vector<double> RS_1__Force_Y;
		vector<double> RS_1__Moment_N;
		vector<double> Wind_Force_X;
		vector<double> Wind_Force_Y;
		vector<double> Wind_Moment_N;
		
		
		SimOpt ( const char* filename ) {
			
			char line[3000];
			
			FILE* fd = fopen(filename,"r");
			int ch = 0;
			nbrsec = 0;
			while(!feof(fd))
			{
				ch = fgetc(fd);
				if(ch == '\n') nbrsec++;
			}
			nbrsec = nbrsec-4;
			fclose(fd);
			
			int varpos;
			int it;
			string buf;
			
			for(int i=0; i<varName.size(); i++) {
				v.clear();
				fd = fopen(filename,"r");
				for(int j=0; j<3; j++) fgets(line, 3000, fd);
				
				it = 0;
				stringstream ss1(line);
				while (ss1 >> buf) {
					if( buf == varName[i] ) {
						varpos = it;
						break;
					}
					it++;
				}
				
				fgets(line, 3000, fd);
				
				v.resize(nbrsec);
				
				for(int i=0; i<nbrsec; i++) {
					fgets(line, 3000, fd);
					it = 0;
					stringstream ss2(line); 
					while (ss2 >> buf) {
						if( it == varpos ) {
							v[i] = stod(buf,0);
							break;
						}
						it++;
					}
				}
				
				if(varName[i] == "Pos_x"         ) Pos_x          = v;
				if(varName[i] == "Pos_y"         ) Pos_y          = v;
				if(varName[i] == "Heading"       ) Heading        = v;
				if(varName[i] == "VelocLin_tw_u" ) VelocLin_tw_u  = v;
				if(varName[i] == "AccelLin_tw_u" ) AccelLin_tw_u  = v;
				if(varName[i] == "VelocLin_tw_v" ) VelocLin_tw_v  = v;
				if(varName[i] == "VelocRot_tw_r" ) VelocRot_tw_r  = v;
				if(varName[i] == "PU_0__RPMAct"  ) PU_0__RPMAct   = v;
				if(varName[i] == "PU_0__PitchAct") PU_0__PitchAct = v;
				if(varName[i] == "PU_0__EngineTorque") PU_0__EngineTorque = v;
				if(varName[i] == "EOT1"          ) EOT1           = v;
				if(varName[i] == "RS_0__RE_Delta") RS_0__RE_Delta = v;
				if(varName[i] == "RS_0__RE_Cmd"  ) RS_0__RE_Cmd   = v;
				if(varName[i] == "All_Force_X"   ) All_Force_X    = v;
				if(varName[i] == "All_Force_Y"   ) All_Force_Y    = v;
				if(varName[i] == "All_Moment_N"  ) All_Moment_N   = v;
				if(varName[i] == "Hull_Force_X"  ) Hull_Force_X   = v;
				if(varName[i] == "Hull_Force_Y"  ) Hull_Force_Y   = v;
				if(varName[i] == "Hull_Moment_N" ) Hull_Moment_N  = v;
				if(varName[i] == "PU_0__Force_X" ) PU_0__Force_X  = v;
				if(varName[i] == "PU_0__Force_Y" ) PU_0__Force_Y  = v;
				if(varName[i] == "PU_0__Moment_N") PU_0__Moment_N = v;
				if(varName[i] == "PU_1__Force_X" ) PU_1__Force_X  = v;
				if(varName[i] == "PU_1__Force_Y" ) PU_1__Force_Y  = v;
				if(varName[i] == "PU_1__Moment_N") PU_1__Moment_N = v;
				if(varName[i] == "RS_0__Force_X" ) RS_0__Force_X  = v;
				if(varName[i] == "RS_0__Force_Y" ) RS_0__Force_Y  = v;
				if(varName[i] == "RS_0__Moment_N") RS_0__Moment_N = v;
				if(varName[i] == "RS_1__Force_X" ) RS_1__Force_X  = v;
				if(varName[i] == "RS_1__Force_Y" ) RS_1__Force_Y  = v;
				if(varName[i] == "RS_1__Moment_N") RS_1__Moment_N = v;
				if(varName[i] == "Wind_Force_X"  ) Wind_Force_X   = v;
				if(varName[i] == "Wind_Force_Y"  ) Wind_Force_Y   = v;
				if(varName[i] == "Wind_Moment_N" ) Wind_Moment_N  = v;
				
				fclose(fd);
			}
		}
		
		double eval__(vector<double> &v, double t) {
			
			int idx0 = (int) floor(t);
			int idx1 = (int) ceil (t);
			double delta = t - idx0;
			
			if( t < 0 ) return v[0];
			if( t >= nbrsec-1) return v[nbrsec-1];
			if(idx0 == idx1) return v[idx0];
			return v[idx0] + (v[idx1] - v[idx0])*delta;
		}
		
		double tEnd() {
			return nbrsec;
		}
};


class OptTraj : public TransWorhpProblem {
	public:

		XMLNode* scenenode;

		vector <BaseModel*> ship_model;
		int NmodelS;
		int NmodelC;
		int NmodelP;

		
		vector<SimOpt*> SO;
		vector<double> tEnd;

		tw::TWfolder* folder_ptr;
		
		vector<double*> x0Array;
		vector<double*> y0Array;
		vector<double*> psiArray;
		vector<double*> uArray;
		vector<double*> upArray;
		vector<double*> vArray;
		vector<double*> rArray;
		vector<double*> deltaArray;
		vector<double*> nArray;

		OptTraj ( XMLNode *xmlmain,
				  const tw::TWdimension &TWdata,
				  const vector <BaseModel*> &ship_model,
				  TWfolder* folder
				)
			:
			TransWorhpProblem ( TWdata ) ,
			ship_model ( ship_model ) ,
			folder_ptr( folder ){

			if ( xmlmain ) scenenode = xmlmain->GetFirstChild ( "SCENE" );

			NmodelS = ship_model[0]->get_nbrS(  );
			NmodelC = ship_model[0]->get_nbrC(  );
			NmodelP = ship_model[0]->get_nbrIdParams(  );
			
			
			SO.resize(N);
			tEnd.resize(N);
			for(int i=0; i<N; i++) {
				SO[i] = new SimOpt(slog[i]);
				tEnd[i] = SO[i]->tEnd();
			}
			
			x0Array   .resize(N);
			y0Array   .resize(N);
			psiArray  .resize(N);
			uArray    .resize(N);
			vArray    .resize(N);
			rArray    .resize(N);
			deltaArray.resize(N);
			nArray    .resize(N);
			
			for(int i=0; i<N; i++) {
				x0Array   [i] = (double*) realloc(x0Array   [i], n_dis * sizeof(double));
				y0Array   [i] = (double*) realloc(y0Array   [i], n_dis * sizeof(double));
				psiArray  [i] = (double*) realloc(psiArray  [i], n_dis * sizeof(double));
				uArray    [i] = (double*) realloc(uArray    [i], n_dis * sizeof(double));
				vArray    [i] = (double*) realloc(vArray    [i], n_dis * sizeof(double));
				rArray    [i] = (double*) realloc(rArray    [i], n_dis * sizeof(double));
				deltaArray[i] = (double*) realloc(deltaArray[i], n_dis * sizeof(double));
				nArray    [i] = (double*) realloc(nArray    [i], n_dis * sizeof(double));
			}
		}

		void OpenWindows ( tw::Viewer *viewer ) {
			viewer->ThreeD ( "Ship tw::Viewer", scenenode, shipviewer3d );
			
			for(int i=0; i<N; i++) {
				for(int j=0; j<n_dis; j++) {
					x0Array   [i][j] =              SO[i] ->eval__(SO[i]->Pos_x         , tEnd[i]*j/(n_dis-1));
					y0Array   [i][j] =              SO[i] ->eval__(SO[i]->Pos_y         , tEnd[i]*j/(n_dis-1));
					psiArray  [i][j] = M_PI/180*    SO[i] ->eval__(SO[i]->Heading       , tEnd[i]*j/(n_dis-1));
					uArray    [i][j] =              SO[i] ->eval__(SO[i]->VelocLin_tw_u , tEnd[i]*j/(n_dis-1));
					vArray    [i][j] =              SO[i] ->eval__(SO[i]->VelocLin_tw_v , tEnd[i]*j/(n_dis-1));
					rArray    [i][j] = M_PI/180./60*SO[i] ->eval__(SO[i]->VelocRot_tw_r , tEnd[i]*j/(n_dis-1));
					deltaArray[i][j] = M_PI/180.   *SO[i] ->eval__(SO[i]->RS_0__RE_Delta, tEnd[i]*j/(n_dis-1));
					nArray    [i][j] =  1./60      *SO[i] ->eval__(SO[i]->PU_0__RPMAct  , tEnd[i]*j/(n_dis-1));
				}
				viewer->plots[NmodelS*i + 0]->AddCompareCurve2(solver->T.data(),x0Array   [i] , 1, n_dis);
				viewer->plots[NmodelS*i + 1]->AddCompareCurve2(solver->T.data(),y0Array   [i] , 1, n_dis);
				viewer->plots[NmodelS*i + 2]->AddCompareCurve2(solver->T.data(),psiArray  [i] , 1, n_dis);
				viewer->plots[NmodelS*i + 3]->AddCompareCurve2(solver->T.data(),uArray    [i] , 1, n_dis);
				viewer->plots[NmodelS*i + 4]->AddCompareCurve2(solver->T.data(),vArray    [i] , 1, n_dis);
				viewer->plots[NmodelS*i + 5]->AddCompareCurve2(solver->T.data(),rArray    [i] , 1, n_dis);
// 				viewer->plots[NmodelS*i + 6]->AddCompareCurve2(solver->T.data(),deltaArray[i] , 1, n_dis);
// 				viewer->plots[NmodelS*i + 6]->AddCompareCurve2(solver->T.data(),nArray    [i] , 1, n_dis);
			}
		}

		void selectWindows ( tw::Viewer *viewer ) {

			for ( int i=0; i<N*NmodelS; i++ ) {
				viewer->AddStateView ( i , "" );
			}
		}

		void p_init ( double *p ) {
		}

		double obj() {
			double objective = 0;
			
			for(int j=0; j<N; j++) {
				for(int i=0; i<n_dis; i++) objective += 1.e-5*( (x(i,NmodelS*j + 0) -              SO[j]->eval__(SO[j]->Pos_x        , tEnd[j]*i/(n_dis-1)))*(x(i,NmodelS*j + 0) -              SO[j]->eval__(SO[j]->Pos_x        , tEnd[j]*i/(n_dis-1))));
				for(int i=0; i<n_dis; i++) objective += 1.e-2*( (x(i,NmodelS*j + 1) -              SO[j]->eval__(SO[j]->Pos_y        , tEnd[j]*i/(n_dis-1)))*(x(i,NmodelS*j + 1) -              SO[j]->eval__(SO[j]->Pos_y        , tEnd[j]*i/(n_dis-1))));
				for(int i=0; i<n_dis; i++) objective += 1.e+3*( (x(i,NmodelS*j + 3) -              SO[j]->eval__(SO[j]->VelocLin_tw_u, tEnd[j]*i/(n_dis-1)))*(x(i,NmodelS*j + 3) -              SO[j]->eval__(SO[j]->VelocLin_tw_u, tEnd[j]*i/(n_dis-1))));
				for(int i=0; i<n_dis; i++) objective += 1.e+3*( (x(i,NmodelS*j + 4) -              SO[j]->eval__(SO[j]->VelocLin_tw_v, tEnd[j]*i/(n_dis-1)))*(x(i,NmodelS*j + 4) -              SO[j]->eval__(SO[j]->VelocLin_tw_v, tEnd[j]*i/(n_dis-1))));
				for(int i=0; i<n_dis; i++) objective += 1.e+7*( (x(i,NmodelS*j + 5) - M_PI/180./60*SO[j]->eval__(SO[j]->VelocRot_tw_r, tEnd[j]*i/(n_dis-1)))*(x(i,NmodelS*j + 5) - M_PI/180./60*SO[j]->eval__(SO[j]->VelocRot_tw_r, tEnd[j]*i/(n_dis-1))));
// 				for(int i=0; i<n_dis; i++) objective += 1.e+1*( (x(i,NmodelS*j + 6) - M_PI/180.   *SO[j]->eval__(SO[j]->RS_0__RE_Delta,tEnd[j]*i/(n_dis-1)))*(x(i,NmodelS*j + 6) - M_PI/180.   *SO[j]->eval__(SO[j]->RS_0__RE_Delta,tEnd[j]*i/(n_dis-1))));
// 				for(int i=0; i<n_dis; i++) objective += 1.e+2*( (x(i,NmodelS*j + 7) - 1./60       *SO[j]->eval__(SO[j]->PU_0__RPMAct          ,tEnd[j]*i/(n_dis-1)))*(x(i,NmodelS*j + 7) - 1./60     *SO[j]->eval__(SO[j]->PU_0__RPMAct          ,tEnd[j]*i/(n_dis-1))));
			}
			objective *= 1.e-4;
			return objective;
		}

		bool obj_structure ( tw::DiffStructure &s ) {
			return false;
		}

		bool obj_diff ( tw::DiffStructure &s ) {
			return false;
		}

		void ode ( double *dx, double t, const double *x, const double *u, const double *p ) {

			double ctrl  [ NmodelC ];
			double state [ NmodelS ];
			double dstate[ NmodelS ];
			
			double Xrest;
			double Yrest;
			double Nrest;

			for(int i=0; i<N; i++) {
				
				for ( int j=0; j<NmodelS; j++ ) state[j] = x[NmodelS*i + j];
				
				double fac = 1;
				Xrest = SO[i]->eval__(SO[i]->All_Force_X ,t*tEnd[i]) - SO[i]->eval__(SO[i]->Hull_Force_X ,t*tEnd[i]) - SO[i]->eval__(SO[i]->PU_0__Force_X ,t*tEnd[i]) - SO[i]->eval__(SO[i]->PU_1__Force_X ,t*tEnd[i]) - SO[i]->eval__(SO[i]->RS_0__Force_X ,t*tEnd[i]) - SO[i]->eval__(SO[i]->RS_1__Force_X ,t*tEnd[i]) - SO[i]->eval__(SO[i]->Wind_Force_X ,t*tEnd[i]);
				Yrest = SO[i]->eval__(SO[i]->All_Force_Y ,t*tEnd[i]) - SO[i]->eval__(SO[i]->Hull_Force_Y ,t*tEnd[i]) - SO[i]->eval__(SO[i]->PU_0__Force_Y ,t*tEnd[i]) - SO[i]->eval__(SO[i]->PU_1__Force_Y ,t*tEnd[i]) - SO[i]->eval__(SO[i]->RS_0__Force_Y ,t*tEnd[i]) - SO[i]->eval__(SO[i]->RS_1__Force_Y ,t*tEnd[i]) - SO[i]->eval__(SO[i]->Wind_Force_Y ,t*tEnd[i]);
				Nrest = SO[i]->eval__(SO[i]->All_Moment_N,t*tEnd[i]) - SO[i]->eval__(SO[i]->Hull_Moment_N,t*tEnd[i]) - SO[i]->eval__(SO[i]->PU_0__Moment_N,t*tEnd[i]) - SO[i]->eval__(SO[i]->PU_1__Moment_N,t*tEnd[i]) - SO[i]->eval__(SO[i]->RS_0__Moment_N,t*tEnd[i]) - SO[i]->eval__(SO[i]->RS_1__Moment_N,t*tEnd[i]) - SO[i]->eval__(SO[i]->Wind_Moment_N,t*tEnd[i]);
				
				ctrl[0] = 1000* Xrest;
				ctrl[1] = 1000* Yrest;
				ctrl[2] = 1000* Nrest;
				ctrl[3] = SO[i]->eval__(SO[i]->RS_0__RE_Delta ,t*tEnd[i])*M_PI/180;
				ctrl[4] = SO[i]->eval__(SO[i]->PU_0__RPMAct ,t*tEnd[i])/60;
// 				ctrl[3] = SO[i]->eval__(SO[i]->RS_0__RE_Cmd ,t*tEnd[i])*M_PI/180;
// 				ctrl[3] = SO[i]->eval__(SO[i]->EOT1         ,t*tEnd[i])/1000;
				ctrl[5] = SO[i]->eval__(SO[i]->PU_0__PitchAct ,t*tEnd[i])/100;
			
				
				double q[9] = {0,0,0,NAN,0,NAN,NAN,NAN,1000.};
				
				ship_model[i]->ode ( dstate, state, ctrl, p, q );

				for ( int j=0; j<NmodelS; j++ ) dx[NmodelS*i + j] = dstate[j];

				for ( int j=0; j<NmodelS; j++ ) dx[NmodelS*i + j] *= tEnd[i];
			}

		}

		bool ode_structure ( tw::DiffStructure &s ) {
			return false;
		}

		void p_boundary ( double *p_low, double *p_upp ) {
			ship_model[0]->get_param_bounds ( p_low, p_upp );
		}

		void x_boundary ( double *x_low, double *x_upp ) {

			double low[NmodelS+NmodelC];
			double upp[NmodelS+NmodelC];
			ship_model[0]->get_bounds ( low, upp );

			for(int i=0; i<N; i++) {
				for ( int j=0; j<NmodelS; j++ ) {
					x_low[NmodelS*i + j] = low[j];
					x_upp[NmodelS*i + j] = upp[j];
				}
			}
		}

		void var_boundary ( double *x_low, double *x_upp ) {

// 			for ( int i=0; i<NmodelS; i++ ) {
// 				x_low[x_index ( 0,i ) ] = x_upp[x_index ( 0,i ) ] = ship_model[0]->initialS[i];
// 			}
			
// 			x_low[x_index ( 0,0 ) ] = x_upp[x_index ( 0,0 ) ] = x0->eval(0);
// 			x_low[x_index ( 0,1 ) ] = x_upp[x_index ( 0,1 ) ] = y0->eval(0);
// 			x_low[x_index ( 0,2 ) ] = x_upp[x_index ( 0,2 ) ] = M_PI/180.*psi->eval(0);
// 			x_low[x_index ( 0,3 ) ] = x_upp[x_index ( 0,3 ) ] = u->eval(0);
// 			x_low[x_index ( 0,4 ) ] = x_upp[x_index ( 0,4 ) ] = v->eval(0);
// 			x_low[x_index ( 0,5 ) ] = x_upp[x_index ( 0,5 ) ] = M_PI/180./60*r->eval(0);
		}


		void terminate() {
// 			double ctrl  [ NmodelC ];
// 			double state [ NmodelS ];
// 			double dstate[ NmodelS ];
// 			
// 			for(int i=0; i<NmodelP; i++) cout << p(i) << endl;
// 			double m = (p(1)*1e-0)/(p(0)*1e-7);
// 			double Xup = m - 1/(p(0)*1e-7);
// 			double xG = (p(2)*1e+1)*(m - Xup)/m;
		}

};





int main ( int argv, char* argc[] ) {

	tw::TWparameter twparameter ( "transworhp.xml" );
	map<string,string> args = twparameter.Arguments ( argv,argc );

	XMLNode *xml_shipviewer = tw::TWparameter::ReadParams ( "shipviewer.xml" );

	Viewer *viewer = 0;

	if ( twparameter.PLOT ) viewer = new tw::Viewer ( &twparameter );

	vector <BaseModel*> ship_model ( N );
	for(int i=0; i<N; i++) ship_model[i] = new Model1v4 ( params );

	int NmodelS;
	int NmodelC;
	int NmodelP;
	NmodelS = ship_model[0] -> get_nbrS       (  );
	NmodelC = ship_model[0] -> get_nbrC       (  );
    NmodelP = ship_model[0] -> get_nbrIdParams(  );
	
	tw::TWdimension TWdim;
	TWdim.ID      = "Optimal Trajectories";
	TWdim.n_dis   = twparameter.NDIS;
	TWdim.n_ode   = N*NmodelS;
	TWdim.n_ctrl  = 0;
	TWdim.n_param = NmodelP;
	TWdim.n_rand  = 0;
	

	{
		tw::TWfolder folder ( &twparameter,0 );
		
		OptTraj ph ( xml_shipviewer, TWdim, ship_model, &folder );

		ph0 = &ph;

		ph.setSolver ( &twparameter );
		folder.Add ( &ph );
		folder.Init();
		
		vector<SimOpt*> SO;
		vector<double> tEnd;
		SO.resize(N);
		tEnd.resize(N);
		for(int i=0; i<N; i++) {
			SO[i] = new SimOpt(slog[i]);
			tEnd[i] = SO[i]->tEnd();
		}

		viewer->SetFloatTime ( tEnd[0] /60 , 0 );

		for(int i=0; i<N; i++) {
			ph.solver->X[ph.x_index ( 0,NmodelS*i + 0 )] =               SO[i] ->eval__(SO[i]->Pos_x         ,0);
			ph.solver->X[ph.x_index ( 0,NmodelS*i + 1 )] =               SO[i] ->eval__(SO[i]->Pos_y         ,0);
			ph.solver->X[ph.x_index ( 0,NmodelS*i + 2 )] = M_PI/180.*    SO[i] ->eval__(SO[i]->Heading       ,0);
			ph.solver->X[ph.x_index ( 0,NmodelS*i + 3 )] =               SO[i] ->eval__(SO[i]->VelocLin_tw_u ,0);
			ph.solver->X[ph.x_index ( 0,NmodelS*i + 4 )] =               SO[i] ->eval__(SO[i]->VelocLin_tw_v ,0);
			ph.solver->X[ph.x_index ( 0,NmodelS*i + 5 )] = M_PI/180./60* SO[i] ->eval__(SO[i]->VelocRot_tw_r ,0);
// 			ph.solver->X[ph.x_index ( 0,NmodelS*i + 6 )] = M_PI/180.   * SO[i] ->eval__(SO[i]->RS_0__RE_Delta,0);
// 			ph.solver->X[ph.x_index ( 0,NmodelS*i + 7 )] = 1./60       * SO[i] ->eval__(SO[i]->PU_0__RPMAct  ,0);
		}

		int j=0;
		for(int i=0; i<ship_model[0]->nbrP; i++) {
			if(ship_model[0]->paramIdfy[i]) { 
				ph.solver->X[ph.p_index ( j )] = *(ship_model[0]->paramAddr[i])/ship_model[0]->paramScales[i]; 
				j++; 
			}
		}
		
		int steps = ph.solver->Integrate ( twparameter.butchertableau );
		cout << "Integrationsschritte: " << steps << endl;

		folder.Init ( viewer );
		folder.Loop(3000);
		
// 		CopyFile(params, "ParamXML/s1v1_m1v4_p2_backup.xml");
// 		saveParams ( &ph, params );

		ph.solver->ToMATLAB ( "opt_traj.m" );
	}

	if ( viewer ) viewer->CloseAll();
	
	delete viewer;
	delete xml_shipviewer;
	delete ship_model[0];
	

	return 0;
}



