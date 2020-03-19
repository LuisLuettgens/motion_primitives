#include <iostream>
#include <fstream>
#include <cmath>
#include "clipper.hpp"
 
using namespace std;
using namespace ClipperLib;
 
// void GetPolygonSub(Paths& p)
// {
// 	p[0].push_back(IntPoint(1,1));
// 	p[0].push_back(IntPoint(1,3));
// 	p[0].push_back(IntPoint(3,3));
// 	p[0].push_back(IntPoint(3,1));
// }
//  
// void GetPolygonClp(Paths& p)
// {
// 	p[0].push_back(IntPoint(2,2));
// 	p[0].push_back(IntPoint(2,4));
// 	p[0].push_back(IntPoint(4,4));
// 	p[0].push_back(IntPoint(4,2));
// }

double scaling_fac = 10000000000;

void loadPolygonsFromFile(Paths& p, string filename){
	ifstream f;
	string str;
	f.open ( filename.c_str(), ios::in );
	if ( !f ) {
		throw runtime_error ( "Error Reading File in loadPolygonsFromFile; Wrong Filename/Path?" );
	}
	
	string xcoord;
	string ycoord;
	
	getline ( f,str );
	while ( !f.eof() ) {
		
		if(str.size()>0){
			Path poly;
			while( str.size()>0 ){
			
				unsigned int i=0;
				while ( i<str.size() && str.at ( i ) ==' ' ) {
					i++;
				}
				xcoord  = "";
				ycoord  = "";
				while ( i<str.size() && str.at ( i ) !=' ' ) {
					xcoord += str.at ( i );
					i++;
				}
				while ( i<str.size() && str.at ( i ) ==' ' ) {
					i++;
				}
				while ( i<str.size() && str.at ( i ) !=' ' ) {
					ycoord += str.at ( i );
					i++;
				}
				poly.push_back ( IntPoint( (cInt) (atof(ycoord.c_str())*scaling_fac),
										   (cInt) (atof(xcoord.c_str())*scaling_fac) ) );
				getline ( f,str );
			}
			p.push_back(poly);
			getline ( f,str );
		}
	}
	f.close();
}

void savePolygonsToFile(Paths& p, string filename){
	ofstream f;
	string str;
	f.open ( filename.c_str(), ios::out );
	if ( !f ) {
		throw runtime_error ( "Error Reading File in savePolygonsToFile; Wrong Filename/Path?" );
	}
	f.precision(14);
	
	for (unsigned int i=0; i<p.size(); i++) {
		for (unsigned int j=0; j<p[i].size(); j++) {
			double x=p[i][j].X/scaling_fac;
			double y=p[i][j].Y/scaling_fac;
			f<<y<<" "<<x<<endl; // x and y are interchanged!?
		}
		f<<endl;
	}
	
}
 
int main( int argc, char *argv[] )
{
	
	if ( argc!=4 ) {
		cout<<"Error: Wrong number of arguments, must be 3!"<<endl;
		exit(1);
	}

	string filename_sub = argv[1];
	string filename_clp = argv[2];
	string filename_out = argv[3];
	
	
  //set up the subject and clip polygons ...
  Paths sub;
//   GetPolygonSub(sub);
  loadPolygonsFromFile(sub, filename_sub);
   
  Paths clp;
//   GetPolygonClp(clp);
  loadPolygonsFromFile(clp, filename_clp);
   
  //get the intersection of the subject and clip polygons ...
  Clipper clpr;
  clpr.AddPaths(sub, ptSubject, true);
  clpr.AddPaths(clp, ptClip, true);
  Paths solution;
  clpr.Execute(ctIntersection, solution, pftEvenOdd, pftEvenOdd);
  
  savePolygonsToFile(solution, filename_out);
  
//   for (unsigned int i=0; i<solution.size(); i++) {
// 	  cout<<"Poly "<<i<<":"<<endl;
// 	  for (unsigned int j=0; j<solution[i].size(); j++) {
// 		  cout<<"x="<<solution[i][j].X<<", y="<<solution[i][j].Y<<endl;
// 	  }
// 	  cout<<endl;
//   }
  
}