#include "pathcost.h"

namespace dsl
{
	
PathCost2d ::PathCost2d(){};
	
PathCost2d ::~PathCost2d(){};

double PathCost2d ::Heur ( const TypedCell2d& a, const TypedCell2d& b ) const
{	
	/*
	double x = fabs(a.c[0]-b.c[0]);
	double y = fabs(a.c[1]-b.c[1]);
	x=fabs(x);
	y=fabs(y);
	if(x>y) return sqrt(2.*y*y)+x-y;
	else return sqrt(2.*x*x)+y-x;*/
	return sqrt( (a.c[0]-b.c[0])*(a.c[0]-b.c[0]) + (a.c[1]-b.c[1])*(a.c[1]-b.c[1]) );
};

double PathCost2d ::Real ( const TypedCell2d& a, const TypedCell2d& b ) const
{
	std::cout<<"Err PathCost2d::Real not implemented"<<std::endl;
	assert(false);
	return -1.;
};

double PathCost2d ::EdgeCost ( const TypedCell2d& a, const TypedCell2d& b ) const
{
	return sqrt( (a.c[0]-b.c[0])*(a.c[0]-b.c[0]) + (a.c[1]-b.c[1])*(a.c[1]-b.c[1]) );
};




	
PathCost3d ::PathCost3d()
	:head_scale(500./(2.*M_PI))
{};

PathCost3d ::~PathCost3d(){};

double PathCost3d ::Heur ( const TypedCell3d& a, const TypedCell3d& b ) const
{	
	/*
	double radius=400.0;
	double normal_down[2] = {-1,0};
	double diff[2] = {a.c[0]-b.c[0],a.c[1]-b.c[1]};
	double frac = -1*diff[0]/sqrt(pow(diff[0],2)+pow(diff[1],2));
	//cout << frac << endl;
	double local_a2 = diff[1]<0? acos(frac):2*M_PI-acos(frac);
	if (!isnormal(local_a2)){
		local_a2 =0;
	}
	double q0[3] = {a.c[0], a.c[1], local_a2};//x,y,psi
	//double q0[3] = {a.c[0], a.c[1], a.c[2]};//x,y,psi
	//std::cout << a.c[2] << std::endl;
	double q1[3] = {b.c[0], b.c[1], b.c[2]};
	DubinsPath path;
	dubins_init( q0, q1, radius, &path);
	//std::cout << "diff: " << diff[0] << " " << diff[1] << endl; 
	//std::cout <<"a:    " << a.c[0]<< " " << a.c[1]<< " " << a.c[2] << std::endl;
	//cout << 180*local_a2/M_PI << endl;
	return  dubins_path_length( &path );
	*/
	//std::cout << b.c[0]<< " " << b.c[1]<< " " << b.c[2] << std::endl;
	/*bool once = true;
	if (once){
		std::cout << b.c[0]<< " " << b.c[1]<< " " << b.c[2] << std::endl;
		once = false;
	}*/
	return sqrt( (a.c[0]-b.c[0])*(a.c[0]-b.c[0]) + (a.c[1]-b.c[1])*(a.c[1]-b.c[1]));
};

double PathCost3d ::Real ( const TypedCell3d& a, const TypedCell3d& b ) const
{
	std::cout<<"Err PathCost3d::Real not implemented"<<std::endl;
	assert(false);
	return -1.;
};

double PathCost3d ::EdgeCost ( const TypedCell3d& a, const TypedCell3d& b, double dist ) const
{
	return dist;
};

// neu 4D extension

PathCost4d ::PathCost4d()
	:head_scale(500./(2.*M_PI))
{};

PathCost4d ::~PathCost4d(){};

double PathCost4d ::Heur ( const TypedCell4d& a, const TypedCell4d& b ) const
{	
	return sqrt( (a.c[0]-b.c[0])*(a.c[0]-b.c[0]) + (a.c[1]-b.c[1])*(a.c[1]-b.c[1]))*1.1;
};

double PathCost4d ::Real ( const TypedCell4d& a, const TypedCell4d& b ) const
{
	std::cout<<"Err PathCost4d::Real not implemented"<<std::endl;
	assert(false);
	return -1.;
};

double PathCost4d ::EdgeCost ( const TypedCell4d& a, const TypedCell4d& b, double dist ) const
{
	return dist;
};


}