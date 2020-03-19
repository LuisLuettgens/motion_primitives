#ifndef PATHCOST_H
#define PATHCOST_H

#include "utils/functions.h"
#include "dubins.h"
#include "cost.h"
#include "data_types.h"
#include <fstream>
#include <math.h> 
namespace dsl
{


class PathCost2d : public Cost<TypedCell2d>{
public:
	PathCost2d();
	~PathCost2d();
	double Heur ( const TypedCell2d& a, const TypedCell2d& b ) const;
	double Real ( const TypedCell2d& a, const TypedCell2d& b ) const;
	double EdgeCost ( const TypedCell2d& a, const TypedCell2d& b ) const;
};


class PathCost3d : public Cost<TypedCell3d>{
public:
	PathCost3d();
	~PathCost3d();
	double Heur ( const TypedCell3d& a, const TypedCell3d& b ) const;
	double Real ( const TypedCell3d& a, const TypedCell3d& b ) const;
	double EdgeCost ( const TypedCell3d& a, const TypedCell3d& b, double dist ) const;

	double head_scale;
};
// neu 4d extension
class PathCost4d : public Cost<TypedCell4d>{
public:
	PathCost4d();
	~PathCost4d();
	double Heur ( const TypedCell4d& a, const TypedCell4d& b ) const;
	double Real ( const TypedCell4d& a, const TypedCell4d& b ) const;
	double EdgeCost ( const TypedCell4d& a, const TypedCell4d& b, double dist ) const;

	double head_scale;
};

}

#endif