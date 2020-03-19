#ifndef DSL_PATHGRID_H
#define DSL_PATHGRID_H

#include "data_types.h"
#include "grid.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>

namespace dsl
{
	
struct PathGrid2d : public Grid<Eigen::Vector2d, VertData2d>{
public:
    PathGrid2d ( const Eigen::Vector2d& xlb, const Eigen::Vector2d& xub, const Vectorni& gs );
    TypedCell2d* Get ( const Eigen::Vector2d& x, bool checkValid = true ) const;
    TypedCell2d* Get ( int id ) const;

};

struct PathGrid3d : public Grid<Eigen::Vector3d, VertData3d>{
public:
PathGrid3d ( const Eigen::Vector3d& xlb, const Eigen::Vector3d& xub, const Vectorni& gs );
TypedCell3d* Get ( const Eigen::Vector3d& x, bool checkValid = true ) const;
TypedCell3d* Get ( int id ) const;

};

// neu 4D extension

struct PathGrid4d : public Grid<Eigen::Vector4d, VertData4d>{
public:
PathGrid4d ( const Eigen::Vector4d& xlb, const Eigen::Vector4d& xub, const Vectorni& gs );
TypedCell4d* Get ( const Eigen::Vector4d& x, bool checkValid = true ) const;
TypedCell4d* Get ( int id ) const;

};

}

#endif
