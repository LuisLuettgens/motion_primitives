#ifndef DSL_UTILS_H
#define DSL_UTILS_H

#include <sys/time.h>
#include <Eigen/Dense>
#include "map.h"
#include <vector>

namespace dsl
{

Map<bool, 2> load ( const char* filename, const Eigen::Vector2d &cs );

void save ( const dsl::Map<bool, 2> &map, const char* filename, const std::vector<Eigen::Vector2d> *path = 0 );

void save_map ( const char* map, int width, int height, const char* filename );

char* load_map ( int &width, int &height, const char* filename );

void timer_start ( struct timeval* time );

long timer_us ( struct timeval* time );

double fangle ( double a );

void se2_q2g ( Eigen::Matrix3d& m, const Eigen::Vector3d& q );

void se2_g2q ( Eigen::Vector3d& q, const Eigen::Matrix3d& m );

void se2_inv ( Eigen::Matrix3d& mi, const Eigen::Matrix3d& m );

void se2_exp ( Eigen::Matrix3d& m, const Eigen::Vector3d& v, double tol = 1e-16 );

void replaceExt ( std::string& s, const std::string& newExt );

/**
 * Function that returns a sign for object of any class as long as the operator
 * - and operator < are defined
 * for that class
 * @param val the object whose sign we need to check
 * @return sign of the object
 */
template < typename T >
int sgn ( T val )
{
    return ( T ( 0 ) < val ) - ( val < T ( 0 ) );
}
}

#endif
