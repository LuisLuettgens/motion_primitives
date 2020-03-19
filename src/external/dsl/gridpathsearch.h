// This file is part of libdsl, a library for heuristic graph search
//
// Copyright (C) 2004 Marin Kobilarov <marin@jhu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DSL_GRIDPATHSEARCH_H
#define DSL_GRIDPATHSEARCH_H

#include "search.h"
#include "pathcost.h"
#include "pathgrid.h"
#include "data_types.h"
#include "prim.h"
#include "utils.h"
#include "map.h"
#include <omp.h>
#include <Eigen/Geometry>
#include <vector>
#include <iostream>
#include "dirent.h"
#include <fstream>
#include <map>

#include "utils/functions.h"
#include "geometric_constr.h"

/**
 *  Grid based D*-Lite, extends graph-based D* Lite.
 *  The class maps grid cell costs and the transition costs b/n two cells
 *  into graph edge costs of the base Search class
 *
 *  Supports n-dimensional grid with abstract "connectivity" interface
 *
 *  For example in 2D, the cell transition costs can be encoded as the euclidean
 *distance
 *  b/n the centers of the cells (i.e. each cell has 8 neighbors:
 *  the neighbors at N,S,W,E have transition cost of 1, and the
 *  neighbors at NE, SE, NW, SW have transition costs of sqrt(2)
 *  These transition costs are added to the maximum of the values of
 *  two neighboring cells to compute the cost of their connecting edge.
 *
 *
 *  The planner is used as follows:
 *
 *  1) GridPathSearch2d(grid, connectivity, cost)
 *  2) SetStart(Vector2d(x, y))
 *  3) SetGoal(Vector2d(x, y))
 *  4) Plan(path)
 *  5) follow path until map changes are observed
 *  6) for each change: SetCost(Vector2d(x, y), cost)
 *  7) SetStart(Vector2d(x,y)) -- change the start to the current robot position
 *  8) goto 4
 *
 *
 *  Author: Marin Kobilarov
 */

namespace dsl
{


const double dt = 25.; 
const bool VERBOSE = false;

using std::vector;

class SimplePathSearch{
public:
    SimplePathSearch(GeometricConstr& geom_constr, Eigen::Vector2d d_xy_ned_grid,
                    Eigen::Vector2d start_lla, Eigen::Vector2d goal_lla);

    vector< Vec2d > solution;

};


class GridPathSearch2d : public Search< TypedCell2d, EdData2d >
{
public:
	
    GridPathSearch2d(TypedGraph2d& graph_, TypedMap& map_, PathCost2d& cost_, PathGrid2d& grid_,
                     double lat0_, double lon0_, double hgt0_, double d_occ_map_);

    virtual ~GridPathSearch2d();

    /**
     * Set start location in the map
     * @param x euclidean point vector
     * @return true on success
     */
    bool SetStart ( const Eigen::Vector2d& x, VertData2d data_ );

    /**
     * Set goal location in the map
     * @param x euclidean point vector
     * @return true on success
     */
    bool SetGoal ( const Eigen::Vector2d& x, VertData2d data_ );

    /**
     * Expand successors or predecessors of a given vertex. This is mostly used
     * internally
     * during search. In some cases can be "pre-called" externally in advance of
     * searching in this area of
     * the grid, for efficiency purposes.
     * @param from the given vertex
     * @param fwd if true then expand forward in time, i.e. successors, otherwise
     * expand predecessors
     * @return true on success
     */
    bool Expand ( TypedVertex2d& from, bool fwd = true );
	bool Flow ( TypedCell2d* from_cell, TypedCell2d*& to_cell, EdData2d* EdData2d, int motion_type, bool fwd );

    const PathGrid2d& grid; ///< the grid
    
	const TypedMap& map; ///<occupancy map from ppm file

    TypedVertex2d** vertexMap; ///< vertex grid array
	
	std::vector< Eigen::Vector2d > visited_pos;
    
    double lat0;
	double lon0;
    double hgt0;
    
    double d_occ_map;
    
    bool save_sol_to_file;
};



class PrimPathSearch{
public:
    PrimPathSearch(GeometricConstr& geom_constr, Eigen::Vector3d d_nepsi_grid,
                    Eigen::Vector3d start_lla, Eigen::Vector3d goal_lla);
    void create_init_guess(	vector< dsl::TypedEdge3d* > sol);
    void create_init_guess_with_holes(	vector< dsl::TypedEdge3d* > sol,vector<vector <double> > intersection_positions);
    void create_init_guess_without_reparsing(	vector< dsl::TypedEdge3d* > sol,vector<vector <double> > intersection_positions, std::vector< prim > prims);
    
    vector< vector<double> > solution;
        vector<double> start_pos;
        vector<double> final_pos;
        vector<vector <double> > intersection_positions_w_gaps;
        vector<vector <double> > intersection_positions_wo_gaps;
};

class ArcPathSearch{
public:
    ArcPathSearch(GeometricConstr& geom_constr, Eigen::Vector3d d_nepsi_grid,
                    Eigen::Vector3d start_lla, Eigen::Vector3d goal_lla);

    vector< vector<double> > solution;

};


class PrimGridPathSearch3d : public Search< TypedCell3d, EdData3d >
{
public:
	
    PrimGridPathSearch3d(TypedGraph3d& graph_, TypedMap& map_, PathCost3d& cost_,
                     PathGrid3d& grid_, double lat0_, double lon0_, double hgt0_, double d_xy_,
                     Eigen::Vector3d d_nepsi_grid_, double d_occ_map_, int n_trims_);

    virtual ~PrimGridPathSearch3d();

    /**
     * Set start location in the map
     * @param x euclidean point vector
     * @return true on success
     */
    bool SetStart ( const Eigen::Vector3d& x, VertData3d data_ );

    /**
     * Set goal location in the map
     * @param x euclidean point vector
     * @return true on success
     */
    bool SetGoal ( const Eigen::Vector3d& x, VertData3d data_ );

    /**
     * Expand successors or predecessors of a given vertex. This is mostly used
     * internally
     * during search. In some cases can be "pre-called" externally in advance of
     * searching in this area of
     * the grid, for efficiency purposes.
     * @param from the given vertex
     * @param fwd if true then expand forward in time, i.e. successors, otherwise
     * expand predecessors
     * @return true on success
     */
    bool Expand ( TypedVertex3d& from, bool fwd = true );
   
    bool Flow ( TypedCell3d* from_cell, TypedCell3d*& to_cell, EdData3d* EdData3d, int motion_type,
                bool fwd);
    //template<bool save_path>
   
    template <typename T>
    void zero(T** buf,std::size_t size);

    double CropPsi ( double psi );

    template<bool save_path_bool, bool crop_psi_bool>
    bool ExecutePrim( const Eigen::Vector3d& v0, Eigen::Vector3d& vf, double& d, const int prim_nr, const double dt, const bool fwd);
    
    template<bool save_path_bool, bool crop_psi_bool>
    bool ExecuteManeuver( const Eigen::Vector3d& v0, Eigen::Vector3d& vf, double& d, const int man_nr, const bool fwd);
    
    template<bool save_path_bool, bool crop_psi_bool>
    bool ExecuteTrim( const Eigen::Vector3d& v0, Eigen::Vector3d& vf, double& d, const int trim_nr, const double dt, const bool fwd);

    bool SetPrimitives ( std::string prim_path );
    
    const PathGrid3d& grid; ///< the grid
        
	const TypedMap& map; ///<occupancy map from ppm file

    TypedVertex3d** vertexMap; ///< vertex grid array
	
    std::vector< Eigen::Vector3d > visited_pos;
    std::vector< std::vector<double> > saved_path;
    std::vector< std::vector<double> > saved_ship_state;
    
    double lat0;
	double lon0;
    double hgt0;
    
    std::vector< prim > prims;
    int n_mans;
    int n_prims;
    int n_trims;
    int nmbr_motion_primitives;
    double d_xy;
    std::vector<std::vector<int >> prim_links_fwd;
    std::vector<std::vector<int >> prim_links_bwd;
    std::vector<int> result_trim_fwd;
    std::vector<int> result_trim_bwd;
    std::vector<double> dt_vec;

    bool luis_once = true;

    //arc-motion parameters
    //double radius_arc=518.3283638044612;//Thats the radius of motion trim 5 resp. 7
    //double angle_arc=10.*M_PI/180.;
    double straight_mov_fac=4.;

    Eigen::Vector3d d_nepsi_grid;

    TypedCell3d* no_heading_goal_cell;

    double d_occ_map;                                         

    
    bool save_sol_to_file;
};

class GridPathSearch3d : public Search< TypedCell3d, EdData3d >
{
public:
	
    GridPathSearch3d(TypedGraph3d& graph_, TypedMap& map_, PathCost3d& cost_,
                     PathGrid3d& grid_, double lat0_, double lon0_, double hgt0_,
                     Eigen::Vector3d d_nepsi_grid_, double d_occ_map_);

    virtual ~GridPathSearch3d();

    /**
     * Set start location in the map
     * @param x euclidean point vector
     * @return true on success
     */
    bool SetStart ( const Eigen::Vector3d& x, VertData3d data_ );

    /**
     * Set goal location in the map
     * @param x euclidean point vector
     * @return true on success
     */
    bool SetGoal ( const Eigen::Vector3d& x, VertData3d data_ );

    /**
     * Expand successors or predecessors of a given vertex. This is mostly used
     * internally
     * during search. In some cases can be "pre-called" externally in advance of
     * searching in this area of
     * the grid, for efficiency purposes.
     * @param from the given vertex
     * @param fwd if true then expand forward in time, i.e. successors, otherwise
     * expand predecessors
     * @return true on success
     */

    bool Expand ( TypedVertex3d& from, bool fwd = true );
    bool Flow ( TypedCell3d* from_cell, TypedCell3d*& to_cell, EdData3d* EdData3d, int motion_type,
                bool fwd);
    template<bool save_path>
    bool ExecuteMotion(Eigen::Vector3d v0_nepsi, Eigen::Vector3d& vf_nepsi,
                       int motion_type, double& dist, bool fwd);
    double CropPsi ( double psi );

    const PathGrid3d& grid; ///< the grid
    
	const TypedMap& map; ///<occupancy map from ppm file

    TypedVertex3d** vertexMap; ///< vertex grid array
	
    std::vector< Eigen::Vector3d > visited_pos;
    std::vector< std::vector<double> > saved_path;

    double lat0;
	double lon0;
    double hgt0;
    
    //arc-motion parameters
    double radius_arc=518.3283638044612;//Thats the radius of motion trim 5 resp. 7
    double angle_arc=10.*M_PI/180.;
    double straight_mov_fac=4.;

    Eigen::Vector3d d_nepsi_grid;

    TypedCell3d* no_heading_goal_cell;

    double d_occ_map;

    //values are: u, v, r, delta, n
    vector< vector<double> > motion_values = {  // straight (trim 6)
                                                {4.2480797, 0., 0., 0., 0.8333333},
                                                //short-left (trim 5)
                                                {2.0639419, 0.5631852, -0.0041275, 0.5235988, 0.8333333},
                                                //long-left (trim 5)
                                                {2.0639419, 0.5631852, -0.0041275, 0.5235988, 0.8333333},
                                                //short-right (trim 7)
                                                {2.0639419, -0.5631852, 0.0041275, -0.5235988, 0.8333333},
                                                //long-right (trim 7)
                                                {2.0639419, -0.5631852, 0.0041275, -0.5235988, 0.8333333} };

    bool save_sol_to_file;
};

}

#endif
