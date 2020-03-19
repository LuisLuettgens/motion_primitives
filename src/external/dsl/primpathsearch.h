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


const double dt = 25.; // execution time of trim primitives
const bool VERBOSE = false; // variable that can be used for output information

using std::vector;
class NEW_PrimPathSearch{
public:
    NEW_PrimPathSearch(GeometricConstr& geom_constr, Eigen::Vector4d d_nepsi_grid,
                    Eigen::Vector3d start_lla, Eigen::Vector3d goal_lla, string path2prims, string a_star_result_file);
    void create_init_guess(	vector< dsl::TypedEdge4d* > sol);
    void create_init_guess_with_holes(	vector< dsl::TypedEdge4d* > sol,vector<vector <double> > intersection_positions);
    void create_init_guess_without_reparsing(	vector< dsl::TypedEdge4d* > sol,vector<vector <double> > intersection_positions, std::vector< prim > prims);
    
    vector< vector<double> > solution;
    vector<double> start_pos;
    vector<double> final_pos;
    vector<vector <double> > intersection_positions_w_gaps;
    vector<vector <double> > intersection_positions_wo_gaps;
    double total_length;

};

class NEW_PrimGridPathSearch4d : public Search< TypedCell4d, EdData4d >
{
public:
	
    NEW_PrimGridPathSearch4d(TypedGraph4d& graph_, TypedMap& map_, PathCost4d& cost_,
                     PathGrid4d& grid_, double lat0_, double lon0_, double hgt0_, double d_xy_,
                     Eigen::Vector4d d_nepsi_grid_, double d_occ_map_, int n_trims_, string path2prims, string a_star_result_file);

    virtual ~NEW_PrimGridPathSearch4d();

    /**
     * Set start location in the map
     * @param x euclidean point vector
     * @return true on success
     */
    bool SetStart ( const Eigen::Vector4d& x, VertData4d data_ );

    /**
     * Set goal location in the map
     * @param x euclidean point vector
     * @return true on success
     */
    bool SetGoal ( const Eigen::Vector4d& x, VertData4d data_ );

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
    bool Expand ( TypedVertex4d& from, bool fwd = true );
   
    bool Flow ( TypedCell4d* from_cell, TypedCell4d*& to_cell, EdData4d* EdData4d, int motion_type,
                bool fwd);
    
    /**
     * Crops the input angle to a value between 0 and 2*pi
     * @param the angle that shall be cropped
     * @return psi mod 2pi
     */
    double CropPsi ( double psi );

    /**
     * This function executes a motion primitive and calls depending on the prim_nr either ExecuteTrim or ExecuteManeuver
     * @param v0 the initial state of the object
     * @param vf a 3D-vector to which the final state of the object is written after the Prim was executed
     * @param d for distance this variable carries the spacial length of a primitive 
     * @param prim_id is the ID of the prim and is used to distiguish between trims and maneuvers
     * @param dt is the total time for how long a trim is executed
     * @param fwd is a bool that decides whether the a* search is executed start to goal or goal to start, due to the missing
     *        goal heading fwd is always false.
     * 
     * @return returns a bool that indicates whether a prim is feasible. That includes the geometric constraints and if the cell
     *         was visited before.
     */ 
    template<bool save_path_bool, bool crop_psi_bool>
    bool ExecutePrim( const Eigen::Vector4d& v0, Eigen::Vector4d& vf, double& d, const int prim_nr, const double dt, const bool fwd);
    

    /**
     * This function executes a maneuver. For m intermediate points the position is computed and checked whether they are feasible
     * if all are feasible true is returned false otherwise. The position is interpolated between the corresponding timestamps in the prim struct.
     * @param see docu of @ExecutePrim
     * @return returns a bool that indicates whether a primitive is feasible or not.
     */
    template<bool save_path_bool, bool crop_psi_bool>
    bool ExecuteManeuver( const Eigen::Vector4d& v0, Eigen::Vector4d& vf, double& d, const int man_nr, const bool fwd);
    
     /**
     * This function executes a trims. For m intermediate points the position is computed and checked whether they are feasible
     * if all are feasible true is returned false otherwise. Due to the steady motion of trims their trajectory does not need to be stored,
     * but can be computed.
     * @param see docu of @ExecutePrim
     * @return returns a bool that indicates whether a primitive is feasible or not.
     */
    template<bool save_path_bool, bool crop_psi_bool>
    bool ExecuteTrim( const Eigen::Vector4d& v0, Eigen::Vector4d& vf, double& d, const int trim_nr, const double dt, const bool fwd);
    
    /**
     * SetPrimitives parses all the text files that specify motion primitives and stores them in a the struct prim, it also distinguishes
     * between trim primitives, trim because their behaviour can be specified by just one pair of states and controls, and maneuvers,
     * maneuvers are finite trajectories that are used to concatinate two trims.
     * @param prim_path it the path the directory where all primitives are stored
     */ 
    bool SetPrimitives ( std::string prim_path );
    

    /**
     * An important distinguish that has to be made is the diffrence between the map and the grid. The map only takes the position into account
     * the grid here is a 3D-object which indicates that is also takes the heading into account.
     */
    const PathGrid4d& grid; ///< the grid
    /**
     * the map is used to decide, if a cell is blocked by the habour or free, for that purpose the heading is not needed, therefore it is a 2D-object.
     */ 
	const TypedMap& map; ///<occupancy map from ppm file

    TypedVertex4d** vertexMap; ///< vertex grid array
	
    std::vector< Eigen::Vector4d > visited_pos; //stores all the positions that have been visited already
    std::vector< std::vector<double> > saved_path; //when a solution was found the it is saved to this vector
    
    double lat0; // lat of the local origin
	double lon0; // lon of the local origin
    double hgt0; // height of the local origin
    

    std::vector< prim > prims; // the vector that contains all the parsed primitives
    int n_prims; // total number of primitives
    int n_mans;  // total number of maneuvers
    int n_trims; // total number of trims
    double d_xy; // diameter of a grid cell

    string a_star_result_file; // the path to the file where evaluation data is written to after an a* search.
    
    Eigen::Vector4d d_nepsi_grid;

    TypedCell4d* no_heading_goal_cell;

    double d_occ_map;                                         
    
    int only_steps = 1;
    
    std::map<int, std::vector<int>> trim_man_map_fwd;
    std::map<int, std::vector<int>> trim_man_map_bwd;
	

    bool save_sol_to_file;
};


}

#endif
