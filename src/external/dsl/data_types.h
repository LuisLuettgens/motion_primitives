#ifndef DATA_TYPES2_H
#define DATA_TYPES2_H

#include "edge.h"
#include "vertex.h"
#include "graph.h"
#include "cell.h"
#include "grid.h"
#include "map.h"
#include <Eigen/Geometry>

namespace dsl
{

// 2d case: simple path search
struct EdData2d{
	int from_cell_id;
	int to_cell_id;
	int motion_type;
	double costs;
};
struct VertData2d{
	double path_costs;
	int nr_of_motions;
	Vertex<Cell<Eigen::Vector2d, VertData2d>,EdData2d>* vert;
};
using TypedCell2d = Cell<Eigen::Vector2d, VertData2d>;
using TypedEdge2d = Edge<TypedCell2d,EdData2d>;
using TypedVertex2d = Vertex<TypedCell2d,EdData2d>;
using TypedGraph2d = Graph<TypedCell2d,EdData2d>;
using Vector2i =  Eigen::Matrix< int, 2, 1 >;

// 3d case: path search with straight lines and arcs and one dim primitive search
struct EdData3d{
	int from_cell_id;
	int to_cell_id;
	int motion_type;
	double costs;
};
struct VertData3d{
	double path_costs;
	int nr_of_motions;
	Vertex<Cell<Eigen::Vector3d, VertData3d>,EdData3d>* vert;
};
using TypedCell3d = Cell<Eigen::Vector3d, VertData3d>;
using TypedEdge3d = Edge<TypedCell3d,EdData3d>;
using TypedVertex3d = Vertex<TypedCell3d,EdData3d>;
using TypedGraph3d = Graph<TypedCell3d,EdData3d>;
using Vector3i =  Eigen::Matrix< int, 3, 1 >;

// 4d case:multi dimensional prims
struct EdData4d{
	int from_cell_id;
	int to_cell_id;
	int motion_type;
	double costs;
};
struct VertData4d{
	double path_costs;
	int nr_of_motions;
	Vertex<Cell<Eigen::Vector4d, VertData4d>,EdData4d>* vert;
};
using TypedCell4d = Cell<Eigen::Vector4d, VertData4d>;
using TypedEdge4d = Edge<TypedCell4d,EdData4d>;
using TypedVertex4d = Vertex<TypedCell4d,EdData4d>;
using TypedGraph4d = Graph<TypedCell4d,EdData4d>;
using Vector4i =  Eigen::Matrix< int, 4, 1 >;

// occupancy map is always in 2d
using TypedMap = Map<bool, 2>;


}

#endif