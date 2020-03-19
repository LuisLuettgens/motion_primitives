#include "pathgrid.h"

namespace dsl {

using std::vector;

PathGrid2d::PathGrid2d ( const Eigen::Vector2d& xlb, const Eigen::Vector2d& xub, const Vectorni& gs )
    : Grid<Eigen::Vector2d, VertData2d>( xlb, xub, gs )
{}

TypedCell2d* PathGrid2d::Get ( const Eigen::Vector2d& x, bool checkValid) const {
	if ( checkValid )
		if ( !Valid ( x ) ) {
			return 0;
		}
	
	int id = Id ( x );
	TypedCell2d* out = Get ( id );
	if(!out){
		Eigen::Vector2d x_cell = xlb + Eigen::Vector2d ( ( Index ( x, 0 ) + 0.5 ) * cs[0],
														 ( Index ( x, 1 ) + 0.5 ) * cs[1] );
		cells[id] = new TypedCell2d ( id, x_cell );
		cells[id]->data.path_costs = std::numeric_limits<double>::infinity();
		cells[id]->data.nr_of_motions = 0;
		out=cells[id];
	}
	return out;
}

TypedCell2d* PathGrid2d::Get ( int id ) const {
	assert ( id >= 0 );
	if ( id >= nc ) {
		return 0;
	}
	return cells[id];
}





PathGrid3d::PathGrid3d ( const Eigen::Vector3d& xlb, const Eigen::Vector3d& xub, const Vectorni& gs )
	: Grid<Eigen::Vector3d, VertData3d>( xlb, xub, gs )
{}

TypedCell3d* PathGrid3d::Get ( const Eigen::Vector3d& x, bool checkValid) const {
	if ( checkValid )
		if ( !Valid ( x ) ) {
			return 0;
		}

	int id = Id ( x );
	TypedCell3d* out = Get ( id );
	if(!out){
		Eigen::Vector3d x_cell = xlb + Eigen::Vector3d ( ( Index ( x, 0 ) + 0.5 ) * cs[0],
														( Index ( x, 1 ) + 0.5 ) * cs[1],
														( Index ( x, 2 ) + 0.5 ) * cs[2] );
		cells[id] = new TypedCell3d ( id, x_cell );
		cells[id]->data.path_costs = std::numeric_limits<double>::infinity();
		cells[id]->data.nr_of_motions = 0;
		out=cells[id];
	}
	return out;
}

TypedCell3d* PathGrid3d::Get ( int id ) const {
	assert ( id >= 0 );
	if ( id >= nc ) {
		return 0;
	}
	return cells[id];
}

// neu 4D extension

PathGrid4d::PathGrid4d ( const Eigen::Vector4d& xlb, const Eigen::Vector4d& xub, const Vectorni& gs )
	: Grid<Eigen::Vector4d, VertData4d>( xlb, xub, gs )
{}

TypedCell4d* PathGrid4d::Get ( const Eigen::Vector4d& x, bool checkValid) const {
	if ( checkValid )
		if ( !Valid ( x ) ) {
			return 0;
		}

	int id = Id ( x );
	TypedCell4d* out = Get ( id );
	if(!out){
		Eigen::Vector4d x_cell = xlb + Eigen::Vector4d ( ( Index ( x, 0 ) + 0.5 ) * cs[0],
														 ( Index ( x, 1 ) + 0.5 ) * cs[1],
														 ( Index ( x, 2 ) + 0.5 ) * cs[2],
														// neu geadded, vorher null
														 ( Index ( x, 3 ) + 0.5 ) * cs[3]
														 );
		cells[id] = new TypedCell4d ( id, x_cell );
		cells[id]->data.path_costs = std::numeric_limits<double>::infinity();
		cells[id]->data.nr_of_motions = 0;
		out=cells[id];
	}
	return out;
}

TypedCell4d* PathGrid4d::Get ( int id ) const {
	assert ( id >= 0 );
	if ( id >= nc ) {
		return 0;
	}
	return cells[id];
}

}
