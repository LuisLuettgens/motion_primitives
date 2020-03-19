#include "gridpathsearch.h"

namespace dsl
{

using std::vector;

SimplePathSearch::SimplePathSearch(GeometricConstr& geom_constr, Eigen::Vector2d d_xy_ned_grid, Eigen::Vector2d start_lla, Eigen::Vector2d goal_lla)
{
	
	dsl::TypedGraph2d graph;
	
	double x_l_ned = geom_constr.occ_map->xlb[0];
	double x_u_ned = geom_constr.occ_map->xub[0];
	double y_l_ned = geom_constr.occ_map->xlb[1];
	double y_u_ned = geom_constr.occ_map->xub[1];
	
	int n_x = floor( (x_u_ned-x_l_ned)/d_xy_ned_grid[0] );
	int n_y = floor( (y_u_ned-y_l_ned)/d_xy_ned_grid[1] );
	double d_occ_map=geom_constr.occ_map_cell_dim_ned[0];
	if(d_occ_map>geom_constr.occ_map_cell_dim_ned[1]) d_occ_map=geom_constr.occ_map_cell_dim_ned[1];
	
	dsl::PathGrid2d grid( Eigen::Vector2d(x_l_ned, y_l_ned), Eigen::Vector2d(x_u_ned, y_u_ned), dsl::Vector2i(n_x, n_y) );
	
	dsl::PathCost2d cost = dsl::PathCost2d();
	
	dsl::GridPathSearch2d search(graph, *(geom_constr.occ_map), cost, grid, geom_constr.lat0, geom_constr.lon0, geom_constr.hgt0, d_occ_map);
	
	
	dsl::VertData2d vert_data;
	vert_data.path_costs = std::numeric_limits<double>::infinity();
	vert_data.nr_of_motions = 0;
	

	Eigen::Vector2d start_ned(0.,0.);
	Eigen::Vector2d goal_ned(0.,0.);
	double dummy=0.;
	LLAtoNED(start_lla[0], start_lla[1], 0., geom_constr.lat0, geom_constr.lon0, geom_constr.hgt0, &start_ned[0], &start_ned[1], &dummy);
	start_ned[2]=start_lla[2];
	LLAtoNED(goal_lla[0], goal_lla[1], 0., geom_constr.lat0, geom_constr.lon0, geom_constr.hgt0, &goal_ned[0], &goal_ned[1], &dummy);
	goal_ned[2]=goal_lla[2];
	// Goal has to be set before Start
	if( !search.SetGoal(goal_ned, vert_data) ) exit(1);
	if( !search.SetStart(start_ned, vert_data) ) exit(1);
	
	vector< dsl::TypedEdge2d* > sol;
	
	struct timeval timer;
	long time = dsl::timer_us ( &timer );
	dsl::timer_start ( &timer );
	search.Plan(sol);
	time = dsl::timer_us ( &timer );
	printf ( "plan path time= %ld  us\n", time );
	cout << "Graph has " << search.Vertices() << " vertices and " << search.Edges() << " edges. " << endl;

	if(search.visited_pos.size()>0)
	{
		string filename2 = "visited_positions.m";
		ofstream myfile2 (filename2);
		if (!myfile2.is_open())
		{
			cout << "Unable to open file: " << filename2 << endl;
			exit(1);
		}
		myfile2 << setprecision(16);
		myfile2 << "vis_pos=[";
		for(unsigned int i=0;i<search.visited_pos.size();i++)
		{
			myfile2<<search.visited_pos[i][0]<<" "<<search.visited_pos[i][1]<<";\n";
		}
		myfile2 << "];";
		myfile2.close();
	}

	int n=sol.size();
	if(n>0)
	{
		solution.resize(0);
		solution.push_back(Vec2d(sol[n-1]->to->data.c[0], sol[n-1]->to->data.c[1]));
		for(int i=n-1; i>=0; i--)
		{
			solution.push_back(Vec2d(sol[i]->from->data.c[0], sol[i]->from->data.c[1]));
		}
		
		if(search.visited_pos.size()>0)
		{
			string filename1 = "solution.m";
			ofstream myfile1 (filename1);
			if (!myfile1.is_open())
			{
				cout << "Unable to open file: " << filename1 << endl;
				exit(1);
			}
			myfile1 << setprecision(16);
			myfile1 << "sol=["<<sol[n-1]->to->data.c[0]<<" "<<sol[n-1]->to->data.c[1]<<";\n";
			for(int i=n-1; i>=0; i--)
			{
				myfile1<<sol[i]->from->data.c[0]<<" "<<sol[i]->from->data.c[1]<<";\n";
			}
			myfile1 << "];";
			myfile1.close();
		}
	}else
	{
		cout << "Warning: No solution found!" << endl;
	}

};





GridPathSearch2d::GridPathSearch2d (TypedGraph2d& graph_, TypedMap& map_, PathCost2d& cost_, PathGrid2d& grid_,double lat0_, double lon0_, double hgt0_, double d_occ_map_)
		 :Search< TypedCell2d, EdData2d > ( graph_, cost_ ),
		  grid(grid_),
		  map(map_),
		  lat0(lat0_),
		  lon0(lon0_),
		  hgt0(hgt0_),
		  d_occ_map(d_occ_map_),
		  save_sol_to_file(false)
{

    vertexMap = new TypedVertex2d* [grid.nc];
    memset ( vertexMap, 0, grid.nc * sizeof ( TypedVertex2d* ) );
	
}

GridPathSearch2d::~GridPathSearch2d()
{

	for ( auto edge : graph.edges )
	{
		delete edge.second;
    }
	for ( auto vertex : graph.vertices )
	{
	    delete vertex.second;
    }

    delete[] vertexMap;
}

bool GridPathSearch2d::Expand ( TypedVertex2d& from, bool fwd )
{
    // if this is true then all vertices have already been expanded at
    // construction
	if ( fwd && from.succExpanded )
	{
        return true;
    }

	if ( !fwd && from.predExpanded ) 
	{
        return true;
    }

    // cell must exist
    TypedCell2d* from_cell = grid.Get ( from.data.id );
    assert ( from_cell );

	// there 8 possible motions (for each neighbor-cell):
	//	0: north-west
	//	1: north
	//	2: north-east
	//	3: west
	//	4: east
	//	5: south-west
	//	6: south
	//	7: south-east
	for ( unsigned int motion_type=0; motion_type<8; motion_type++ )
	{

		TypedCell2d *to_cell;
		EdData2d* eddata = new EdData2d();
		if ( !Flow ( from_cell, to_cell, eddata, motion_type, fwd ) )
		{
			continue;
		}
		assert ( to_cell );

		int id = to_cell->id;

		// if this vertex doesn't exist, create it and add to graph
		if ( !vertexMap[id] ) 
		{
			vertexMap[id] = new TypedVertex2d ( *to_cell );
			graph.AddVertex ( *vertexMap[id] );
		}
		TypedVertex2d* to = vertexMap[id];

		// if fwd and incoming edge from->to exists, then do not create a new one
		if ( from.Find ( *to, !fwd ) )
		{
			continue;
		}

		TypedEdge2d* edge = fwd ? new TypedEdge2d ( *eddata, &from, to, eddata->costs ) : new TypedEdge2d ( *eddata, to, &from, eddata->costs );
		graph.AddEdge ( *edge );
		
		if(save_sol_to_file) visited_pos.push_back(to->data.c);
		
    }
    
	if(Vertices()>1.0E6)
	{
		std::cout << "Search exceeds search limit: abort" << std::endl;
		
		if(visited_pos.size()>0)
		{
			std::string filename2 = "visited_positions.m";
			std::ofstream myfile2 (filename2);
			if (!myfile2.is_open())
			{
				std::cout << "Unable to open file: " << filename2 << std::endl;
				exit(1);
			}
			myfile2 << std::setprecision(16);
			myfile2 << "vis_pos=[";
			for(unsigned int i=0;i<visited_pos.size();i++)
			{
				myfile2<<visited_pos[i][0]<<" "<<visited_pos[i][1]<<";\n";
			}
			myfile2 << "];";
			myfile2.close();
		}
		
		exit(1);
	}

	if ( fwd )
	{
        from.succExpanded = true;
	} else
	{
        from.predExpanded = true;
    }

    return true;
}

bool GridPathSearch2d::Flow ( TypedCell2d* from_cell, TypedCell2d*& to_cell, EdData2d* eddata, int motion_type, bool fwd )
{
	
	Eigen::Vector2d v0(from_cell->c);
	Eigen::Vector2d d_v(0., 0.);
	
	// choose motion via i
	switch ( motion_type ){
		//	0: north-west
		case 0:
			d_v[0]=-grid.cs[0];
			d_v[1]= grid.cs[1];
			break;
		//	1: north
		case 1:
			d_v[0]=0.;
			d_v[1]=grid.cs[1];
			break;
		//	2: north-east
		case 2:
			d_v[0]=grid.cs[0];
			d_v[1]=grid.cs[1];
			break;
		//	3: west
		case 3:
			d_v[0]=-grid.cs[0];
			d_v[1]= 0.;
			break;
		//	4: east
		case 4:
			d_v[0]=grid.cs[0];
			d_v[1]=0.;
			break;
		//	5: south-west
		case 5:
			d_v[0]=-grid.cs[0];
			d_v[1]=-grid.cs[1];
			break;
		//	6: south
		case 6:
			d_v[0]= 0.;
			d_v[1]=-grid.cs[1];
			break;
		//	7: south-east
		case 7:
			d_v[0]= grid.cs[0];
			d_v[1]=-grid.cs[1];
			break;
		default:
			std::cout<<"Error inGridPathSearch2d::Flow: motion "<<motion_type<<" not implemented"<<std::endl;
			assert(false);
	}
	//dsl searches backwards!
	if(fwd) d_v=-1.*d_v;
	
	Eigen::Vector2d vf(0., 0.);

	double dist=sqrt(d_v[0]*d_v[0]+d_v[1]*d_v[1]);

	int m = ceil(dist/d_occ_map);
	for (int j=1; j<m; j++)
	{
		vf=v0+((double)j)/((double)m)*d_v;
		if( map.Get(Eigen::Vector2d(vf[0],vf[1])) ) return false;
	}
	
	vf=v0+d_v;
	if( map.Get(Eigen::Vector2d(vf[0],vf[1])) ) return false;
	to_cell = grid.Get ( vf );
	if(!to_cell) return false;
	if( map.Get(Eigen::Vector2d(to_cell->c[0],to_cell->c[1])) ) return false;
	
	eddata->costs = dist;//cost.EdgeCost(*to_cell, *from_cell);
	
	if(to_cell->data.path_costs > from_cell->data.path_costs+eddata->costs)
	{
		to_cell->data.path_costs=from_cell->data.path_costs+eddata->costs;
		to_cell->data.nr_of_motions=from_cell->data.nr_of_motions+1;
	}
	
	eddata->motion_type = motion_type;
	
	if(fwd)
	{
		eddata->from_cell_id = from_cell->id;
		eddata->to_cell_id = to_cell->id;
	}else
	{
		eddata->from_cell_id = to_cell->id;
		eddata->to_cell_id = from_cell->id;
	}

    return true;
}

bool GridPathSearch2d::SetStart ( const Eigen::Vector2d& x, VertData2d data_ )
{
	if ( !grid.Valid ( x ) )
	{
        std::cout << "[W] GridPathSearch2d:SetStart: invalid x=" << x.transpose()
                  << std::endl;
        return false;
    }

    TypedCell2d* cell = grid.Get ( x );
	if ( !cell )
	{
        std::cout << "[W] GridPathSearch2d:SetStart: cell at=" << x.transpose() << " does not exist!"
                  << std::endl;
        return false;
    }
	
	if( map.Get(x) )
	{
		std::cout << "[W] GridPathSearch2d:SetStart: Starting Point: " << x.transpose() <<
					 " is occupied!"
                  << std::endl;
        return false;
    }

    int id = cell->id;
	assert ( id >= 0 && id < grid.nc );
	
	cell->data = data_;

    // if it's not added previously add it
	if ( !vertexMap[id] ) 
	{
        vertexMap[id] = new TypedVertex2d ( *cell );
        graph.AddVertex ( *vertexMap[id] );
    }

	//DSL serches backwards!
    Search< TypedCell2d, EdData2d >::SetGoal ( *vertexMap[id] );

    return true;
}

bool GridPathSearch2d::SetGoal ( const Eigen::Vector2d& x, VertData2d data_ )
{
    if ( !grid.Valid ( x ) ) {
        std::cout << "[W] GridPathSearch2d:SetGoal: invalid x=" << x.transpose()
                  << std::endl;
        return false;
    }

    TypedCell2d* cell = grid.Get ( x );
    if ( !cell ) {
        std::cout << "[W] GridPathSearch2d:SetGoal: cell at=" << x.transpose() << " does not exist!"
                  << std::endl;
        return false;
	}
	
	if( map.Get(x) ) {
		std::cout << "[W] GridPathSearch2d:SetGoal: Goal Point: " << x.transpose() <<
					 " is occupied!"
                  << std::endl;
        return false;
    }

	int id = cell->id;
    assert ( id >= 0 && id < grid.nc );
	
	cell->data = data_;

    // if it's not added previously add it
    if ( !vertexMap[id] ) {
        vertexMap[id] = new TypedVertex2d ( *cell );
        graph.AddVertex ( *vertexMap[id] );
    }

	//DSL serches backwards!
    Search< TypedCell2d, EdData2d >::SetStart ( *vertexMap[id] );

    return true;
}


bool PrimGridPathSearch3d::SetPrimitives ( std::string prim_path ) {

	// read all trim files in prim_path
	struct dirent **namelist;
	std::vector<std::string> file_names;
	int n = scandir(prim_path.c_str(), &namelist, 0, versionsort);
	if (n < 0) perror(prim_path.c_str());
	else{
		for(int i =0 ; i < n; ++i){
			file_names.push_back( namelist[i]->d_name );
			free(namelist[i]);
		}
		free(namelist);
	}
	n_mans=0;
	n_trims=0;
	n_prims=0;
	for (unsigned int i=0; i<file_names.size();) {
		std::size_t found_man = file_names[i].find("man");
		std::size_t found_tp = file_names[i].find("tp");
		if (found_man!=std::string::npos){
			n_mans++;
			n_prims++;
			i++;
		}else if(found_tp!=std::string::npos){
			std::string tmp = file_names[i];
			file_names.erase(file_names.begin()+i);
			file_names.insert(file_names.begin()+n_trims, tmp);
			n_trims++;
			n_prims++;
			i++;
		}else{
			file_names.erase(file_names.begin()+i);
		}
	}

	nmbr_motion_primitives = n_mans+n_trims;
	
	// read the data from files
	prim_links_fwd.resize(n_trims);
	prim_links_bwd.resize(n_trims);
	result_trim_fwd.resize(n_prims);
	result_trim_bwd.resize(n_prims);
	
	for (int i=0; i<n_trims;i++) {
		prim pr;
		pr.is_trim=true;
		pr.n_t_steps=0;
		pr.costs=0.;
		std::ifstream f;
		std::string str;
		string namef = prim_path+file_names[i];
			
		f.open ( prim_path+file_names[i].c_str(), std::ios::in );
		if ( !f ) {
			throw std::runtime_error ( "Error Reading File in ShipConnectivity::SetPrimitives; Wrong Filename/Path?" );
		}
		getline ( f,str );
		std::istringstream sin(str);
		double tmp;
		if( sin>>tmp ){
			pr.id=(int)tmp-1;
		}
		getline ( f,str );
		sin.clear();
		sin.str(str);
		if( sin>>tmp ){
			pr.t.push_back(tmp);
		}
		if( sin>>tmp ){
			pr.x.push_back(tmp);
		}
		if( sin>>tmp ){
			pr.y.push_back(tmp);
		}
		if( sin>>tmp ){
			pr.psi.push_back(tmp);
		}
		if( sin>>tmp ){
			pr.u.push_back(tmp);
		}
		if( sin>>tmp ){
			pr.v.push_back(tmp);
		}
		if( sin>>tmp ){
			pr.r.push_back(tmp);
		}
		//neu
		if( sin>>tmp ){
			pr.delta.push_back(tmp);
		}
		if( sin>>tmp ){
			pr.n.push_back(tmp);
		}
		if( sin>>tmp ){
			pr.energy.push_back(tmp);
		}
		if( sin>>tmp ){
			pr.d_delta.push_back(tmp);
		}
		if( sin>>tmp ){
			pr.d_n.push_back(tmp);
		}
		pr.d.push_back(0);
		prims.push_back(pr);
	}
	for (int i=n_trims; i<n_prims;i++) {
		prim pr;
		pr.is_trim=false;
		std::ifstream f;
		std::string str;
		f.open ( prim_path+file_names[i].c_str(), std::ios::in );
		if ( !f ) {
			throw std::runtime_error ( "Error Reading File in ShipConnectivity::SetPrimitives; Wrong Filename/Path?" );
		}
		getline ( f,str );
		std::istringstream sin(str);
		double tmp;
		if( sin>>tmp ){
			pr.id=(int)tmp-1;
		}
		
		getline ( f,str );
		sin.clear();
		sin.str(str);
		if( sin>>tmp ){
			pr.n_t_steps=(int)tmp;
		}
		if( sin>>tmp ){
			pr.costs=tmp;
		}
		
		getline ( f,str );
		while ( !f.eof()){
			sin.clear();
			sin.str(str);
			if( sin>>tmp ){
				pr.t.push_back(tmp);
			}
			if( sin>>tmp ){
				pr.x.push_back(tmp);
			}
			if( sin>>tmp ){
				pr.y.push_back(tmp);
			}
			if( sin>>tmp ){
				pr.psi.push_back(tmp);
			}
			if( sin>>tmp ){
				pr.u.push_back(tmp);
			}
			if( sin>>tmp ){
				pr.v.push_back(tmp);
			}
			if( sin>>tmp ){
				pr.r.push_back(tmp);
			}
			//neu
			
			if( sin>>tmp ){
				pr.delta.push_back(tmp);
			}
			if( sin>>tmp ){
				pr.n.push_back(tmp);
			}
			if( sin>>tmp ){
			pr.energy.push_back(tmp);
			}
			if( sin>>tmp ){
				pr.d_delta.push_back(tmp);
			}
			if( sin>>tmp ){
				pr.d_n.push_back(tmp);
			}
			getline ( f,str );
		}
		pr.d.push_back(0);
		prims.push_back(pr);
/*		
		// set prim_connecs
		int man = atoi(file_names[i].substr(4,2).c_str())-1;
		int from = atoi(file_names[i].substr(7,2).c_str())-1;
		int to = atoi(file_names[i].substr(10,2).c_str())-1;
// 		prim_links_fwd[man].push_back(to);
		prim_links_fwd[from].push_back(man);
		prim_links_bwd[to].push_back(man);
// 		prim_links_bwd[man].push_back(from);
		result_trim_fwd[man]=to;
		result_trim_bwd[man]=from;*/
	}
// 	for (int i=n_trims; i<n_prims;i++) {
// 		int tr=prim_links_fwd[i][0];
// 		for (unsigned int j=0; j<prim_links_fwd[tr].size(); j++) {
// 			prim_links_fwd[i].push_back(prim_links_fwd[tr][j]);
// 		}
// 		tr=prim_links_bwd[i][0];
// 		for (unsigned int j=0; j<prim_links_bwd[tr].size(); j++) {
// 			prim_links_bwd[i].push_back(prim_links_bwd[tr][j]);
// 		}
// 	}
/*
	for (int i=0; i<n_trims;i++) {
		prim_links_fwd[i].insert(prim_links_fwd[i].begin(),i);
		prim_links_bwd[i].insert(prim_links_bwd[i].begin(),i);
		result_trim_fwd[i]=i;
		result_trim_bwd[i]=i;
	}
*/	
	for (int i=n_trims; i<n_prims;i++) {
		prims[i].d.resize(prims[i].x.size());
		for (unsigned int j=1; j<prims[i].d.size(); j++) {
			prims[i].d[j] = prims[i].d[j-1]
							 + sqrt( (prims[i].x[j]-prims[i].x[j-1])*(prims[i].x[j]-prims[i].x[j-1])
								+ (prims[i].y[j]-prims[i].y[j-1])*(prims[i].y[j]-prims[i].y[j-1]) );
		}
	}
	
/*
	double v_max=0.;
	double r_max=0.;
	double psi_max=0.;
	double d_max=0.;
	double radius_min=std::numeric_limits<double>::infinity();
	
	for (int i=0; i<n_trims;i++) { 
		double v = sqrt(prims[i].u[0]*prims[i].u[0]+prims[i].v[0]*prims[i].v[0]);
		if(fabs(v)>v_max) v_max=fabs(v);
		if(fabs(prims[i].r[0])>r_max) r_max=fabs(prims[i].r[0]);
		double radius = v/fabs(prims[i].r[0]);
		if(radius<radius_min) radius_min=radius;
	}
	for (unsigned int i=0; i<dt_vec.size();i++) {
		if(v_max*dt_vec[i]>d_max) d_max=v_max*dt_vec[i];
		if(r_max*dt_vec[i]>psi_max) psi_max=r_max*dt_vec[i];
	}
	for (int i=n_trims; i<n_prims;i++) {
		if(prims[i].d.back()>d_max) d_max=prims[i].d.back();
		if(prims[i].psi.back()>psi_max) psi_max=prims[i].psi.back();
		double radius = 0.;
		double v_mean = 0.;
		double r_mean = 0.;
		for (unsigned int j=0; j<prims[i].u.size(); j++) {
			double v = sqrt(prims[i].u[j]*prims[i].u[j]+prims[i].v[j]*prims[i].v[j]);
			double r = fabs(prims[i].r[j]);
			v_mean+=v;
			r_mean+=r;
			radius += v/r;
		}
		radius=radius/prims[i].u.size();
		if(radius<radius_min){
			radius_min=radius;
		}
		v_mean=v_mean/prims[i].u.size();
		if(fabs(v_mean)>v_max) v_max=fabs(v_mean);
		r_mean=r_mean/prims[i].u.size();
		if(fabs(r_mean)>r_max) r_max=fabs(r_mean);
	}
	
	if( *(cost.dubin_mode)==2) radius_min -= 70.;
	
	*(cost.radius_min) = radius_min;
	*(cost.v_max) = v_max;
	*(cost.d_max) = d_max;
	*(cost.psi_max) = psi_max;
	*(cost.r_max) = r_max;
	
	*(cost.prims) = prims;

	bool is_trim;
	int id;
	int n_t_steps;
	double costs;
	std::vector<double> t;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> psi;
	std::vector<double> u;
	std::vector<double> v;
	std::vector<double> r;
	std::vector<double> d;


	*/

	return true;
}

PrimPathSearch::PrimPathSearch(GeometricConstr& geom_constr, Eigen::Vector3d d_nepsi_grid, Eigen::Vector3d start_llpsi, Eigen::Vector3d goal_llpsi)
{
	dsl::TypedGraph3d graph;

	double x_l_ned = geom_constr.occ_map->xlb[0];
	double x_u_ned = geom_constr.occ_map->xub[0];
	double y_l_ned = geom_constr.occ_map->xlb[1];
	double y_u_ned = geom_constr.occ_map->xub[1];
	// eps = 1.E-9 wegen grid.valid check bei winkeln
	double psi_l = 0.-1.E-9;
	double psi_u = 2*M_PI+1.E-9;
	double d_occ_map=geom_constr.occ_map_cell_dim_ned[0];
	if(d_occ_map>geom_constr.occ_map_cell_dim_ned[1]) d_occ_map=geom_constr.occ_map_cell_dim_ned[1];
	int n_x = floor( (x_u_ned-x_l_ned)/d_nepsi_grid[0] );
	int n_y = floor( (y_u_ned-y_l_ned)/d_nepsi_grid[1] );
	int n_psi = floor( (psi_u-psi_l)/d_nepsi_grid[2] );
	//neu map statt grid
	double d_xy = sqrt(geom_constr.occ_map->cs[0]*geom_constr.occ_map->cs[0]+geom_constr.occ_map->cs[1]*geom_constr.occ_map->cs[1]);
	dsl::PathGrid3d grid( Eigen::Vector3d(x_l_ned, y_l_ned, psi_l),
						  Eigen::Vector3d(x_u_ned, y_u_ned, psi_u),
						  dsl::Vector3i(n_x, n_y, n_psi) );

	dsl::PathCost3d cost = dsl::PathCost3d();
	PrimGridPathSearch3d search(graph, *(geom_constr.occ_map), cost, grid,
								 geom_constr.lat0, geom_constr.lon0, geom_constr.hgt0, d_xy,
								 d_nepsi_grid, d_occ_map, 1);
	dsl::VertData3d vert_data;
	vert_data.path_costs = std::numeric_limits<double>::infinity();
	vert_data.nr_of_motions = 0;
	

	double dummy=0.;
	Eigen::Vector3d start_nepsi(0.,0.,0.);
	Eigen::Vector3d goal_nepsi(0.,0.,0.);
	LLAtoNED(start_llpsi[0], start_llpsi[1], 0., geom_constr.lat0, geom_constr.lon0, geom_constr.hgt0,
			 &start_nepsi[0], &start_nepsi[1], &dummy);
	start_nepsi[2]=start_llpsi[2];
	LLAtoNED(goal_llpsi[0], goal_llpsi[1], 0., geom_constr.lat0, geom_constr.lon0, geom_constr.hgt0,
			 &goal_nepsi[0], &goal_nepsi[1], &dummy);
	goal_nepsi[2]=goal_llpsi[2];
	// Goal has to be set before Start
	if( !search.SetGoal(goal_nepsi, vert_data) ) exit(1);
	if( !search.SetStart(start_nepsi, vert_data) ) exit(1);



	vector< dsl::TypedEdge3d* > sol;

	struct timeval timer;
	long time = dsl::timer_us ( &timer );
	dsl::timer_start ( &timer );
	search.Plan(sol);
	cout << "sol size:" << sol.size() << endl;
	time = dsl::timer_us ( &timer );
	cout << "d_xy main: " << d_xy << endl;
	printf ( "A* search plan time= %ld  us\n", time );
	cout << "Graph has " << search.Vertices() << " vertices and " << search.Edges() << " edges. " << endl;
	double cost_sum=0.;
	for(int i=0; i<sol.size(); i++){
		cost_sum+=sol[i]->data.costs;
	}
	cout << "Final costs=" << cost_sum <<endl;
	int visited_pos_size = search.visited_pos.size();
	if(visited_pos_size>0)
	{
		string filename2 = "visited_positions.m";
		ofstream myfile2 (filename2);
		if (!myfile2.is_open())
		{
			cout << "Unable to open file: " << filename2 << endl;
			exit(1);
		}
		myfile2 << setprecision(16);
		myfile2 << "vis_pos=[";
		for(unsigned int i=0;i<visited_pos_size;i++)
		{
			myfile2<<search.visited_pos[i][0]<<" "<<search.visited_pos[i][1]<<" "<<search.visited_pos[i][2]<<";\n";
		}
		myfile2 << "];";
		myfile2.close();
		start_pos = {search.visited_pos[0][0] , search.visited_pos[0][1] , search.visited_pos[0][2] };
		final_pos = {search.visited_pos[visited_pos_size-1][0] , search.visited_pos[visited_pos_size-1][1] , search.visited_pos[visited_pos_size-1][2] };
	}

	int n=sol.size();
	if(n>0)
	{
		intersection_positions_w_gaps.resize(n);
		intersection_positions_wo_gaps.resize(n);
		vector<Eigen::Vector3d> vec;
		double dist=0.;
		search.saved_path.resize(0);
		search.saved_path.push_back({0., sol[n-1]->to->data.c[0],sol[n-1]->to->data.c[1],sol[n-1]->to->data.c[2]});
		Eigen::Vector3d v0_nepsi_w_gaps(sol[n-1]->to->data.c[2],sol[n-1]->to->data.c[0],sol[n-1]->to->data.c[1]);
		Eigen::Vector3d vf_nepsi_w_gaps(0.,0.,0.);
		Eigen::Vector3d v0_nepsi_wo_gaps(		sol[n-1]->to->data.c[2],sol[n-1]->to->data.c[0],sol[n-1]->to->data.c[1]);
		//Eigen::Vector3d v0_nepsi_wo_gaps(start_nepsi[2],start_nepsi[0],start_nepsi[1]);
		Eigen::Vector3d vf_nepsi_wo_gaps(0.,0.,0.);
		intersection_positions_w_gaps[n-1] = {v0_nepsi_w_gaps[1],v0_nepsi_w_gaps[2],v0_nepsi_w_gaps[0]};
		intersection_positions_wo_gaps[n-1] = {v0_nepsi_wo_gaps[1],v0_nepsi_wo_gaps[2],v0_nepsi_wo_gaps[0]};
		
		// Display of the stating points of each prim
		/*
		cout << "#########################" << endl;
		cout << "prim id: " << sol[n-1]->data.motion_type << endl;
			cout << "x: " << setw(8) << intersection_positions_wo_gaps[n-1][0]
				 <<  " y: " << setw(8) << intersection_positions_wo_gaps[n-1][1]
				 << "  psi: " << setw(8) << intersection_positions_wo_gaps[n-1][2]
				 << endl;
		*/
		for(int i=n-1; i>=1; i--)
		{
			//without gaps
			search.ExecutePrim<false,false>( v0_nepsi_wo_gaps, vf_nepsi_wo_gaps,  dist, sol[i]->data.motion_type, dt, false);
			v0_nepsi_wo_gaps=vf_nepsi_wo_gaps;
			intersection_positions_wo_gaps[i-1] = {v0_nepsi_wo_gaps[1],v0_nepsi_wo_gaps[2],v0_nepsi_wo_gaps[0]};
			
			//Display of the stating and the end points of each prim
			/*	
				cout << "#########################" << endl;
				cout << "without gaps:" << endl;
				cout << "prim id: " << sol[i-1]->data.motion_type << endl;
				cout << "x_vf: " << setw(8) << vf_nepsi_wo_gaps[1]
				 <<  " y_vf: " << setw(8) << vf_nepsi_wo_gaps[2]
				 << "  psi_vf: " << setw(8) << vf_nepsi_wo_gaps[0]
				 << endl;
				 cout << "x_v0: " << setw(8) << intersection_positions_wo_gaps[i-1][0]
				 <<  " y_v0: " << setw(8) << intersection_positions_wo_gaps[i-1][1]
				 << "  psi_v0: " << setw(8) << intersection_positions_wo_gaps[i-1][2]
				 << endl;
			*/	

			//with gaps
			search.ExecutePrim<true,false>( v0_nepsi_w_gaps, vf_nepsi_w_gaps, dist, sol[i]->data.motion_type, dt, false);
			
			double psi_uncropped = vf_nepsi_w_gaps[0];
			double psi_cropped = search.CropPsi(psi_uncropped);
			TypedCell3d* vf_cell = grid.Get ( Eigen::Vector3d(vf_nepsi_w_gaps[1],vf_nepsi_w_gaps[2],psi_cropped) );
			if(!vf_cell)
			{
				 cout<<"Fatal Error: Cell doesn't exist"<<endl; exit(1);
			}
			v0_nepsi_w_gaps=Eigen::Vector3d(vf_cell->c[2]+(psi_uncropped-psi_cropped),vf_cell->c[0],vf_cell->c[1]);
			intersection_positions_w_gaps[i-1] = {v0_nepsi_w_gaps[1],v0_nepsi_w_gaps[2],v0_nepsi_w_gaps[0]};
			

			 //Display of the stating and the end points of each prim
			/*
			 cout << "#########################" << endl;
				cout << "with gaps:" << endl;
				cout << "prim id: " << sol[i-1]->data.motion_type << endl;
				cout << "x_vf: " << setw(8) << vf_nepsi_w_gaps[1]
				 <<  " y_vf: " << setw(8) << vf_nepsi_w_gaps[2]
				 << "  psi_vf: " << setw(8) << vf_nepsi_w_gaps[0]
				 << endl;
				 cout << "x_v0: " << setw(8) << intersection_positions_w_gaps[i-1][0]
				 <<  " y_v0: " << setw(8) << intersection_positions_w_gaps[i-1][1]
				 << "  psi_v0: " << setw(8) << intersection_positions_w_gaps[i-1][2]
				 << endl;
			*/
		}

		solution=search.saved_path;

		if(search.visited_pos.size()>0)
		{
			//create_init_guess(sol);
			//create_init_guess_with_holes(sol,intersection_positions_w_gaps);
			create_init_guess_without_reparsing(sol,intersection_positions_w_gaps, search.prims);
			
			string filename1 = "solution.m";
			ofstream myfile1 (filename1);
			if (!myfile1.is_open())
			{
				cout << "Unable to open file: " << filename1 << endl;
				exit(1);
			}
			myfile1 << "sol=[";
			for(int i=0; i<search.saved_path.size(); i++)
			{
				myfile1<<search.saved_path[i][1]<<" " <<search.saved_path[i][2]<<" "<<search.saved_path[i][3]<<";\n";
			}
			myfile1 << "];";
			myfile1.close();

			string filename2 = "intersection_points.m";
			ofstream myfile2 (filename2);
			if (!myfile2.is_open())
			{
				cout << "Unable to open file: " << filename1 << endl;
				exit(1);
			}
			myfile2.close();
		}

	}else{
		cout << "Warning: No solution found!" << endl;
	}

};

void PrimPathSearch::create_init_guess(	vector< dsl::TypedEdge3d* > sol){
	struct dirent **namelist;
	string path_to_primitives = "../bin/motion_prims_3_cropped/";
	std::vector<std::string> file_names;
	int n = scandir(path_to_primitives.c_str(), &namelist, 0, versionsort);
	if (n < 0) perror(path_to_primitives.c_str());
	else{
		for(int i =0 ; i < n; ++i){
			file_names.push_back( namelist[i]->d_name );
			free(namelist[i]);
		}
		free(namelist);
	}
	// file to which the whole trajectory shall be written
	ofstream myfile;
	myfile.open ("example.txt");

	//man that stores the maneuvertime to a given id 
	map<int,double> maneuvertimes;
	map<int,vector <vector<double>>> motion_prims;
	double t, x, y, psi, u, v, r, delta, rudder_rpm, d_delta, d_rudder_rpm;
					
	// variable to  storer the entire time of the trajectory
	double total_time = 0.0;
	for(int i = sol.size()-1; i>=0;--i){
		string prim_type = (sol[i]->data.motion_type== 0)? "tp_01": "man_0"+to_string(sol[i]->data.motion_type);
		for(auto file: file_names){
			if (file.find(prim_type) == string::npos){
				continue;
			}
			else{
				
				// creating a map containing all occurring motion primitives for this trajectory
				if(motion_prims.find(sol[i]->data.motion_type)== motion_prims.end()){
					ifstream infile;
				
					//cout << path_to_primitives+file << endl;
					infile.open (path_to_primitives+file);	
					while (infile >> t >> x >> y >> psi >> u >> v >> r >> delta >> rudder_rpm >> d_delta >> d_rudder_rpm)
					{
						motion_prims[sol[i]->data.motion_type].push_back({t,x,y,psi,u,v,r,delta, rudder_rpm, d_delta, d_rudder_rpm});
					}
				}
				// chooping of the data-type
				const std::string ext(".txt");
				if ( file != ext && file.size() > ext.size() && file.substr(file.size() - ext.size()) == ".txt" )
				{
   					// if so then strip them off
					   file = file.substr(0, file.size() - ext.size());
				}
				// extract maneuver time from file name
				int underscore_counter =0;
				if(file.substr(0,3)=="man"){
					for(int j=0;j < file.length();++j){
						if(underscore_counter==7){
							total_time += stod(file.substr(j, file.size()));
							maneuvertimes[sol[i]->data.motion_type] =  stod(file.substr(j, file.size()));
							break;
						}
						if(file[j]=='_'){
							underscore_counter++;
						}
					}
				}else{
					total_time +=  dt;
				}
			}
		}
	}
			
	cout << "!wo gaps!The whole trajectory takes: " <<  total_time << " seconds." << endl;
	myfile << "% TransWORHP-Result" << endl;
	myfile << "% " << std::scientific << total_time << endl;
	double current_time = 0;
	double current_x 	= start_pos[0] ;
	double current_y 	= start_pos[1] ;
	double current_psi	= start_pos[2] ;
	double prev_x	 	= start_pos[0] ;
	double prev_y	 	= start_pos[1] ;
	double prev_psi		= start_pos[2] ;
	double prev_time	= 0;
	
	myfile 	<< setw(18) << std::scientific << current_time/total_time
			<< setw(18) << prev_x
			<< setw(18) << prev_y
			<< setw(18) << prev_psi
			<< setw(18) << 3.0
			<< setw(18) << 0.0
			<< setw(18) << 0.0
			<< setw(18) << 0.0
			<< setw(18) << 6.024201667e-01
			<< setw(18) << 0.0
			<< setw(18) << 0.0
			<< setw(18) << 0.0
			<< endl;
	
	for(int i = sol.size()-1; i>=0;--i){
		int prim_id = sol[i]->data.motion_type;
		//cout << "prim id: " << prim_id << endl;
		if(motion_prims.find(prim_id) != motion_prims.end()){
			for(int k=0; k<motion_prims[prim_id].size();++k){
				if(maneuvertimes.find(prim_id)!=maneuvertimes.end()){
					if(k>0){
						current_time = prev_time + motion_prims[prim_id][k][0];
						current_x	 = prev_x    + (motion_prims[prim_id][k][1])*cos(prev_psi)-(motion_prims[prim_id][k][2])*sin(prev_psi);
						current_y 	 = prev_y    + (motion_prims[prim_id][k][1])*sin(prev_psi)+(motion_prims[prim_id][k][2])*cos(prev_psi);
						current_psi  = prev_psi  +0 +motion_prims[prim_id][k][3];
					
						myfile 	<< setw(18) << std::scientific << current_time/total_time
						<< setw(18) << current_x 
						<< setw(18) << current_y 
						<< setw(18) << current_psi
						<< setw(18) << motion_prims[prim_id][k][4]
						<< setw(18) << motion_prims[prim_id][k][5]
						<< setw(18) << motion_prims[prim_id][k][6]
						<< setw(18) << motion_prims[prim_id][k][7]
						<< setw(18) << motion_prims[prim_id][k][8]
						<< setw(18) << motion_prims[prim_id][k][9]
						<< setw(18) << motion_prims[prim_id][k][10]
						<< setw(18) << 0.0
						<< endl;
						if(VERBOSE){
						//	cout << "current time: " << setw(8) << current_time <<  " || x: " << setw(8) << current_x << " || y: " << setw(8) << current_y << " || psi: " << setw(8) << current_psi*180/M_PI<< endl;
						}
					}
				}else{
					current_time = prev_time + dt;
					current_x	 = prev_x + motion_prims[prim_id][k][4]*cos(prev_psi)*dt;
					current_y	 = prev_y + motion_prims[prim_id][k][4]*sin(prev_psi)*dt;
					current_psi	 = prev_psi;
					myfile 	<< setw(18) << std::scientific << current_time/total_time
						<< setw(18) << current_x
						<< setw(18) << current_y
						<< setw(18) << current_psi 
						<< setw(18) << motion_prims[prim_id][k][4]
						<< setw(18) << motion_prims[prim_id][k][5]
						<< setw(18) << motion_prims[prim_id][k][6]
						<< setw(18) << motion_prims[prim_id][k][7]
						<< setw(18) << motion_prims[prim_id][k][8]
						<< setw(18) << motion_prims[prim_id][k][9]
						<< setw(18) << motion_prims[prim_id][k][10]
						<< setw(18) << 0.0
						<< endl;			
				//cout << "current time: " << setw(8) << current_time <<  " || x: " << setw(8) << current_x << " || y: " << setw(8) << current_y << " || psi: " << setw(8) << current_psi*180/M_PI<< endl;
				}
			
					
				if(k ==motion_prims[prim_id].size()-1){
						prev_x 	  = current_x;
						prev_y 	  = current_y;
						prev_psi  = current_psi;
						prev_time = current_time;
					}
			}
		}
	}
	myfile.close();
}

void PrimPathSearch::create_init_guess_with_holes(	vector< dsl::TypedEdge3d* > sol,vector<vector <double> > intersection_positions){
	struct dirent **namelist;
	string path_to_primitives = "../bin/motion_prims_3_cropped/";
	std::vector<std::string> file_names;
	int n = scandir(path_to_primitives.c_str(), &namelist, 0, versionsort);
	if (n < 0) perror(path_to_primitives.c_str());
	else{
		for(int i =0 ; i < n; ++i){
			file_names.push_back( namelist[i]->d_name );
			free(namelist[i]);
		}
		free(namelist);
	}
	// file to which the whole trajectory shall be written
	ofstream myfile;
	myfile.open ("example2.txt");

	//man that stores the maneuvertime to a given id 
	map<int,double> maneuvertimes;
	map<int,vector <vector<double>>> motion_prims;
	double t, x, y, psi, u, v, r, delta, rudder_rpm, d_delta, d_rudder_rpm;
					
	// variable to  storer the entire time of the trajectory
	double total_time = 0.0;
	// rückwärts durch Lösung weil vorwärts von a* gesucht wird, aber a* denkt es ist rückwärts
	for(int i = sol.size()-1; i>=0;--i){
		//cout << "prim id: " << sol[i]->data.motion_type << endl;
		string prim_type = (sol[i]->data.motion_type== 0)? "tp_01": "man_0"+to_string(sol[i]->data.motion_type);
		for(auto file: file_names){
			if (file.find(prim_type) == string::npos){
				continue;
			}
			else{
				
				// creating a map containing all occurring motion primitives for this trajectory
				if(motion_prims.find(sol[i]->data.motion_type)== motion_prims.end()){
					ifstream infile;
				
					//cout << path_to_primitives+file << endl;
					infile.open (path_to_primitives+file);	
					while (infile >> t >> x >> y >> psi >> u >> v >> r >> delta >> rudder_rpm >> d_delta >> d_rudder_rpm)
					{
						motion_prims[sol[i]->data.motion_type].push_back({t,x,y,psi,u,v,r,delta, rudder_rpm, d_delta, d_rudder_rpm});
					}
				}
				// chooping of the data-type
				const std::string ext(".txt");
				if ( file != ext && file.size() > ext.size() && file.substr(file.size() - ext.size()) == ".txt" )
				{
   					// if so then strip them off
					   file = file.substr(0, file.size() - ext.size());
				}
				// extract maneuver time from file name
				int underscore_counter =0;
				if(file.substr(0,3)=="man"){
					for(int j=0;j < file.length();++j){
						if(underscore_counter==7){
							total_time += stod(file.substr(j, file.size()));
							maneuvertimes[sol[i]->data.motion_type] =  stod(file.substr(j, file.size()));
							break;
						}
						if(file[j]=='_'){
							underscore_counter++;
						}
					}
				}else{
					total_time +=  dt;
				}
			}
		}
	}
			
	cout << "!w gaps!The whole trajectory takes: " <<  total_time << " seconds." << endl;
	myfile << "% TransWORHP-Result" << endl;
	myfile << "% " << std::scientific << total_time << endl;
	double current_time = 0;
	//alt
	/*double current_x 	= intersection_positions_w_gaps[sol.size()-1][0] ;
	double current_y 	= intersection_positions_w_gaps[sol.size()-1][1] ;
	double current_psi	= intersection_positions_w_gaps[sol.size()-1][2] ;
	double prev_x	 	= intersection_positions_w_gaps[sol.size()-1][0] ;
	double prev_y	 	= intersection_positions_w_gaps[sol.size()-1][1] ;
	double prev_psi		= intersection_positions_w_gaps[sol.size()-1][2] ;*/
	
	//neu, ansatz mit einer einzigen "create_init_guess"-funktion und intesection_positions als input
	double current_x 	  = intersection_positions[sol.size()-1][0] ;
	double current_y 	  = intersection_positions[sol.size()-1][1] ;
	double current_psi	  = intersection_positions[sol.size()-1][2] ;
	double current_energy = 0;
	double prev_x	 	  = intersection_positions[sol.size()-1][0] ;
	double prev_y	 	  = intersection_positions[sol.size()-1][1] ;
	double prev_psi		  = intersection_positions[sol.size()-1][2] ;
	double prev_energy    = 0;
	
	
	double prev_time	= 0;
	myfile 	<< setw(18) << std::scientific << current_time/total_time
			<< setw(18) << prev_x
			<< setw(18) << prev_y
			<< setw(18) << prev_psi
			<< setw(18) << 3.0
			<< setw(18) << 0.0
			<< setw(18) << 0.0
			<< setw(18) << 0.0
			<< setw(18) << 6.024201667e-01
			<< setw(18) << prev_energy
			<< setw(18) << 0.0
			<< setw(18) << 0.0
			//<< setw(18) << prev_energy
			<< endl;
	
	for(int i = sol.size()-1; i>=0;--i){
		int prim_id = sol[i]->data.motion_type;
		//cout << "prim id: " << prim_id << endl;
		if(motion_prims.find(prim_id) != motion_prims.end()){
			//alt
			/*
			prev_x   = intersection_positions_w_gaps[i][0];
			prev_y   = intersection_positions_w_gaps[i][1];
			prev_psi = intersection_positions_w_gaps[i][2];
			*/

			//neu, ansatz mit einer einzigen "create_init_guess"-funktion und intesection_positions als input
			prev_x   = intersection_positions[i][0];
			prev_y   = intersection_positions[i][1];
			prev_psi = intersection_positions[i][2];
			
			/*Display the difference between end of prim and middle of the cell, expected max diffrence for x,y,psi: d_nepsi_grid/sqrt(2)
			  currently x,y: 7/2*sqrt(2) psi pi/100*sqrt(2)
			
			cout << "###############################" << endl;
			cout << setw(10) << "x_diff: " << setw(8) << pow(prev_x-current_x,2) <<endl;
			cout << setw(10) <<"y_diff: " << setw(8) << pow(prev_y-current_y,2) <<endl;
			cout << setw(10) << "psi_diff: " << setw(8) << pow(prev_psi-current_psi,2) <<endl;
			cout << "###############################" << endl;
			*/
			for(int k=0; k<motion_prims[prim_id].size();++k){
				if(maneuvertimes.find(prim_id)!=maneuvertimes.end()){
				if(k>0){
						current_time   = prev_time + motion_prims[prim_id][k][0];
						current_x	   = prev_x    + (motion_prims[prim_id][k][1])*cos(prev_psi)-(motion_prims[prim_id][k][2])*sin(prev_psi);
						current_y 	   = prev_y    + (motion_prims[prim_id][k][1])*sin(prev_psi)+(motion_prims[prim_id][k][2])*cos(prev_psi);
						current_psi    = prev_psi  + motion_prims[prim_id][k][3];
						current_energy = prev_energy + (motion_prims[prim_id][k][0]-motion_prims[prim_id][k-1][0])/total_time*(motion_prims[prim_id][k][9]*motion_prims[prim_id][k][9]+motion_prims[prim_id][k][10]*motion_prims[prim_id][k][10]);

						myfile 	<< setw(18) << std::scientific << current_time/total_time
						<< setw(18) << current_x 
						<< setw(18) << current_y 
						<< setw(18) << current_psi
						<< setw(18) << motion_prims[prim_id][k][4]
						<< setw(18) << motion_prims[prim_id][k][5]
						<< setw(18) << motion_prims[prim_id][k][6]
						<< setw(18) << motion_prims[prim_id][k][7]
						<< setw(18) << motion_prims[prim_id][k][8]
						//<< setw(18) << current_energy
						//<< setw(18) << 0.0
						<< setw(18) << motion_prims[prim_id][k][9]
						<< setw(18) << motion_prims[prim_id][k][10]
						//<< setw(18) << current_energy
						<< setw(18) << 0.0
						<< endl;

						prev_energy = current_energy;
						//cout << "current time: " << setw(8) << current_time <<  " || x: " << setw(8) << current_x << " || y: " << setw(8) << current_y << " || psi: " << setw(8) << current_psi*180/M_PI<< endl;
					}
				}else{
					//cout << "ich bin ein trim" << endl;
					current_time = prev_time + dt;
					current_x	 = prev_x + motion_prims[prim_id][k][4]*cos(prev_psi)*dt;
					current_y	 = prev_y + motion_prims[prim_id][k][4]*sin(prev_psi)*dt;
					current_psi	 = prev_psi;
					myfile 	<< setw(18) << std::scientific << current_time/total_time
						<< setw(18) << current_x
						<< setw(18) << current_y
						<< setw(18) << current_psi 
						<< setw(18) << motion_prims[prim_id][k][4]
						<< setw(18) << motion_prims[prim_id][k][5]
						<< setw(18) << motion_prims[prim_id][k][6]
						<< setw(18) << motion_prims[prim_id][k][7]
						<< setw(18) << motion_prims[prim_id][k][8]
						//<< setw(18) << current_energy
						//<< setw(18) << 0.0
						<< setw(18) << motion_prims[prim_id][k][9]
						<< setw(18) << motion_prims[prim_id][k][10]
						//<< setw(18) << current_energy
						<< setw(18) << 0.0
						<< endl;			
					//cout << "current time: " << setw(8) << current_time <<  " || x: " << setw(8) << current_x << " || y: " << setw(8) << current_y << " || psi: " << setw(8) << current_psi*180/M_PI<< endl;
				}
			
				
				if(k ==motion_prims[prim_id].size()-1){
						/*prev_x 	  = current_x;
						prev_y 	  = current_y;
						prev_psi  = current_psi;*/
						prev_time = current_time;
					}
			}
		}
	}
	myfile.close();
}

void PrimPathSearch::create_init_guess_without_reparsing(	vector< dsl::TypedEdge3d* > sol,vector<vector <double> > intersection_positions,    std::vector< prim > prims){

	// file to which the whole trajectory shall be written
	ofstream myfile;
	myfile.open ("example3.txt");

	double total_time =0;
	double total_lgt =0;
	for(int i = sol.size()-1; i>=0;--i){
		int prim_id = sol[i]->data.motion_type;
		if(prim_id!=0){
			total_time += prims[prim_id].t[prims[prim_id].t.size()-1];
			total_lgt  += prims[prim_id].d[prims[prim_id].d.size()-1];
		}else{
			total_time +=dt;
			total_lgt  += prims[prim_id].u[prims[prim_id].u.size()-1]*dt;
		}
	}
	cout << "The whole trajectory takes: " <<  total_time << " seconds." << endl;
	cout << "The whole trajectory is: " <<  total_lgt << " meters long." << endl;
	myfile << "% TransWORHP-Result" << endl;
	myfile << "% " << std::scientific << total_time << endl;
	//neu, ansatz mit einer einzigen "create_init_guess"-funktion und intesection_positions als input
	double current_t = 0;
	int last_prim_sol = sol.size()-1; 
	myfile 	<< setw(18) << std::scientific << current_t/total_time
			<< setw(18) << intersection_positions[last_prim_sol][0]
			<< setw(18) << intersection_positions[last_prim_sol][1]
			<< setw(18) << intersection_positions[last_prim_sol][2]
			<< setw(18) << prims[sol[last_prim_sol]->data.motion_type].u[0]
			<< setw(18) << prims[sol[last_prim_sol]->data.motion_type].v[0]
			<< setw(18) << prims[sol[last_prim_sol]->data.motion_type].r[0]
			<< setw(18) << prims[sol[last_prim_sol]->data.motion_type].delta[0]
			<< setw(18) << prims[sol[last_prim_sol]->data.motion_type].n[0]
			<< setw(18) << prims[sol[last_prim_sol]->data.motion_type].energy[0]
			<< setw(18) << prims[sol[last_prim_sol]->data.motion_type].d_delta[0]
			<< setw(18) << prims[sol[last_prim_sol]->data.motion_type].d_n[0]
			<< endl;
	
	double current_x;
	double current_y;
	double current_psi;
	double current_energy;
	double prev_t;
	double prev_x;
	double prev_y;
	double prev_psi;
	double prev_energy = prims[sol[last_prim_sol]->data.motion_type].energy[0];
	for(int i = sol.size()-1; i>=0;--i){
		int prim_id = sol[i]->data.motion_type;
		/*
		cout << "prim id: " << prim_id << endl;
		cout << "intersection_positions:" << endl;
		cout << "current time: " << setw(8) << current_t <<  " || x: " << setw(8) << intersection_positions[i][0] << " || y: " << setw(8) << intersection_positions[i][1] << " || psi: " << setw(8) <<intersection_positions[i][2]*180/M_PI<< endl;	
		cout << "#################################################################" << endl;*/
		
		prev_x   = intersection_positions[i][0];
		prev_y   = intersection_positions[i][1];
		prev_psi = intersection_positions[i][2];


		/*Display the difference between end of prim and middle of the cell, expected max diffrence for x,y,psi: d_nepsi_grid/sqrt(2)
		  currently x,y: 7/2*sqrt(2) psi pi/100*sqrt(2)
		
		cout << "###############################" << endl;
		cout << setw(10) << "x_diff: " << setw(8) << pow(prev_x-current_x,2) <<endl;
		cout << setw(10) <<"y_diff: " << setw(8) << pow(prev_y-current_y,2) <<endl;
		cout << setw(10) << "psi_diff: " << setw(8) << pow(prev_psi-current_psi,2) <<endl;
		cout << "###############################" << endl;
		*/
		for(int k=0; k<prims[prim_id].x.size();++k){
			if(prim_id!=0){
				if(k>0){
					current_t      = prev_t      + prims[prim_id].t[k];
					current_x	   = prev_x      + (prims[prim_id].x[k])*cos(prev_psi)-(prims[prim_id].y[k])*sin(prev_psi);
					current_y 	   = prev_y      + (prims[prim_id].x[k])*sin(prev_psi)+(prims[prim_id].y[k])*cos(prev_psi);
					current_psi    = prev_psi    + prims[prim_id].psi[k];
					current_energy = prev_energy + (prims[prim_id].t[k]-prims[prim_id].t[k-1])*(pow(prims[prim_id].d_delta[k],2)+pow(prims[prim_id].d_n[k],2));
				
					myfile 	<< setw(18) << std::scientific << current_t/total_time
							<< setw(18) << current_x 
							<< setw(18) << current_y 
							<< setw(18) << current_psi
							<< setw(18) << prims[prim_id].u[k]
							<< setw(18) << prims[prim_id].v[k]
							<< setw(18) << prims[prim_id].r[k]
							<< setw(18) << prims[prim_id].delta[k]
							<< setw(18) << prims[prim_id].n[k]
							<< setw(18) << current_energy
							<< setw(18) << prims[prim_id].d_delta[k]
							<< setw(18) << prims[prim_id].d_n[k]
							<< endl;

					prev_energy = current_energy;
					//cout << "current time: " << setw(8) << current_t <<  " || x: " << setw(8) << current_x << " || y: " << setw(8) << current_y << " || psi: " << setw(8) << current_psi*180/M_PI<< endl;
				}
			}else{
				current_t    = prev_t + dt;
				current_x	 = prev_x + prims[prim_id].u[k]*cos(prev_psi)*dt;
				current_y	 = prev_y + prims[prim_id].u[k]*sin(prev_psi)*dt;
				current_psi	 = prev_psi;
				
				myfile 	<< setw(18) << std::scientific << current_t/total_time
						<< setw(18) << current_x 
						<< setw(18) << current_y 
						<< setw(18) << current_psi
						<< setw(18) << prims[prim_id].u[k]
						<< setw(18) << prims[prim_id].v[k]
						<< setw(18) << prims[prim_id].r[k]
						<< setw(18) << prims[prim_id].delta[k]
						<< setw(18) << prims[prim_id].n[k]
						<< setw(18) << prev_energy
						<< setw(18) << prims[prim_id].d_delta[k]
						<< setw(18) << prims[prim_id].d_n[k]
						<< endl;			
				//cout << "current time: " << setw(8) << current_t <<  " || x: " << setw(8) << current_x << " || y: " << setw(8) << current_y << " || psi: " << setw(8) << current_psi*180/M_PI<< endl;
			}
						
			if(k ==prims[prim_id].x.size()-1){
					prev_t = current_t;
			}
		}
	}
	myfile.close();
}

PrimGridPathSearch3d::PrimGridPathSearch3d (TypedGraph3d& graph_, TypedMap& map_, PathCost3d& cost_,
									PathGrid3d& grid_, double lat0_, double lon0_, double hgt0_, double d_xy_,
									Eigen::Vector3d d_nepsi_grid_, double d_occ_map_,int n_trims_)
	:Search< TypedCell3d, EdData3d > ( graph_, cost_ ),
	grid(grid_),
	map(map_),
	lat0(lat0_),
	lon0(lon0_),
	hgt0(hgt0_),
	d_xy(d_xy_),
	d_nepsi_grid(d_nepsi_grid_),
	no_heading_goal_cell(NULL),
	d_occ_map(d_occ_map_),
	save_sol_to_file(true),
	n_trims(n_trims_)
{
	vertexMap = new TypedVertex3d* [(grid.nc+1)];
	memset ( vertexMap, 0, (grid.nc+1) * sizeof ( TypedVertex3d* ) );
	string path_to_primitives = "../bin/motion_prims_5/";
	if (!SetPrimitives(path_to_primitives)) throw std::runtime_error ( "Error while executing function: SetPrimitives!" );
};

PrimGridPathSearch3d::~PrimGridPathSearch3d()
{
	for ( auto edge : graph.edges ) {
		delete edge.second;
	}
	for ( auto vertex : graph.vertices ) {
		delete vertex.second;
	}
	delete[] vertexMap;
	if(no_heading_goal_cell) delete no_heading_goal_cell;
}

template <typename T>
void PrimGridPathSearch3d::zero(T** buf, size_t size)
{
    	size_t my_start, my_size;

	    if (omp_in_parallel())
    	{
        	int id = omp_get_thread_num();
			int num = omp_get_num_threads();

        	my_start = (id*size)/num;
        	my_size = ((id+1)*size)/num - my_start;
			printf("my_start: %ld \n", my_start);
			printf("my_end: %ld\n", my_start+my_size);
			printf("my_size: %ld\n", my_size);
    	}
    	else
    	{
        	my_start = 0;
        	my_size = size;
    	}

	    memset(buf + my_start, 0, my_size);
}


bool PrimGridPathSearch3d::Expand ( TypedVertex3d& from, bool fwd )
{
	// if this is true then all vertices have already been expanded at
	// construction
	if ( fwd && from.succExpanded ) {
		return true;
	}

	if ( !fwd && from.predExpanded ) {
		return true;
	}

	// no-heading goal cell cannot be expanded
	if(from.data.id != grid.nc){
		// cell must exist
		TypedCell3d* from_cell = grid.Get ( from.data.id );
		assert ( from_cell );

		// there  are 3 possible motion primitive (for each neighbor-cell):
		//	0: straight
		//	1: 10 degrees left
		//	2: 10 degrees right

		/*std::vector<int> reachable_prims;
		if(fwd){
			reachable_prims = prim_links_fwd[(int)from.data.c[3]];
		}else{
			reachable_prims = prim_links_bwd[(int)from.data.c[3]];
		}*/
		for ( unsigned int motion_type=0; motion_type<nmbr_motion_primitives; motion_type++ ) {
			TypedCell3d *to_cell;
			EdData3d* eddata = new EdData3d();
			if ( !Flow( from_cell, to_cell, eddata, motion_type, fwd) ) {
				continue;
			}
			assert ( to_cell );

			int id = to_cell->id;
			// if this vertex doesn't exist, create it and add to graph
			if ( !vertexMap[id] ) {
				vertexMap[id] = new TypedVertex3d ( *to_cell );
				graph.AddVertex ( *vertexMap[id] );
			}
			TypedVertex3d* to = vertexMap[id];

			// if fwd and incoming edge from->to exists, then do not create a new one
			if ( from.Find ( *to, !fwd ) ) {
				continue;
			}

			TypedEdge3d* edge = fwd ? new TypedEdge3d ( *eddata, &from, to, eddata->costs ) :
			new TypedEdge3d ( *eddata, to, &from, eddata->costs );
			graph.AddEdge ( *edge );

			if(save_sol_to_file) visited_pos.push_back(to->data.c);
		}

	}

	if(Vertices()>1.0E6){
		std::cout << "Search exceeds search limit: abort" << std::endl;

		if(visited_pos.size()>0){
			std::string filename2 = "visited_positions.m";
			std::ofstream myfile2 (filename2);
			if (!myfile2.is_open()){
				std::cout << "Unable to open file: " << filename2 << std::endl;
				exit(1);
			}
			myfile2 << std::setprecision(16);
			myfile2 << "vis_pos=[";
			for(unsigned int i=0;i<visited_pos.size();i++){
				myfile2<<visited_pos[i][0]<<" "<<visited_pos[i][1]<<" "<<visited_pos[i][2]<<";\n";
			}
			myfile2 << "];";
			myfile2.close();
		}

		exit(1);
	}

	if ( fwd ) {
		from.succExpanded = true;
	} else {
		from.predExpanded = true;
	}
	return true;
}

bool PrimGridPathSearch3d::Flow ( TypedCell3d* from_cell, TypedCell3d*& to_cell, EdData3d* eddata,
							  int motion_type, bool fwd) {
	// meine Vektoren/parameter um excutePrim aufzurufen

	double x0_ned=from_cell->c[0];
	double y0_ned=from_cell->c[1];
	double psi0=from_cell->c[2];
	
	Eigen::Vector3d v0_nepsi=from_cell->c;

	Eigen::Vector3d vf_nepsi(0., 0., 0.);

	const Eigen::Vector3d v0(psi0, x0_ned, y0_ned);
	double dist=0.;
	if(!ExecutePrim<false, true>(v0, vf_nepsi, dist, motion_type, dt, fwd ) ){
		return false;	
	}
	
	// hier habe ich die Reihenfolge der eintraege geändert, die waren vertauscht
	const Eigen::Vector3d vf_nepsi_copy(vf_nepsi[1],vf_nepsi[2],vf_nepsi[0]);
	vf_nepsi = vf_nepsi_copy;
	
	to_cell = grid.Get ( vf_nepsi );
	if(!to_cell) return false;
	if( map.Get(Eigen::Vector2d(to_cell->c[0],to_cell->c[1])) ) return false;

	eddata->costs = dist;

	// special termination for no heading in goal state
	if( no_heading_goal_cell && to_cell->c[0]==no_heading_goal_cell->c[0]
			&& to_cell->c[1]==no_heading_goal_cell->c[1] ){
		to_cell=no_heading_goal_cell;
	}

	if(to_cell->data.path_costs > from_cell->data.path_costs+eddata->costs){
		to_cell->data.path_costs=from_cell->data.path_costs+eddata->costs;
		to_cell->data.nr_of_motions=from_cell->data.nr_of_motions+1;
	}

	eddata->motion_type = motion_type;

	if(fwd){
		eddata->from_cell_id = from_cell->id;
		eddata->to_cell_id = to_cell->id;
	}else{
		eddata->from_cell_id = to_cell->id;
		eddata->to_cell_id = from_cell->id;
	}
	return true;
}

bool PrimGridPathSearch3d::SetStart ( const Eigen::Vector3d& x, VertData3d data_ )
{
	if ( !grid.Valid ( x ) ) {
		std::cout << "[W] PrimGridPathSearch3d:SetStart: invalid x=" << x.transpose()
		<< std::endl;
		return false;
	}

	TypedCell3d* cell = grid.Get ( x );
	if ( !cell ) {
		std::cout << "[W] PrimGridPathSearch3d:SetStart: cell at=" << x.transpose() << " does not exist!"
		<< std::endl;
		return false;
	}

	int id = cell->id;
	assert ( id >= 0 && id < grid.nc );

	if( map.Get( Eigen::Vector2d(x[0], x[1]) )) 
	{ 
		std::cout << "[W] PrimGridPathSearch3d:SetStart: Starting Point: " << x.transpose() <<" is occupied!"<< std::endl;
		return false;
	}

	cell->data = data_;

	// if it's not added previously add it
	if ( !vertexMap[id] ) 
	{
		vertexMap[id] = new TypedVertex3d ( *cell );
		graph.AddVertex ( *vertexMap[id] );
	}

	//DSL serches backwards!
	Search< TypedCell3d, EdData3d >::SetGoal ( *vertexMap[id] );

	return true;
}

bool PrimGridPathSearch3d::SetGoal ( const Eigen::Vector3d& x, VertData3d data_ )
{
	
	if ( !grid.Valid ( x ) ) {
		std::cout << "[W] PrimGridPathSearch3d:SetGoal: invalid x=" << x.transpose()
		<< std::endl;
		return false;
	}

	TypedCell3d* cell;

	// if goal-heading is "nan": create special cell with no heading for termination!
	if( x[2] != x[2] ){

		Eigen::Vector3d x_aux(x);
		x_aux[2]=(grid.xlb[2]+grid.xub[2])/2.;
		cell = grid.Get ( x_aux );
		// x_aux[0]=cell->c[0];
		// x_aux[1]=cell->c[1];
		// x_aux[2]=x[2];
		// no_heading_goal_cell = new TypedCell3d ( grid.nc, x_aux );
		no_heading_goal_cell = new TypedCell3d ( grid.nc, cell->c );
		cell=no_heading_goal_cell;
	// else: proceed as normal
	}else{
		cell = grid.Get ( x );
		if ( !cell ) {
			std::cout << "[W] PrimGridPathSearch3d:SetGoal: cell at=" << x.transpose() <<
						 " does not exist!" << std::endl;
			return false;
		}
	}

	if( map.Get( Eigen::Vector2d(x[0], x[1])) ) {
		std::cout << "[W] PrimGridPathSearch3d:SetGoal: Goal Point: " << x.transpose() <<
		" is occupied!" << std::endl;
		return false;
	}

	if( map.Get( Eigen::Vector2d(cell->c[0], cell->c[1])) ) {
		std::cout << "[W] PrimGridPathSearch3d:SetGoal: Goal Points corresponding cell: " <<
					cell->c.transpose() << " is occupied!" << std::endl;
		return false;
	}
	
	int id = cell->id;
	assert ( id >= 0 && id < grid.nc+1 );

	cell->data = data_;

	// if it's not added previously add it
	if ( !vertexMap[id] ) {
		vertexMap[id] = new TypedVertex3d ( *cell );
		graph.AddVertex ( *vertexMap[id] );
	}

	//DSL serches backwards!
	Search< TypedCell3d, EdData3d >::SetStart ( *vertexMap[id] );

	return true;
}

    template<bool save_path_bool, bool crop_psi_bool>
bool PrimGridPathSearch3d::ExecutePrim( const Eigen::Vector3d& v0, Eigen::Vector3d& vf, double& d, const int prim_nr, const double dt, const bool fwd){
	if( prim_nr < n_trims ) return ExecuteTrim< save_path_bool, crop_psi_bool>( v0, vf, d, prim_nr, dt, fwd);
	else return ExecuteManeuver<save_path_bool, crop_psi_bool>( v0, vf, d, prim_nr, fwd);
	
}

double PrimGridPathSearch3d::CropPsi ( double psi )
{
	if( psi>=0 ) return fmod(psi,2*M_PI);
	else return 2*M_PI+fmod(psi,2*M_PI);
}

//!Vorsicht psi x y und nicht x y psi!
template<bool save_path_bool, bool crop_psi_bool>
bool PrimGridPathSearch3d::ExecuteManeuver( const Eigen::Vector3d& v0, Eigen::Vector3d& vf, double& d, const int man_nr, const bool fwd){

	//number of discretisation points of the maneuver
	int n = prims[man_nr].x.size()-1;
	
	//initital state of the maneuver
	double psi0=v0[0];
	double x0=v0[1];
	double y0=v0[2];
	
	// overall distance sailed during this maneuver
	d = prims[man_nr].d[n];
	
	// number of intermediate points that are check for validity 
	int m = static_cast<int>(round(d/d_xy));
	
	//fraction of the maneuver sailed to get from one intermediate point to the next
	double d_aux = d/((double) m);

	// NEU //
	double d_d_vec[m];
	double x_vec[m];
	double y_vec[m];
	double psi_vec[m];
	Eigen::Vector3d vf_vec[m]; 

	//fwd negiert, weil durch einen Trick die Suche vorwärts anstatt rückwärts ausgeführt wirdhttps://9gag.com/
	if(!fwd){
			
		int i =1;
		for (int j=1; j<m; j++) {
			for (; i<=n && ((double) j)*d_aux>prims[man_nr].d[i]; i++);
			double d_d = ( ((double) j)*d_aux - prims[man_nr].d[i-1])/(prims[man_nr].d[i]-prims[man_nr].d[i-1]);
			double x       = prims[man_nr].x[i-1]       + d_d * (prims[man_nr].x[i]-prims[man_nr].x[i-1]);
			double y       = prims[man_nr].y[i-1]       + d_d * (prims[man_nr].y[i]-prims[man_nr].y[i-1]);
			double psi     = prims[man_nr].psi[i-1]     + d_d * (prims[man_nr].psi[i]-prims[man_nr].psi[i-1]);
						
			// in case a solution was found the whole trajectory is executed again, but this time, headings of more than 2pi are feasible
			if(crop_psi_bool){
				
				vf[0] = CropPsi( psi0+psi );
			}
			else{
				vf[0] = psi0+psi;
			}
			//computation of the intermediate positions from this maneuver
			vf[1]=x0+x*cos(psi0)-y*sin(psi0);
			vf[2]=y0+x*sin(psi0)+y*cos(psi0);
			
			//check the spacial coordinates for validity
			if( map.Get(Eigen::Vector2d(vf[1],vf[2])) ){
				return false;
			}
			//like crop_psi_cool, after successful termination of a* the trajectory is executed again but this time stored
			if(save_path_bool){
				
				saved_path.push_back({0, vf[1], vf[2],vf[0]});
			}
		}
		//computation of the final positions from this maneuver
		double x       = prims[man_nr].x[n];
		double y       = prims[man_nr].y[n];
		double psi     = prims[man_nr].psi[n];
		if(crop_psi_bool){
				vf[0] = CropPsi( psi0+psi );
			}
			else{
				vf[0] = psi0+psi;
			}
		vf[1]=x0+x*cos(psi0)-y*sin(psi0);
		vf[2]=y0+x*sin(psi0)+y*cos(psi0);
		if( map.Get(Eigen::Vector2d(vf[1],vf[2])) ){
			return false;
		}
		if(save_path_bool){
				saved_path.push_back({0, vf[1],vf[2],vf[0]});
			}
	}else{
		//same procedure as in the case of not fwd, but this time everythis is executed backwards, in our case not important due to the missing final heading 
		double xf = prims[man_nr].x[n];
		double yf = prims[man_nr].y[n];
		double psif = prims[man_nr].psi[n];
		double psi1 = psi0 - psif;
		int i=n-1;
		for (int j=1; j<m; j++) {
			for (; i>=0 && (d-((double) j)*d_aux)<prims[man_nr].d[i]; i--);
			double d_d     = (prims[man_nr].d[i+1]      - (d-((double) j)*d_aux))/(prims[man_nr].d[i+1]-prims[man_nr].d[i]);
			double x       = prims[man_nr].x[i+1]       - d_d * (prims[man_nr].x[i+1]-prims[man_nr].x[i]) - xf;
			double y       = prims[man_nr].y[i+1]       - d_d * (prims[man_nr].y[i+1]-prims[man_nr].y[i]) - yf;
			double psi     = prims[man_nr].psi[i+1]     - d_d * (prims[man_nr].psi[i+1]-prims[man_nr].psi[i]) - psif;
			if(crop_psi_bool){
				vf[0] = CropPsi( psi0+psi );
			}
			else{
				vf[0] = psi0+psi;
			}
			vf[1]=x0+x*cos(psi1)-y*sin(psi1);
			vf[2]=y0+x*sin(psi1)+y*cos(psi1);
			if( map.Get(Eigen::Vector2d(vf[1],vf[2])) ){
				return false;
			}
			if(save_path_bool){
				saved_path.push_back({0, vf[1],vf[2],vf[0]});
			}
		}
		double x       = prims[man_nr].x[0] - xf;
		double y       = prims[man_nr].y[0] - yf;
		if(crop_psi_bool){
			vf[0]=CropPsi( psi0 - psif );
		}
		else{
			vf[0] = psi0 - psif; 
		}
		vf[1]=x0+x*cos(psi1)-y*sin(psi1);
		vf[2]=y0+x*sin(psi1)+y*cos(psi1);
		if( map.Get(Eigen::Vector2d(vf[1],vf[2])) ){
			return false;
		}
		if(save_path_bool){
				saved_path.push_back({0, vf[1],vf[2],vf[0]});
		}
	}
	return true;
}

template<bool save_path_bool, bool crop_psi_bool>
bool PrimGridPathSearch3d::ExecuteTrim( const Eigen::Vector3d& v0, Eigen::Vector3d& vf, double& d, const int trim_nr, const double dt, const bool fwd){
	
	double x0=v0[1];
	double y0=v0[2];
	double psi0=v0[0];
	double dt_aux;
	//fwd negiert, weil durch einen Trick die Suche vorwärts anstatt rückwärts ausgeführt wird
	if(!fwd){
		dt_aux=dt;
	}	
	else {
		dt_aux=-dt;
	}
	
	
	// distinction between no motion, straight trims and circular trims

	// straight (motion id 0)
	if(trim_nr==0){
		double u = prims[trim_nr].u[0];
		double dx = u*cos(psi0);
		double dy = u*sin(psi0);
		vf[0] = psi0;
		d = u*dt;
		int m = round(d/d_xy);
		for (int j=1; j<m; j++) {
			vf[1] = x0 + ((double)j)/((double)m)*dt_aux*dx;
			vf[2] = y0 + ((double)j)/((double)m)*dt_aux*dy;
			if( map.Get(Eigen::Vector2d(vf[1],vf[2])) ){
				return false;
			}
			if(save_path_bool){
				saved_path.push_back({0, vf[1],vf[2],vf[0],});
			}
		}
		vf[1] = x0 + dt_aux*dx;
		vf[2] = y0 + dt_aux*dy;
		if( map.Get(Eigen::Vector2d(vf[1],vf[2])) ){
			return false;
		}
		if(save_path_bool){
			saved_path.push_back({0, vf[1],vf[2],vf[0]});
		}
		
	}
	return true;
}

ArcPathSearch::ArcPathSearch(GeometricConstr& geom_constr, Eigen::Vector3d d_nepsi_grid,
							 Eigen::Vector3d start_llpsi, Eigen::Vector3d goal_llpsi){

	dsl::TypedGraph3d graph;

	double x_l_ned = geom_constr.occ_map->xlb[0];
	double x_u_ned = geom_constr.occ_map->xub[0];
	double y_l_ned = geom_constr.occ_map->xlb[1];
	double y_u_ned = geom_constr.occ_map->xub[1];
	double psi_l = 0.;
	double psi_u = 2*M_PI;
	double d_occ_map=geom_constr.occ_map_cell_dim_ned[0];
	if(d_occ_map>geom_constr.occ_map_cell_dim_ned[1]) d_occ_map=geom_constr.occ_map_cell_dim_ned[1];
	int n_x = floor( (x_u_ned-x_l_ned)/d_nepsi_grid[0] );
	int n_y = floor( (y_u_ned-y_l_ned)/d_nepsi_grid[1] );
	int n_psi = floor( (psi_u-psi_l)/d_nepsi_grid[2] );

	dsl::PathGrid3d grid( Eigen::Vector3d(x_l_ned, y_l_ned, psi_l),
	Eigen::Vector3d(x_u_ned, y_u_ned, psi_u),
	dsl::Vector3i(n_x, n_y, n_psi) );

	dsl::PathCost3d cost = dsl::PathCost3d();

	dsl::GridPathSearch3d search(graph, *(geom_constr.occ_map), cost, grid,
								 geom_constr.lat0, geom_constr.lon0, geom_constr.hgt0,
								 d_nepsi_grid, d_occ_map);

	dsl::VertData3d vert_data;
	vert_data.path_costs = std::numeric_limits<double>::infinity();
	vert_data.nr_of_motions = 0;
	

	double dummy=0.;
	Eigen::Vector3d start_nepsi(0.,0.,0.);
	Eigen::Vector3d goal_nepsi(0.,0.,0.);
	LLAtoNED(start_llpsi[0], start_llpsi[1], 0., geom_constr.lat0, geom_constr.lon0, geom_constr.hgt0,
			 &start_nepsi[0], &start_nepsi[1], &dummy);
	start_nepsi[2]=start_llpsi[2];
	LLAtoNED(goal_llpsi[0], goal_llpsi[1], 0., geom_constr.lat0, geom_constr.lon0, geom_constr.hgt0,
			 &goal_nepsi[0], &goal_nepsi[1], &dummy);
	goal_nepsi[2]=goal_llpsi[2];
	// Goal has to be set before Start
	if( !search.SetGoal(goal_nepsi, vert_data) ) exit(1);
	if( !search.SetStart(start_nepsi, vert_data) ) exit(1);



	vector< dsl::TypedEdge3d* > sol;

	struct timeval timer;
	long time = dsl::timer_us ( &timer );
	dsl::timer_start ( &timer );
	search.Plan(sol);
	time = dsl::timer_us ( &timer );
	printf ( "A* search plan time= %ld  us\n", time );
	cout << "Graph has " << search.Vertices() << " vertices and " << search.Edges() << " edges. " << endl;

	if(search.visited_pos.size()>0){
		string filename2 = "visited_positions.m";
		ofstream myfile2 (filename2);
		if (!myfile2.is_open()){
			cout << "Unable to open file: " << filename2 << endl;
			exit(1);
		}
		myfile2 << setprecision(16);
		myfile2 << "vis_pos=[";
		for(unsigned int i=0;i<search.visited_pos.size();i++){
			myfile2<<search.visited_pos[i][0]<<" "<<search.visited_pos[i][1]<<" "<<search.visited_pos[i][2]<<";\n";
		}
		myfile2 << "];";
		myfile2.close();
	}

	int n=sol.size();
	if(n>0){

		Eigen::Vector3d vf_nepsi(0.,0.,0.);
		double dist=0.;
		search.saved_path.resize(0);
		search.saved_path.push_back({0., sol[n-1]->to->data.c[0],sol[n-1]->to->data.c[1],sol[n-1]->to->data.c[2],
									 search.motion_values[sol[n-1]->data.motion_type][0],
									 search.motion_values[sol[n-1]->data.motion_type][1],
									 search.motion_values[sol[n-1]->data.motion_type][2],
									 search.motion_values[sol[n-1]->data.motion_type][3],
									 search.motion_values[sol[n-1]->data.motion_type][4]});
		for(int i=n-1; i>=0; i--){
			search.ExecuteMotion<true>(sol[i]->to->data.c, vf_nepsi,
								 	   sol[i]->data.motion_type, dist, false);
		}

		solution=search.saved_path;

		if(search.visited_pos.size()>0){
			string filename1 = "solution.m";
			ofstream myfile1 (filename1);
			if (!myfile1.is_open()){
				cout << "Unable to open file: " << filename1 << endl;
				exit(1);
			}
			myfile1 << "sol=[";
			for(int i=0; i<search.saved_path.size(); i++){
				myfile1<<search.saved_path[i][1]<<" "
					<<search.saved_path[i][2]<<" "
					<<search.saved_path[i][3]<<";\n";
			}
			myfile1 << "];";
			myfile1.close();
		}

	}else{
		cout << "Warning: No solution found!" << endl;
	}

};





GridPathSearch3d::GridPathSearch3d (TypedGraph3d& graph_, TypedMap& map_, PathCost3d& cost_,
									PathGrid3d& grid_, double lat0_, double lon0_, double hgt0_,
									Eigen::Vector3d d_nepsi_grid_, double d_occ_map_)
	:Search< TypedCell3d, EdData3d > ( graph_, cost_ ),
	grid(grid_),
	map(map_),
	lat0(lat0_),
	lon0(lon0_),
	hgt0(hgt0_),
	d_nepsi_grid(d_nepsi_grid_),
	no_heading_goal_cell(NULL),
	d_occ_map(d_occ_map_),
	save_sol_to_file(true)
{
	vertexMap = new TypedVertex3d* [grid.nc+1];
	memset ( vertexMap, 0, (grid.nc+1) * sizeof ( TypedVertex3d* ) );
}

GridPathSearch3d::~GridPathSearch3d()
{
	for ( auto edge : graph.edges ) {
		delete edge.second;
	}
	for ( auto vertex : graph.vertices ) {
		delete vertex.second;
	}
	delete[] vertexMap;
	if(no_heading_goal_cell) delete no_heading_goal_cell;
}

bool GridPathSearch3d::Expand ( TypedVertex3d& from, bool fwd )
{
	// if this is true then all vertices have already been expanded at
	// construction
	if ( fwd && from.succExpanded ) {
		return true;
	}

	if ( !fwd && from.predExpanded ) {
		return true;
	}

	// no-heading goal cell cannot be expanded
	if(from.data.id != grid.nc){
		// cell must exist
		TypedCell3d* from_cell = grid.Get ( from.data.id );
		assert ( from_cell );

		// there 8 possible motions (for each neighbor-cell):
		//	0: straight
		//	1: short-left
		//	2: long-left
		//	3: short-right
		//	4: long-right
		for ( unsigned int motion_type=0; motion_type<5; motion_type++ ) {

			TypedCell3d *to_cell;
			EdData3d* eddata = new EdData3d();
			if ( !Flow( from_cell, to_cell, eddata, motion_type, fwd) ) {
				continue;
			}
			assert ( to_cell );

			int id = to_cell->id;

			// if this vertex doesn't exist, create it and add to graph
			if ( !vertexMap[id] ) {
				vertexMap[id] = new TypedVertex3d ( *to_cell );
				graph.AddVertex ( *vertexMap[id] );
			}
			TypedVertex3d* to = vertexMap[id];

			// if fwd and incoming edge from->to exists, then do not create a new one
			if ( from.Find ( *to, !fwd ) ) {
				continue;
			}

			TypedEdge3d* edge = fwd ? new TypedEdge3d ( *eddata, &from, to, eddata->costs ) :
			new TypedEdge3d ( *eddata, to, &from, eddata->costs );
			graph.AddEdge ( *edge );

			if(save_sol_to_file) visited_pos.push_back(to->data.c);

		}

	}

	if(Vertices()>1.0E6){
		std::cout << "Search exceeds search limit: abort" << std::endl;

		if(visited_pos.size()>0){
			std::string filename2 = "visited_positions.m";
			std::ofstream myfile2 (filename2);
			if (!myfile2.is_open()){
				std::cout << "Unable to open file: " << filename2 << std::endl;
				exit(1);
			}
			myfile2 << std::setprecision(16);
			myfile2 << "vis_pos=[";
			for(unsigned int i=0;i<visited_pos.size();i++){
				myfile2<<visited_pos[i][0]<<" "<<visited_pos[i][1]<<" "<<visited_pos[i][2]<<";\n";
			}
			myfile2 << "];";
			myfile2.close();
		}

		exit(1);
	}

	if ( fwd ) {
		from.succExpanded = true;
	} else {
		from.predExpanded = true;
	}

	return true;
}

bool GridPathSearch3d::Flow ( TypedCell3d* from_cell, TypedCell3d*& to_cell, EdData3d* eddata,
							  int motion_type, bool fwd) {

	double dist=0.;
	double dummy=0.;

	double x0_ned=from_cell->c[0];
	double y0_ned=from_cell->c[1];
	double psi0=from_cell->c[2];
	Eigen::Vector3d v0_nepsi=from_cell->c;
	Eigen::Vector3d vf_nepsi(0., 0., 0.);

	if(!ExecuteMotion<false>(v0_nepsi, vf_nepsi, motion_type, dist, fwd) ) return false;

	to_cell = grid.Get ( vf_nepsi );
	if(!to_cell) return false;
	if( map.Get(Eigen::Vector2d(to_cell->c[0],to_cell->c[1])) ) return false;

	eddata->costs = dist;//cost.EdgeCost(*to_cell, *from_cell, dist);

	// special termination for no heading in goal state
	if( no_heading_goal_cell && to_cell->c[0]==no_heading_goal_cell->c[0]
			&& to_cell->c[1]==no_heading_goal_cell->c[1] ){
		to_cell=no_heading_goal_cell;
	}

	if(to_cell->data.path_costs > from_cell->data.path_costs+eddata->costs){
		to_cell->data.path_costs=from_cell->data.path_costs+eddata->costs;
		to_cell->data.nr_of_motions=from_cell->data.nr_of_motions+1;
	}

	eddata->motion_type = motion_type;

	if(fwd){
		eddata->from_cell_id = from_cell->id;
		eddata->to_cell_id = to_cell->id;
	}else{
		eddata->from_cell_id = to_cell->id;
		eddata->to_cell_id = from_cell->id;
	}

	return true;
}

template<bool save_path>
bool GridPathSearch3d::ExecuteMotion(Eigen::Vector3d v0_nepsi, Eigen::Vector3d& vf_nepsi,
									 int motion_type, double& dist, bool fwd){

	double psi0 = v0_nepsi[2];
	// choose motion via motion_type:
	//	0: straight
	if(motion_type==0){
		Eigen::Vector3d d_v(0., 0., 0.);
		d_v[0]= cos(psi0)*straight_mov_fac*d_nepsi_grid[0];
		d_v[1]= sin(psi0)*straight_mov_fac*d_nepsi_grid[1];
		d_v[2]= 0.;
		dist=sqrt(d_v[0]*d_v[0]+d_v[1]*d_v[1]);

		//dsl searches backwards!
		if(fwd) d_v=-1.*d_v;

		int m = ceil(dist/d_occ_map);
		double d_t = (dist/((double)m))/(sqrt(motion_values[motion_type][0]*motion_values[motion_type][0]
											  + motion_values[motion_type][1]*motion_values[motion_type][1]));
		for (int j=1; j<m; j++) {
			vf_nepsi=v0_nepsi+((double)j)/((double)m)*d_v;
			if( map.Get(Eigen::Vector2d(vf_nepsi[0],vf_nepsi[1])) ) return false;
			if( save_path ) saved_path.push_back({d_t, vf_nepsi[0],vf_nepsi[1],vf_nepsi[2],
												  motion_values[motion_type][0], motion_values[motion_type][1],
												  motion_values[motion_type][2], motion_values[motion_type][3],
												  motion_values[motion_type][4]});
		}

		vf_nepsi=v0_nepsi+d_v;
		if( map.Get(Eigen::Vector2d(vf_nepsi[0],vf_nepsi[1])) ) return false;
		if( save_path ) saved_path.push_back({d_t, vf_nepsi[0],vf_nepsi[1],vf_nepsi[2],
											  motion_values[motion_type][0], motion_values[motion_type][1],
											  motion_values[motion_type][2], motion_values[motion_type][3],
											  motion_values[motion_type][4]});

	//	1-4: arc
	}else if(motion_type>=1 && motion_type<=4){

		double radius = 0.;
		double angle = 0.;
		double beta = 0.;
		double gamma = 0.;

		//specify arc motion:
		switch ( motion_type ){
			//	1: short-left
			case 1:
				angle = -1.*angle_arc;
				beta = 1./2.*M_PI;
				gamma = 3./2.*M_PI;
				break;
			//	2: long-left
			case 2:
				angle = -2.*angle_arc;
				beta = 1./2.*M_PI;
				gamma = 3./2.*M_PI;
				break;
			//	3: short-right
			case 3:
				angle = angle_arc;
				beta = 3./2.*M_PI;
				gamma = 1./2.*M_PI;
				break;
			//	4: long right
			case 4:
				angle = 2.*angle_arc;
				beta = 3./2.*M_PI;
				gamma = 1./2.*M_PI;
				break;
		}

		dist = radius_arc*abs(angle);
		int m = ceil(dist/d_occ_map);
		double d_t = (dist/((double)m))/(sqrt(motion_values[motion_type][0]*motion_values[motion_type][0]
											  + motion_values[motion_type][1]*motion_values[motion_type][1]));
		
		// DSL searches backwards
		if(fwd) angle=-1.*angle;

		double xc_ned = v0_nepsi[0] + radius_arc*cos(psi0-beta);
		double yc_ned = v0_nepsi[1] + radius_arc*sin(psi0-beta);

		for (int j=1; j<m; j++) {
			vf_nepsi[0] = xc_ned + radius_arc*cos( ((double)j)/((double)m)*angle + psi0-gamma) ;
			vf_nepsi[1] = yc_ned + radius_arc*sin( ((double)j)/((double)m)*angle + psi0-gamma) ;
			vf_nepsi[2] = CropPsi( psi0 + ((double)j)/((double)m)*angle );

			if( map.Get(Eigen::Vector2d(vf_nepsi[0],vf_nepsi[1])) ) return false;
			if( save_path ) saved_path.push_back({d_t, vf_nepsi[0],vf_nepsi[1],vf_nepsi[2],
												  motion_values[motion_type][0], motion_values[motion_type][1],
												  motion_values[motion_type][2], motion_values[motion_type][3],
												  motion_values[motion_type][4]});
		}
		vf_nepsi[0] = xc_ned + radius_arc*cos( angle + psi0-gamma) ;
		vf_nepsi[1] = yc_ned + radius_arc*sin( angle + psi0-gamma) ;
		vf_nepsi[2] = CropPsi( psi0 + angle );

		if( map.Get(Eigen::Vector2d(vf_nepsi[0],vf_nepsi[1])) ) return false;
		if( save_path ) saved_path.push_back({d_t, vf_nepsi[0],vf_nepsi[1],vf_nepsi[2],
											  motion_values[motion_type][0], motion_values[motion_type][1],
											  motion_values[motion_type][2], motion_values[motion_type][3],
											  motion_values[motion_type][4]});

	//	>=5: not implemented
	}else{
		std::cout<<"Error inGridPathSearch3d::Flow: motion "<<motion_type<<
				   " not implemented"<<std::endl;
		assert(false);
	}
	return true;
}

bool GridPathSearch3d::SetStart ( const Eigen::Vector3d& x, VertData3d data_ )
{
	if ( !grid.Valid ( x ) ) {
		std::cout << "[W] GridPathSearch3d:SetStart: invalid x=" << x.transpose()
		<< std::endl;
		return false;
	}

	TypedCell3d* cell = grid.Get ( x );
	if ( !cell ) {
		std::cout << "[W] GridPathSearch3d:SetStart: cell at=" << x.transpose() << " does not exist!"
		<< std::endl;
		return false;
	}

	int id = cell->id;
	assert ( id >= 0 && id < grid.nc );

	if( map.Get( Eigen::Vector2d(x[0], x[1]) )) {
		std::cout << "[W] GridPathSearch3d:SetStart: Starting Point: " << x.transpose() <<
		" is occupied!"
		<< std::endl;
		return false;
	}

	cell->data = data_;

	// if it's not added previously add it
	if ( !vertexMap[id] ) {
		vertexMap[id] = new TypedVertex3d ( *cell );
		graph.AddVertex ( *vertexMap[id] );
	}

	//DSL serches backwards!
	Search< TypedCell3d, EdData3d >::SetGoal ( *vertexMap[id] );

	return true;
}

bool GridPathSearch3d::SetGoal ( const Eigen::Vector3d& x, VertData3d data_ )
{
	
	if ( !grid.Valid ( x ) ) {
		std::cout << "[W] GridPathSearch3d:SetGoal: invalid x=" << x.transpose()
		<< std::endl;
		return false;
	}

	TypedCell3d* cell;

	// if goal-heading is "nan": create special cell with no heading for termination!
	if( x[2] != x[2] ){

		Eigen::Vector3d x_aux(x);
		x_aux[2]=(grid.xlb[2]+grid.xub[2])/2.;
		cell = grid.Get ( x_aux );
		// x_aux[0]=cell->c[0];
		// x_aux[1]=cell->c[1];
		// x_aux[2]=x[2];
		// no_heading_goal_cell = new TypedCell3d ( grid.nc, x_aux );
		no_heading_goal_cell = new TypedCell3d ( grid.nc, cell->c );
		cell=no_heading_goal_cell;
	// else: proceed as normal
	}else{
		cell = grid.Get ( x );
		if ( !cell ) {
			std::cout << "[W] GridPathSearch3d:SetGoal: cell at=" << x.transpose() <<
						 " does not exist!" << std::endl;
			return false;
		}
	}

	if( map.Get( Eigen::Vector2d(x[0], x[1])) ) {
		std::cout << "[W] GridPathSearch3d:SetGoal: Goal Point: " << x.transpose() <<
		" is occupied!" << std::endl;
		return false;
	}

	if( map.Get( Eigen::Vector2d(cell->c[0], cell->c[1])) ) {
		std::cout << "[W] GridPathSearch3d:SetGoal: Goal Points corresponding cell: " <<
					cell->c.transpose() << " is occupied!" << std::endl;
		return false;
	}
	
	int id = cell->id;
	assert ( id >= 0 && id < grid.nc+1 );

	cell->data = data_;

	// if it's not added previously add it
	if ( !vertexMap[id] ) {
		vertexMap[id] = new TypedVertex3d ( *cell );
		graph.AddVertex ( *vertexMap[id] );
	}

	//DSL serches backwards!
	Search< TypedCell3d, EdData3d >::SetStart ( *vertexMap[id] );

	return true;
}

double GridPathSearch3d::CropPsi ( double psi )
{
	if( psi>=0 ) return fmod(psi,2*M_PI);
	else return 2*M_PI+fmod(psi,2*M_PI);
}

}
