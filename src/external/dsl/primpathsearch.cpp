#include "primpathsearch.h"

namespace dsl
{

using std::vector;
/**
 * this function is designed to create prims, whicha re defined in the prim header. Input of the function is a path to the
 * diretory containing all primitives. To avoid undefined behaviour one has to follow a simple syntax while naming the prim
 * itives. The syntax is "typeOfPrimitve_ID_describtion.txt" the description is optional, (%TODO connector snytax). If the
 * primitives were created with transWORHP one may use the "parse2primitive.py" script to obtain the needed name and format
 * of the primitves.
 * The type of the primitive is either "man" or "tp" this function counts the number of trims and maneuvers respectively.
 * Afterwards the textfiles are stored line by line into the prim struct. An information that is not provided by the textfiles
 * expicitly is the length of the prims, but that is calculated afterwards via the displacement during the primitive.
 */

bool NEW_PrimGridPathSearch4d::SetPrimitives ( std::string prim_path ) {

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

	
	for (int i=0; i<n_trims;i++) {
		prim pr;
		pr.is_trim=true;
		pr.n_t_steps=0;
		pr.costs=0.;
		std::ifstream f;
		std::string str;
		string namef = prim_path+file_names[i];
		cout << prim_path<<file_names[i] << endl;
		f.open ( prim_path+file_names[i].c_str(), std::ios::in );
		if ( !f ) {
			throw std::runtime_error ( "Error Reading File in ShipConnectivity::SetPrimitives; Wrong Filename/Path?" );
		}
		getline ( f,str );
		std::istringstream sin(str);
		double tmp;
		if( sin>>tmp ){
			pr.id=(int)tmp;
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
	try
  		{
			cout << "prim id im parser: " << pr.id << endl;
			trim_man_map_fwd[pr.id].push_back(pr.id);
			trim_man_map_bwd[pr.id].push_back(pr.id);
			pr.from = pr.id;
			pr.to   = pr.id;
	
    		throw 20;
  		}
  		catch (int e)
  		{
    		//cout << "An exception occurred. Exception Nr. " << e << '\n';
  		}
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
			pr.id=(int)tmp;
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

		try
  		{
			pr.from = std::stoi(file_names[i].substr(7,2));
			pr.to   = std::stoi(file_names[i].substr(10,2));

			trim_man_map_fwd[pr.from].push_back(pr.id);
			trim_man_map_bwd[pr.to].push_back(pr.id);
			throw 20;
  		}
  		catch (int e)
  		{
    		//cout << "An exception occurred. Exception Nr. " << e << '\n';
  		}
		

		pr.d.push_back(0);
		prims.push_back(pr);
	}
	for (int i=n_trims; i<n_prims;i++) {
		prims[i].d.resize(prims[i].x.size());
		for (unsigned int j=1; j<prims[i].d.size(); j++) {
			prims[i].d[j] = prims[i].d[j-1]
							 + sqrt( (prims[i].x[j]-prims[i].x[j-1])*(prims[i].x[j]-prims[i].x[j-1])
								+ (prims[i].y[j]-prims[i].y[j-1])*(prims[i].y[j]-prims[i].y[j-1]) );
		}
	}

	cout << "Printing for each trim the concatiable prims for backwards search:" << endl;

	for(std::map<int, vector<int> >::const_iterator it = trim_man_map_bwd.begin(); it != trim_man_map_bwd.end(); ++it)	
	{
		int n = it->second.size();
		std::cout << "trim id: " <<  it->first << " => ";
		for(int k=0; k< n; ++k){
		std::cout << it->second[k] << " ";
		}
		cout << endl;
	}

	return true;
}

/**
 * The class NEW_PrimpathSearch is used to initialize the a* search and  to handle (un-)successful termination of the search.
 * The needed parameters to start the search are, the upper and lower bounds of the grid (including the angular bounds). From
 * the information, and the dimension of angrid cell the number of grid cells can be obtained. Then the grid is initialized.
 * Then a object of the class NEW_PrimGridPathSearch4d is created, this will executes the a* search. Finally, start and end
 * position of the a* search are transformed from global spherical coordinates into local plain coordinates via the function
 * LLAtoNED. Then the a* search is performed.
 * Once the a* search is done the visited positions are stored into a file, this information is only required purposes for 
 * visualisation. More important is the next step, the primitives that forms the solution are executed a second time, this
 * is done to compute the displacements that are caused by the discretisation of the geometry. This is an inportant step, because
 * otherwise it would come to a discrepancy between the a* search solution and the generated by the motion primitives (because a
 * concatination of them contains no displacements)
 * 
 * Inputs:
 * 		geom_constr is a bnary map containing all blocked and all free cells
 * 		d_nepsi_grid contains the dimensions of a grid 
 * 		start_llpsi starting state, position in lat lon
 * 		goal_llpsi final state, position in lat lon
 * 		path2prims path to the primitive used in the a* search
 * 		a_star_result_file path to a file where data about the a* search is stored to.
 */ 
NEW_PrimPathSearch::NEW_PrimPathSearch(GeometricConstr& geom_constr, Eigen::Vector4d d_nepsi_grid, Eigen::Vector3d start_llpsi, Eigen::Vector3d goal_llpsi, string path2prims, string a_star_result_file)
{
	dsl::TypedGraph4d graph;

	double x_l_ned = geom_constr.occ_map->xlb[0]; //lower bound of the habour in x direction
	double x_u_ned = geom_constr.occ_map->xub[0]; //upper bound of the habour in x direction
	double y_l_ned = geom_constr.occ_map->xlb[1]; //lower bound of the habour in y direction
	double y_u_ned = geom_constr.occ_map->xub[1]; //upper bound of the habour in y direction
	// eps = 1.E-9 wegen grid.valid check bei winkeln
	double psi_l = 0.-1.E-9;					  //lower bound of the habour in psi direction
	double psi_u = 2*M_PI+1.E-9;				  //upper bound of the habour in psi direction

	// neu 4D prim bounds
	double trim_l = -0.5;					  //lower bound of the habour in prim direction
	double trim_u =  1.5;				 	  //upper bound of the habour in prim direction
	
	//setting the dimension of the occupancy gridmap to the min of its dimensions
	double d_occ_map=geom_constr.occ_map_cell_dim_ned[0];
	if(d_occ_map>geom_constr.occ_map_cell_dim_ned[1]) d_occ_map=geom_constr.occ_map_cell_dim_ned[1];

	//computing the number of grid cells for each dimension
	int n_x    = floor( (x_u_ned-x_l_ned)/d_nepsi_grid[0] );
	int n_y    = floor( (y_u_ned-y_l_ned)/d_nepsi_grid[1] );
	int n_psi  = floor( (psi_u-psi_l)/d_nepsi_grid[2] );
	int n_trim = round( (trim_u-trim_l)/d_nepsi_grid[3] );


	// dimension of the occcupancy grid map
	double d_xy = sqrt(geom_constr.occ_map->cs[0]*geom_constr.occ_map->cs[0]+geom_constr.occ_map->cs[1]*geom_constr.occ_map->cs[1]);
	
	//creating a grid with the above defined variables
	dsl::PathGrid4d grid( Eigen::Vector4d(x_l_ned, y_l_ned, psi_l, trim_l),
						  Eigen::Vector4d(x_u_ned, y_u_ned, psi_u, trim_u),
						  dsl::Vector4i(n_x, n_y, n_psi, n_trim) );

	//initialising the cost (%TODO recherchieren wofür)
	dsl::PathCost4d cost = dsl::PathCost4d();

	//initializing the a* search
	NEW_PrimGridPathSearch4d search(graph, *(geom_constr.occ_map), cost, grid,
								 geom_constr.lat0, geom_constr.lon0, geom_constr.hgt0, d_xy,
								 d_nepsi_grid, d_occ_map, 1, path2prims, a_star_result_file);
	dsl::VertData4d vert_data;
	vert_data.path_costs = std::numeric_limits<double>::infinity();
	vert_data.nr_of_motions = 0;
	
	// conversion from LLA to NED
	double dummy=0.;
	// neu 4D start und end grid cell dimension increased
	Eigen::Vector4d start_nepsi(0.,0.,0.,1.);
	Eigen::Vector4d goal_nepsi(0.,0.,0.,0.);
	LLAtoNED(start_llpsi[0], start_llpsi[1], 0., geom_constr.lat0, geom_constr.lon0, geom_constr.hgt0,
			 &start_nepsi[0], &start_nepsi[1], &dummy);
	start_nepsi[2]=start_llpsi[2];
	LLAtoNED(goal_llpsi[0], goal_llpsi[1], 0., geom_constr.lat0, geom_constr.lon0, geom_constr.hgt0,
			 &goal_nepsi[0], &goal_nepsi[1], &dummy);
	goal_nepsi[2]=goal_llpsi[2];
	// Goal has to be set before Start
	if( !search.SetGoal(goal_nepsi, vert_data) ) exit(1);
	if( !search.SetStart(start_nepsi, vert_data) ) exit(1);


	// vector containing the ID's of all primitives forming a solution
	vector< dsl::TypedEdge4d* > sol;

	// launching the a* search
	struct timeval timer;
	long time = dsl::timer_us ( &timer );
	dsl::timer_start ( &timer );
	search.Plan(sol);
	
	
	time = dsl::timer_us ( &timer );
	printf ( "A* search plan time= %ld  us\n", time );
	cout << "Graph has " << search.Vertices() << " vertices and " << search.Edges() << " edges. " << endl;
	
	// computing the total cost of the solution
	double cost_sum=0.;
	for(int i=0; i<sol.size(); i++){
		cost_sum+=sol[i]->data.costs;
	}
	cout << "Final costs=" << cost_sum <<endl;
	
	// storing all the visited position to a matlab-file, needed for visualisation in matlab
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

	/** in case of a successfull termination of the a* search the primitives forming a solution are executed again and all the endpoints
	 *  of prims are stored to a vector, this is done with and without repositioning. Repositioning guarantees, that all geometric constraints are satisfied.
	 */ 



	int n=sol.size();
	if(n>0)
	{
		intersection_positions_w_gaps.resize(n);
		intersection_positions_wo_gaps.resize(n);
		vector<Eigen::Vector3d> vec;
		double dist=0.;
		search.saved_path.resize(0);
		search.saved_path.push_back({0., sol[n-1]->to->data.c[0],sol[n-1]->to->data.c[1],sol[n-1]->to->data.c[2]});
		// neu 4D appended zeros to vertor
		Eigen::Vector4d v0_nepsi_w_gaps(sol[n-1]->to->data.c[2],sol[n-1]->to->data.c[0],sol[n-1]->to->data.c[1],0.);
		Eigen::Vector4d vf_nepsi_w_gaps(0.,0.,0.,0.);
		Eigen::Vector4d v0_nepsi_wo_gaps(start_nepsi[2],start_nepsi[0],start_nepsi[1],0.);
		Eigen::Vector4d vf_nepsi_wo_gaps(0.,0.,0.,0.);
		intersection_positions_w_gaps[n-1] = {v0_nepsi_w_gaps[1],v0_nepsi_w_gaps[2],v0_nepsi_w_gaps[0],v0_nepsi_w_gaps[3]};
		intersection_positions_wo_gaps[n-1] = {v0_nepsi_wo_gaps[1],v0_nepsi_wo_gaps[2],v0_nepsi_wo_gaps[0],v0_nepsi_wo_gaps[3]};
		
		for(int i=n-1; i>=1; i--)
		{

			cout << sol[i]->data.motion_type << endl;
			//without gaps
			search.ExecutePrim<false,false>( v0_nepsi_wo_gaps, vf_nepsi_wo_gaps,  dist, sol[i]->data.motion_type, dt, false);
			v0_nepsi_wo_gaps=vf_nepsi_wo_gaps;
			//todo hier fehlt evtl eine dimension
			intersection_positions_wo_gaps[i-1] = {v0_nepsi_wo_gaps[1],v0_nepsi_wo_gaps[2],v0_nepsi_wo_gaps[0]};

			//Display of the stating and the end points of each prim
			if(VERBOSE){
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
			}

			//with gaps
			search.ExecutePrim<true,false>( v0_nepsi_w_gaps, vf_nepsi_w_gaps, dist, sol[i]->data.motion_type, dt, false);
			double psi_uncropped = vf_nepsi_w_gaps[0];
			double psi_cropped = search.CropPsi(psi_uncropped);
			// neu 4D null im vierten Eintrag geadded

			// in den letzten eintrag muss nepsi w gap 3 dafür setcell trim von flow in exeprim
			TypedCell4d* vf_cell = grid.Get ( Eigen::Vector4d(vf_nepsi_w_gaps[1],vf_nepsi_w_gaps[2],psi_cropped, vf_nepsi_w_gaps[3]) );
			if(!vf_cell)
			{
				 cout<<"Fatal Error: Cell doesn't exist"<<endl; exit(1);
			}
			// neu 4d added dimension
			v0_nepsi_w_gaps=Eigen::Vector4d(vf_cell->c[2]+(psi_uncropped-psi_cropped),vf_cell->c[0],vf_cell->c[1],vf_cell->c[3]);
			intersection_positions_w_gaps[i-1] = {v0_nepsi_w_gaps[1],v0_nepsi_w_gaps[2],v0_nepsi_w_gaps[0]};
			
			 //Display of the stating and the end points of each prim
			if(VERBOSE){
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
			}
		}

		solution=search.saved_path;

		if(search.visited_pos.size()>0)
		{
			// using the computed intermediate positions to create the initial guess
			create_init_guess_without_reparsing(sol,intersection_positions_w_gaps, search.prims);
			
			//parsing hte solution path to a matlab-file for visualisation
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

	// parsing data about the a* search into the file passed to this object
	ofstream output;
	output.open(a_star_result_file, ios::app);
	output  << setw(8) << total_length  << " " << setw(3) << sol.size() << " " << setw(8)<< time << " " << setw(7) << search.Vertices() << " " << setw(7) << search.Edges()<< " " << endl;
	output.close();
	

	}else{
		cout << "Warning: No solution found!" << endl;
		ofstream output;
		output.open(a_star_result_file, ios::app);
		// parsing data into the file with the right format, but identifying this search to be failed
		output  << setw(8) << 0  << " " << setw(3) << 0 << " " << setw(8)<< time << " " << setw(7) << 0 << " " << setw(7) << 0 << " " << endl;
		output.close();
		// exit code used to tell test_main that the a* search failed
		exit(1);
	}

};


/*  Old version of create_init guess, not used anymore*/
void NEW_PrimPathSearch::create_init_guess(	vector< dsl::TypedEdge4d* > sol){
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
							cout << "current time: " << setw(8) << current_time <<  " || x: " << setw(8) << current_x << " || y: " << setw(8) << current_y << " || psi: " << setw(8) << current_psi*180/M_PI<< endl;
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


/* second version of create init guess, also outdated*/
void NEW_PrimPathSearch::create_init_guess_with_holes(	vector< dsl::TypedEdge4d* > sol,vector<vector <double> > intersection_positions){
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

/**
 * this function generates an initial guess based on the a* search. The inputs for this functions are:
 * 		 sol a vector containing the primitive ID's
 * 		 intersection_positions vector containing the new starting positions of a prim 
 * 		 prims vector containg all the prim's (can be removed eventually)
 * 
 * First the total length and total time of the primitive are computed, the total time is needed, since transWORHP requieres an input
 * scaled to zero to one.  The final length of the prim is needed for the computation of a suitable amont of discretisation points.
 * After that, the file is created, fist the header is generated, then iteratively the prims are executed in the correct order. The final
 * position is stored and a new prim is shifted and rotated by the final position of the previous prim. Then, the file is written.
 */
void NEW_PrimPathSearch::create_init_guess_without_reparsing(	vector< dsl::TypedEdge4d* > sol,vector<vector <double> > intersection_positions,    std::vector< prim > prims){

	// file to which the whole trajectory shall be written
	ofstream myfile;
	myfile.open ("example3.txt");

	double total_time =0;
	total_length = 0;
	for(int i = sol.size()-1; i>=0;--i){
		int prim_id = sol[i]->data.motion_type;
		if(prim_id!=0){
			total_time += prims[prim_id].t[prims[prim_id].t.size()-1];
			total_length += prims[prim_id].d[prims[prim_id].d.size()-1];
		}else{
			total_time +=dt;
			total_length += dt*sqrt(pow(prims[prim_id].u[prims[prim_id].u.size()-1],2)+pow(prims[prim_id].v[prims[prim_id].v.size()-1],2));
		}
	}
	//cout << "The whole trajectory takes: " <<  total_time << " seconds." << endl;
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
		
		if(VERBOSE)
		{
			cout << "prim id: " << prim_id << endl;
			cout << "intersection_positions:" << endl;
			cout << "current time: " << setw(8) << current_t <<  " || x: " << setw(8) << intersection_positions[i][0] << " || y: " << setw(8) << intersection_positions[i][1] << " || psi: " << setw(8) <<intersection_positions[i][2]*180/M_PI<< endl;	
			cout << "#################################################################" << endl;
		}
		
		prev_x   = intersection_positions[i][0];
		prev_y   = intersection_positions[i][1];
		prev_psi = intersection_positions[i][2];


		/*Display the difference between end of prim and middle of the cell, expected max diffrence for x,y,psi: d_nepsi_grid/sqrt(2)
		  currently x,y: 7/2*sqrt(2) psi pi/100*sqrt(2)*/
		if(VERBOSE){
			cout << "###############################" << endl;
			cout << setw(10) << "x_diff: "   << setw(8) << pow(prev_x-current_x,2)     << endl;
			cout << setw(10) << "y_diff: "   << setw(8) << pow(prev_y-current_y,2)     << endl;
			cout << setw(10) << "psi_diff: " << setw(8) << pow(prev_psi-current_psi,2) << endl;
			cout << "###############################" << endl;
		}
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

NEW_PrimGridPathSearch4d::NEW_PrimGridPathSearch4d (TypedGraph4d& graph_, TypedMap& map_, PathCost4d& cost_,
									PathGrid4d& grid_, double lat0_, double lon0_, double hgt0_, double d_xy_,
									Eigen::Vector4d d_nepsi_grid_, double d_occ_map_,int n_trims_, string path2prims, string a_star_result_file)
	:Search< TypedCell4d, EdData4d > ( graph_, cost_ ),
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
	n_trims(n_trims_),
	a_star_result_file(a_star_result_file)
{
	vertexMap = new TypedVertex4d* [grid.nc+1];
	std::cout << "memset size: " << (grid.nc+1) * sizeof ( TypedVertex4d* ) << std::endl;
	memset ( vertexMap, 0, (grid.nc+1) * sizeof ( TypedVertex4d* ) );
	if(path2prims.back() != '/'){
		path2prims+= '/';
	}
	if (!SetPrimitives(path2prims)) throw std::runtime_error ( "Error while executing function: SetPrimitives!" );
};

NEW_PrimGridPathSearch4d::~NEW_PrimGridPathSearch4d()
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

bool NEW_PrimGridPathSearch4d::Expand ( TypedVertex4d& from, bool fwd )
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
		TypedCell4d* from_cell = grid.Get ( from.data.id );
		assert ( from_cell );

		std::vector<int> reachable_prims;
		
		//o2c sucht vorwärts weil kein ziel heading, daher ist zielcell die start Zelle für den A* algorithmus, dieser sucht rückwärts
		if(!fwd){
			reachable_prims = trim_man_map_fwd[from_cell->c[3]];
		}else{
			reachable_prims = trim_man_map_bwd[from_cell->c[3]];
		}

		/*
		std::cout << "reachable prims" << std::endl;
		for(const auto & prim : reachable_prims)
		{
			std::cout << prim << " ";
		}
		std::cout << std::endl;
		*/
		//cout << from_cell->motion_type << endl;
		//cout << "n_prims "<< n_prims << endl;
		// TODO diesen for-loop durch TRIM spezifischen loop ersetzen
		//for ( unsigned int motion_type=0; motion_type<n_prims; motion_type++ )
		for(std::vector<int>::iterator prim_it = reachable_prims.begin(); prim_it != reachable_prims.end(); ++prim_it) 
		 {
			//cout << "try prim: " << *prim_it << endl;
			TypedCell4d *to_cell;

			// ich glaube hier ist ein memory leak
			EdData4d* eddata = new EdData4d();

			if ( !Flow( from_cell, to_cell, eddata, *prim_it, fwd) ) {
				// delete eddata;  // missing here?
				continue;
			}
			

			assert ( to_cell );

			double lat_print;
			double lon_print;
			double h_print;


			int id = to_cell->id;
			//std::cout << "prim: " << *prim_it << std::endl;
			//NEDtoLLA(to_cell->c[0],to_cell->c[1],0,54.17057475, 12.10074142,0,&lat_print,&lon_print,&h_print);
			//std::cout << lat_print << "," << lon_print << std::endl;
			//std::cout << to_cell->c[0] << " " << to_cell->c[1] << " " << to_cell->c[2] << std::endl;
			// if this vertex doesn't exist, create it and add to graph
			if ( !vertexMap[id] ) {
				vertexMap[id] = new TypedVertex4d ( *to_cell );
				graph.AddVertex ( *vertexMap[id] );
			}
			TypedVertex4d* to = vertexMap[id];
			// if fwd and incoming edge from->to exists, then do not create a new one
			if ( from.Find ( *to, !fwd ) ) {
				continue;
			}
			TypedEdge4d* edge = fwd ? new TypedEdge4d ( *eddata, &from, to, eddata->costs ) :
			new TypedEdge4d ( *eddata, to, &from, eddata->costs );
			graph.AddEdge ( *edge );

			if(save_sol_to_file) visited_pos.push_back(to->data.c);
		}
	// do just one step (for testing)
	only_steps++;
	}

	if(Vertices()>1.0E6){
		std::cout << "Search exceeds search limit: abort" << std::endl;
		ofstream output;
		output.open(a_star_result_file, ios::app);
		output  << setw(8) << 0  << " " << setw(3) << 0 << " " << setw(8)<< 0 << " " << setw(7) << 0 << " " << setw(7) << 0 << " " << endl;
		output.close();
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

bool NEW_PrimGridPathSearch4d::Flow ( TypedCell4d* from_cell, TypedCell4d*& to_cell, EdData4d* eddata,
							  int motion_type, bool fwd) {
	// meine Vektoren/parameter um excutePrim aufzurufen
	if (only_steps == INT64_MAX) exit(1);
	double x0_ned = from_cell->c[0];
	double y0_ned = from_cell->c[1];
	double psi0   = from_cell->c[2];
	// neu 4D cell entry for prim if
	double prim0  = from_cell->c[3];
	

	/***   TODO: USElESS?   ***/
	Eigen::Vector4d v0_nepsi=from_cell->c;

	// cneu 4D dimension increasesed and zero added
	Eigen::Vector4d vf_nepsi(0., 0., 0.,0.);

	const Eigen::Vector4d v0(psi0, x0_ned, y0_ned, prim0);
	double dist=0.;
	if(!ExecutePrim<false, true>(v0, vf_nepsi, dist, motion_type, dt, fwd )){
		return false;	
	}


	// TODO liste die maneuver (und trims) auf to und from (trim) mappt

	
	//TODO vf_nepsi[3] da muss die "ZIEL"TRIM id rein
	//cout << "executing prim: " << motion_type << endl;
	//cout << "trim id of cell: " << vf_nepsi[3] << endl;	
	
	// hier habe ich die Reihenfolge der eintraege geändert, die waren vertauscht
	// neu 4d und dimension vergrößert
	//cout << vf_nepsi[1] << " " << vf_nepsi[2] << " "<< vf_nepsi[0] << " " << vf_nepsi[3] << endl;
	
	/***   TODO: inplace change of order   ***/
	const Eigen::Vector4d vf_nepsi_copy(vf_nepsi[1],vf_nepsi[2],vf_nepsi[0],vf_nepsi[3]);
	vf_nepsi = vf_nepsi_copy;
	
	to_cell = grid.Get ( vf_nepsi );
	if(!to_cell) return false;
	/***   TODO: USELESS check(checked in executeprim too)   ***/
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
	// TODO motion_type => prim_nr
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

bool NEW_PrimGridPathSearch4d::SetStart ( const Eigen::Vector4d& x, VertData4d data_ )
{
	cout << "start cell: " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << endl;
	if ( !grid.Valid ( x ) ) {
		std::cout << "[W] PrimGridPathSearch4d:SetStart: invalid x=" << x.transpose()
		<< std::endl;
		return false;
	}

	TypedCell4d* cell = grid.Get ( x );
	if ( !cell ) {
		std::cout << "[W] PrimGridPathSearch4d:SetStart: cell at=" << x.transpose() << " does not exist!"
		<< std::endl;
		return false;
	}

	int id = cell->id;
	assert ( id >= 0 && id < grid.nc );

	if( map.Get( Eigen::Vector2d(x[0], x[1]) )) 
	{ 
		std::cout << "[W] PrimGridPathSearch4d:SetStart: Starting Point: " << x.transpose() <<" is occupied!"<< std::endl;
		return false;
	}

	cell->data = data_;

	// if it's not added previously add it
	if ( !vertexMap[id] ) 
	{
		vertexMap[id] = new TypedVertex4d ( *cell );
		graph.AddVertex ( *vertexMap[id] );
	}

	//DSL serches backwards!
	Search< TypedCell4d, EdData4d >::SetGoal ( *vertexMap[id] );

	return true;
}

bool NEW_PrimGridPathSearch4d::SetGoal ( const Eigen::Vector4d& x, VertData4d data_ )
{
		cout << "goal prim: " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << endl;	
	if ( !grid.Valid ( x ) ) {
		std::cout << "[W] PrimGridPathSearch4d:SetGoal: invalid x=" << x.transpose()
		<< std::endl;
		return false;
	}

	TypedCell4d* cell;

	// if goal-heading is "nan": create special cell with no heading for termination!
	if( x[2] != x[2] ){

		Eigen::Vector4d x_aux(x);
		x_aux[2]=(grid.xlb[2]+grid.xub[2])/2.;
		x_aux[3]=(grid.xlb[3]+grid.xub[3])/2.;

		cell = grid.Get ( x_aux );
		//cout << "cell 3: " << cell->c[3] << endl;
		// x_aux[0]=cell->c[0];
		// x_aux[1]=cell->c[1];
		// x_aux[2]=x[2];
		// no_heading_goal_cell = new TypedCell3d ( grid.nc, x_aux );
		no_heading_goal_cell = new TypedCell4d ( grid.nc, cell->c );
		cell=no_heading_goal_cell;
		//cout << "cell 3: " << cell->c[3] << endl;
	// else: proceed as normal
	}else{
		cell = grid.Get ( x );
		if ( !cell ) {
			std::cout << "[W] PrimGridPathSearch4d:SetGoal: cell at=" << x.transpose() <<
						 " does not exist!" << std::endl;
			return false;
		}
	}

	if( map.Get( Eigen::Vector2d(x[0], x[1])) ) {
		std::cout << "[W] PrimGridPathSearch4d:SetGoal: Goal Point: " << x.transpose() <<
		" is occupied!" << std::endl;
		return false;
	}

	if( map.Get( Eigen::Vector2d(cell->c[0], cell->c[1])) ) {
		std::cout << "[W] PrimGridPathSearch4d:SetGoal: Goal Points corresponding cell: " <<
					cell->c.transpose() << " is occupied!" << std::endl;
		return false;
	}
	
	int id = cell->id;
	assert ( id >= 0 && id < grid.nc+1 );

	cell->data = data_;

	// if it's not added previously add it
	if ( !vertexMap[id] ) {
		vertexMap[id] = new TypedVertex4d ( *cell );
		graph.AddVertex ( *vertexMap[id] );
	}

	//DSL serches backwards!
	Search< TypedCell4d, EdData4d >::SetStart ( *vertexMap[id] );

	return true;
}

double NEW_PrimGridPathSearch4d::CropPsi ( double psi )
{
	if( psi>=0 ) return fmod(psi,2*M_PI);
	else return 2*M_PI+fmod(psi,2*M_PI);
}

template<bool save_path_bool, bool crop_psi_bool>
bool NEW_PrimGridPathSearch4d::ExecutePrim( const Eigen::Vector4d& v0, Eigen::Vector4d& vf, double& d, const int prim_nr, const double dt, const bool fwd){
	if( prim_nr < n_trims ) return ExecuteTrim< save_path_bool, crop_psi_bool>( v0, vf, d, prim_nr, dt, fwd);
	else return ExecuteManeuver<save_path_bool, crop_psi_bool>( v0, vf, d, prim_nr, fwd);
	
}
//!Vorsicht psi x y und nicht x y psi!
template<bool save_path_bool, bool crop_psi_bool>
bool NEW_PrimGridPathSearch4d::ExecuteManeuver( const Eigen::Vector4d& v0, Eigen::Vector4d& vf, double& d, const int man_nr, const bool fwd){

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

	//fwd negiert, weil durch einen Trick die Suche vorwärts anstatt rückwärts ausgeführt wird
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
			//comtodoputation of the intermediate positions from this maneuver
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

		// hier haben wir wieder wegen der gedrehten suchrichtung from und to vertauscht
		for(int i =0; i < prims.size(); ++i) 
		{
			if(man_nr == prims[i].id)
			{
				vf[3] = (double) prims[i].to;
			}
		}
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
				// hier haben wir wieder wegen der gedrehten suchrichtung from und to vertauscht
		for(int i =0; i < prims.size(); ++i) 
		{
			if(man_nr == prims[i].id)
			{
				vf[3] = (double) prims[i].from;
			}
		}
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
bool NEW_PrimGridPathSearch4d::ExecuteTrim( const Eigen::Vector4d& v0, Eigen::Vector4d& vf, double& d, const int trim_nr, const double dt, const bool fwd)
{	
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
				saved_path.push_back({0, vf[1],vf[2],vf[0]});
			}
		}
		vf[1] = x0 + dt_aux*dx;
		vf[2] = y0 + dt_aux*dy;
				// hier haben wir wieder wegen der gedrehten suchrichtung from und to vertauscht
		for(int i =0; i < prims.size(); ++i) 
		{
			if(trim_nr == prims[i].id)
			{
				vf[3] = fwd ? (double) prims[i].from  : (double) prims[i].to;
			}
		}
		if( map.Get(Eigen::Vector2d(vf[1],vf[2])) ){
			return false;
		}
		if(save_path_bool){
			saved_path.push_back({0, vf[1],vf[2],vf[0]});
		}
		
	
	return true;
}
}