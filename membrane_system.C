#include "membrane_system.H" 

/*
int FormFunctionGradient(TAO_APPLICATION taoapp,Vec X,double *f, Vec G,void *ptr)
{
  AppCtx *user = (AppCtx *) ptr;  
  int    info;

	PetscInt nn=user->n/2;
	double ff=0;

  PetscScalar *x,*g;

  // Get pointers to vector data 
  
  info = VecGetArray(X,&x); CHKERRQ(info);
  info = VecGetArray(G,&g); CHKERRQ(info);

// user->sys->pull_solution();
	//std::cout<<"Someone calls me: \n";
	//VecView(X,	PETSC_VIEWER_STDOUT_WORLD );

	//THIS IS THE MAIN GRADIENT CALCULATION FUNCTION
	user->sys->function_gradient(&x, &ff, &g);


	
  // Restore vectors 
  info = VecRestoreArray(X,&x); CHKERRQ(info);
  info = VecRestoreArray(G,&g); CHKERRQ(info);
  *f=ff;
  
  // user->sys->pull_solution();

	//FIX
  info = PetscLogFlops(nn*15); CHKERRQ(info);
  return 0;
}

//we dont use hessian, it is empty
int FormHessian(TAO_APPLICATION taoapp,Vec X,Mat *HH, Mat *Hpre, MatStructure *flag,void *ptr)
{return 0;}

*/

MembraneSystem::MembraneSystem(int n_memb){
	this->n_membranes = n_memb;
	membranes = new Membrane[this->n_membranes];
	
	wallX=10.0;
	wallY=10.0;
	kinetic_energy = 0.;
	potential_energy = 0.;

	
	Estretch = 0.0;
	Ebending = 0.0;
	Ebend4 = 0.0;
	Ecoop = 0.0;
	Elandau = 0.0;
	Etotal = 0.0;
	
	
	tao_iter = 0; tao_iter0=0; vtk_iter=1;
}

void MembraneSystem::read_parameters(Reader* params){
	
	stretching_on = params->get_bool("stretching_on");
	bending_on = params->get_bool("bending_on");
	Srigidity1 = params->get_double("stretching_rigidity1");
	Srigidity2 = params->get_double("stretching_rigidity2");
	Srigidity3 = params->get_double("stretching_rigidity3");
	Srigidity4 = params->get_double("stretching_rigidity4");
	Srigidity5 = params->get_double("stretching_rigidity5");
	Brigidity1 = params->get_double("bending_rigidity1");
	Brigidity2 = params->get_double("bending_rigidity2");
	Brigidity3 = params->get_double("bending_rigidity3");
	Brigidity4 = params->get_double("bending_rigidity4");
	Brigidity5 = params->get_double("bending_rigidity5");
	d_bond = params->get_double("bond_length");
	d_bond1 = params->get_double("bond_length1");
	d_bond2 = params->get_double("bond_length2");
	alpha_1 = params->get_double("angle_1");
	alpha_2 = params->get_double("angle_2");
	alpha_3 = params->get_double("angle_3");
	alpha_4 = params->get_double("angle_4");
	alpha_5 = params->get_double("angle_5");
	mesh_file = params->get_char("mesh_file");
	mesh_type = params->get_int("mesh_type");
	has_seam = params->get_bool("has_seam");
	
	pA = Brigidity2*params->get_double("pA");
	pB = Brigidity2*params->get_double("pB");
	pC = Brigidity2*params->get_double("pC");
	limit_angle = params->get_double("barrier");
	harmonic = params->get_bool("harmonic");
		
	coop_coef = params->get_double("coop_coef");
	
	noise_strength = params->get_double("noise_strength");
	
	nodal_on = params->get_bool("nodal_on");
	mu_nodal = params->get_double("nodal_strength");
	max_force_strength = params->get_double("force_strength");
	boundary0 = params->get_double("boundary0");
	boundary1 = params->get_double("boundary1");
	boundary2 = params->get_double("boundary2");
	
	is_clamped = params->get_bool("is_clamped");
	// if not clamped set boundary1 to negative
	if(!is_clamped) boundary1 = -1.0;
	
	printeach = params->get_int("print_each");
	
	update_coef = 1.0/double(params->get_int("number_of_steps"));
	force_strength = max_force_strength*update_coef;
	force_is_gradual = params->get_bool("force_is_gradual");
	
	moment_on = params->get_bool("moment_on");
	moment_is_gradual = params->get_bool("moment_is_gradual");
	moment_boundary_up = params->get_double("moment_boundary_up");
	moment_boundary_down = params->get_double("moment_boundary_down");
	max_moment_strength = params->get_double("moment_strength");
	moment_strength = max_moment_strength*update_coef;
	
	container_on = params->get_bool("container_on");
	mu_cont = params->get_double("container_rigidity");
	cradius = params->get_double("container_radius");

	spin_on = params->get_bool("spin_on");
	deltaE = params->get_double("deltaE");
	JJ = params->get_double("JJ");
	
	gamma = params->get_double("gamma");
	mass = params->get_double("mass");
	temperature = params->get_double("temperature");
	boltzmann = params->get_double("boltzmann");
	diffusion = params->get_double("diffusion");
	dt = params->get_double("delta_t");
	
	//2 k_B T \gamma \Delta t
	noise_coef = sqrt(2.0*boltzmann*temperature*gamma*mass/dt);
	
	//diffusion_coef = boltzmann*temperature/gamma;
	diffusion_coef = diffusion/(boltzmann*temperature);
	brownian_coef = sqrt(2.0*diffusion*dt);
	
	//Initialize the random number generator
	r1 = new Ran(params->get_int("random_seed"));
	
	//Initialize the directory name
	this->directory = params->get_directory();
	
	//initilaize measurement file
	std::stringstream measure_name;
	measure_name<<directory<<"/measure.dat";	
	measurements.open(measure_name.str().c_str(),std::ios::out);
	measurements<<"time"<<"\t"<<"kinetic\t"<<"\t"<<"stretching"<<"\t"<<"bending\t"<<"bending4\t"<<"\t"<<"external\t"
	<<"ecoop\t"<<"total"<<std::endl;	
	

	std::cout<<"bending: "<<bending_on<<"\n";
	std::cout<<"str: "<<stretching_on<<"\n";
	std::cout<<"Srig: "<<Srigidity1<<"\n";
	std::cout<<"Brig: "<<Brigidity1<<"\n";
	std::cout<<"d_bond: "<<d_bond<<"\n";
	std::cout<<"a_1: "<<alpha_1<<"\n";
	std::cout<<"a_2: "<<alpha_2<<"\n";
	std::cout<<"a_3: "<<alpha_3<<"\n";
	std::cout<<"mass: "<<mass<<"\n";
	std::cout<<"gamma: "<<gamma<<"\n";
	std::cout<<"boltzmann: "<<boltzmann<<"\n";
	std::cout<<"temperature: "<<temperature<<"\n";
	std::cout<<"delta t: "<<dt<<"\n";
	std::cout<<"noise coef: "<<noise_coef<<"\n";
	std::cout<<"diffusion coef: "<<diffusion_coef<<"\n";
	std::cout<<"brownian coef: "<<brownian_coef<<"\n";
	std::cout<<std::endl;

}

//Thanks goes to libmesh off reader
void MembraneSystem::read_off(const char* off_file){
	std::ifstream offstream;
	offstream.open(off_file);
	if(offstream.is_open()) std::cout<<"mesh file is opened\n";
	else{ 
		std::cout<<"Unable to open mesh file\n";
		exit(0);
	}
	std::string label; // Read the first string.  It should say "OFF"
	offstream >> label;
	
	unsigned int nn, ne, nf;
	offstream >> nn >> nf >> ne; //number of nodes, faces, edges
	sample.resize(nn, nf);

	double x=0., y=0., z=0.;
  	// Read the nodes
	for (unsigned int n=0; n<nn; n++){
		offstream >> x >> y >> z;
		RealVectorValue pos(x,y,z);
		sample.get_node(n).get_xyz() = pos; 
		sample.get_node(n).get_xyz0() =pos; 
		sample.get_node(n).set_mass(this->mass);
    }
	
	unsigned int dummy, n0, n1, n2;
	// Read the triangles
	for (unsigned int e=0; e<nf; e++){
      // The number of nodes in the object
      offstream >> dummy;
      offstream >> n0 >> n1 >> n2;
      sample.get_triangle(e).set(n0,n1,n2);
    }
}


//Thanks goes to libmesh off reader
void MembraneSystem::read_rectangular_off(const char* off_file){
	std::ifstream offstream;
	offstream.open(off_file);
	if(offstream.is_open()) std::cout<<"mesh file is opened\n";
	else{ 
		std::cout<<"Unable to open mesh file\n";
		exit(0);
	}
	std::string label; // Read the first string.  It should say "OFF"
	offstream >> label;
	
	unsigned int nn, ne, nf;
	offstream >> nn >> nf >> ne; //number of nodes, faces, edges
	sample.resize(nn, nf);

	double x=0., y=0., z=0.;
  	// Read the nodes
	for (unsigned int n=0; n<nn; n++){
		offstream >> x >> y >> z;
		RealVectorValue pos(x,y,z);
		sample.get_node(n).get_xyz() = pos; 
		sample.get_node(n).get_xyz0() =pos; 
		sample.get_node(n).set_mass(this->mass);
    }
	
	unsigned int dummy, n0, n1, n2, n3;
	// Read the triangles
	for (unsigned int e=0; e<nf; e++){
      // The number of nodes in the object
      offstream >> dummy;
      offstream >> n0 >> n1 >> n2 >> n3;
      sample.get_rectangle(e).set(e,n0,n1,n2,n3);
    }
}


void MembraneSystem::write_off(){

	std::stringstream output_mesh;
	output_mesh<<directory;
	output_mesh<<"/result.off";
	std::ofstream omesh(output_mesh.str().c_str(),std::ios::out);
	
	omesh<<std::setprecision(10);
	omesh<<"OFF\n";
	omesh<<n_nodes<<" "<<n_triangles<<" 0"<<"\n";
	
	RealVectorValue pos;
	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over different membranes
		Membrane& m1 = membranes[mit];
		int nsize = m1.get_n_nodes();
		
		for (int nit = 0; nit<nsize; nit++){ //iterate over nodes
			Node& n1 = m1.get_node(nit);
			pos = n1.get_xyz();
			omesh<<pos(0)<<" "<<pos(1)<<" "<<pos(2)<<"\n";
		}	
		
		int esize = m1.get_n_triangles();
		for (int eit = 0; eit<esize; eit++){
			Triangle& e1 = m1.get_triangle(eit);
			omesh<<"3 "<<e1.get(0)<<" "<<e1.get(1)<<" "<<e1.get(2)<<"\n";
		}
	}
		omesh.close();
}

void MembraneSystem::input_solution(const char* in_file){
	//get the first membrane
	Membrane& m1 = membranes[0];
	std::ifstream instream;
	instream.open(in_file);
	if(instream.is_open()) std::cout<<"input file is opened\n";
	else{ 
		std::cout<<"Unable to open input file\n";
		exit(0);
	}
	std::string label; 
	
	// Read unused header
	char st[256];
	int nn = 0; //number of nodes
	do{
		instream.getline (st, 256);
		label = st;
		if(label.find("n_nodes") ==1){
			int n0 = label.find_first_of(" ");
			int n1 = label.size();
			nn = atoi(label.substr(n0+1,n1).c_str());
		}
	}while(st[0] == '#');

	kinetic_energy=0.0;
	double x=0., y=0., z=0., vx=0., vy=0., vz=0., sp = 0.;
  	// Read the nodes
	for (int n=0; n<nn; n++){
		instream >> x >> y >> z >> vx >> vy >> vz >> sp;
		RealVectorValue pos(x,y,z);
		RealVectorValue vel(vx,vy,vz);
		m1.get_node(n).get_xyz() += pos; 
		m1.get_node(n).get_velocity() = vel; 
		m1.get_node(n).spin = sp; 
		
		kinetic_energy += 0.5*m1.get_node(n).mass*vel*vel; //i update kinetic energy here for the 0th output
		//sample.get_node(n).get_xyz0() =pos; 
    }
	
}


void MembraneSystem::output_solution(){

	std::stringstream output;
	output<<directory;
	output<<"/solution.mdf";
	std::ofstream omesh(output.str().c_str(),std::ios::out);
	
	omesh<<std::setprecision(10);
	omesh<<"#Solution File\n";
	omesh<<"#Mesh File: "<<this->mesh_file<<"\n";
	omesh<<"#n_nodes "<<n_nodes<<"\n";
	omesh<<"# data order: u v w velx vely velz spin\n";
	omesh<<"\n"; //write an empty line to ease reading
	RealVectorValue pos;
	RealVectorValue vel;
	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over different membranes
		Membrane& m1 = membranes[mit];
		int nsize = m1.get_n_nodes();
		
		for (int nit = 0; nit<nsize; nit++){ //iterate over nodes
			Node& n1 = m1.get_node(nit);
			//output just the displacement with respect to reference mesh
			pos = n1.get_xyz()-n1.get_xyz0();
			vel = n1.get_velocity();
			omesh<<pos(0)<<" "<<pos(1)<<" "<<pos(2)<<" "<<vel(0)<<" "<<vel(1)<<" "<<vel(2)<<" "<<n1.spin<<"\n";
		}	
		
		//int esize = m1.get_n_faces();
		//for (int eit = 0; eit<esize; eit++){
			//Element& e1 = m1.get_face(eit);
			//omesh<<"3 "<<e1.get(0)<<" "<<e1.get(1)<<" "<<e1.get(2)<<"\n";
		//}
	}
		omesh.close();
}

void MembraneSystem::vtk_writer(int step){
	
	double ener=0.0, tener=0.0, sp=0.0;
	
	int dum=step;
	int n_dec=0;
	//find out how many digit is in time
	while(dum!=0){
		dum /=10; 
		n_dec++;
		}
	
	std::stringstream file_name;
	std::stringstream file_name2;

	//FIX: if data folder does not exist, it does nothing 
	//FIXED in parameter reader class
	
	file_name<<directory<<"/results-";
	file_name2<<directory<<"/results-";
	std::string vtu_source = "results-";
	//add enough 0's to the file name
	if (n_dec == 0) n_dec++; //otherwise after adding t=0 to file name, we have an extra digit
	while(n_dec != 8){
		file_name<<"0";
		file_name2<<"0";
		vtu_source +="0";
		n_dec++;
		}
	file_name<<step;
	file_name2<<step;
	vtu_source += to_string(step);
	file_name<<".vtu";
	file_name2<<".pvtu";
	vtu_source += ".vtu";
	//std::cout<<file_name.str()<<std::endl;

	//get the string, then make it char
	std::ofstream results(file_name.str().c_str(),std::ios::out);
	std::ofstream presults(file_name2.str().c_str(),std::ios::out);
	
	//write the header
	results<<"<?xml version=\"1.0\"?>\n";
	results<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
	results<<"\t<UnstructuredGrid>\n";

	//write the number of points
	results<<"\t<Piece NumberOfPoints=\""<<n_nodes<<"\" NumberOfCells=\""<<n_triangles<<"\">\n";
	
	//write some point data
	results<<"\t\t<PointData>\n";
		
		//ENERGY
		results<< "<DataArray type=\"Float64\" Name=\""<<"energy"<<"\" format=\"ascii\" RangeMin=\""<<-10000<<"\" RangeMax=\""<<10000<<"\">\n";
		for(int mem_iter = 0; mem_iter<n_membranes; mem_iter++){
			int msize = membranes[mem_iter].get_n_nodes();
			for(int part_iter = 0; part_iter<msize; part_iter++){
				ener = membranes[mem_iter].get_node(part_iter).energy;
			results<<"\t\t\t"<<ener<<std::endl;					
			}
		}		
		results<<"</DataArray>\n";
		
		//TOTAL ENERGY
		results<< "<DataArray type=\"Float64\" Name=\""<<"total_energy"<<"\" format=\"ascii\" RangeMin=\""<<-10000<<"\" RangeMax=\""<<10000<<"\">\n";
		for(int mem_iter = 0; mem_iter<n_membranes; mem_iter++){
			int msize = membranes[mem_iter].get_n_nodes();
			for(int part_iter = 0; part_iter<msize; part_iter++){
				tener = membranes[mem_iter].get_node(part_iter).total_energy;
			results<<"\t\t\t"<<tener<<std::endl;					
			}
		}		
		results<<"</DataArray>\n";

		//SPIN
		results<< "<DataArray type=\"Float64\" Name=\""<<"spin"<<"\" format=\"ascii\" RangeMin=\""<<-10<<"\" RangeMax=\""<<10<<"\">\n";
		for(int mem_iter = 0; mem_iter<n_membranes; mem_iter++){
			int msize = membranes[mem_iter].get_n_nodes();
			for(int part_iter = 0; part_iter<msize; part_iter++){
				sp = membranes[mem_iter].get_node(part_iter).spin;
			results<<"\t\t\t"<<sp<<std::endl;					
			}
		}		
		results<<"</DataArray>\n";
		
	results<<"\t\t</PointData>\n";
	
	//Cell data
	results<<"\t\t<CellData>\n"; 
	results<<"\t\t</CellData>\n";	
	
	
	
	results<<"\t\t<Points>\n";
	
	//Position data are coming
	results<<"\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"-100.0\" RangeMax=\"100.0\">\n";
	for(int mem_iter = 0; mem_iter<n_membranes; mem_iter++){
		int msize = membranes[mem_iter].get_n_nodes();
		for(int part_iter = 0; part_iter<msize; part_iter++){
			RealVectorValue p1 = membranes[mem_iter].get_node(part_iter).get_xyz();
			results<<"\t\t\t"<<p1(0)<<" "<<p1(1)<<" "<<p1(2)<<std::endl;
			//std::cout<<"vtk: "<<p1<<std::endl;					
		}
	}
	results<<"\t\t\t</DataArray>\n";
	results<<"\t\t</Points>\n";
	
	//Now Cell data are coming, for regular particles nothing interesting
	//FUTUR: may be interesting for bounded particles
	results<<"\t\t<Cells>\n";
	//connectivity
	results<<"\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\""<<n_triangles<<"\">\n";
	
	int sofar=0;
	for(int mem_iter = 0 ; mem_iter < n_membranes; ++mem_iter){
		int esize = membranes[mem_iter].get_n_triangles();
		for(int e_iter=0; e_iter < esize; e_iter++){
			Triangle el = membranes[mem_iter].get_triangle(e_iter);
			results<<"\t\t\t"<<el.get(0)+sofar<<" "<<el.get(1)+sofar<<" "<<el.get(2)+sofar<<std::endl;
		}
		sofar += membranes[mem_iter].get_n_nodes(); //shift by the sum of the n_nodes of prev membranes
	}
    results<<"\t\t\t</DataArray>\n";
    
    //offsets
    results<<"\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" RangeMin=\"1\" RangeMax=\""<<n_nodes<<"\">\n";
	for(int iter =0 ; iter < n_triangles; ++iter){ 
		results<<"\t\t\t"<<3*(iter+1)<<std::endl;
	}  
	results<<"\t\t\t</DataArray>\n";
	//types: 1 for regular particles
	// 5 for triangles
	results<<"\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\" RangeMin=\"1\" RangeMax=\"1\">\n";
	for(int iter =0 ; iter < n_triangles; ++iter){ 
		results<<"\t\t\t"<<5<<std::endl;
	}  
    results<<"\t\t\t</DataArray>\n";
	results<<"\t\t</Cells>\n";

//end of file
    results<<"\t</Piece>\n";
	results<<"\t</UnstructuredGrid>\n";
	results<<"</VTKFile>";
	results.close();
	

//WRITE THE PVTU file
		presults<<"<?xml version=\"1.0\"?>\n";
		presults<<"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
		presults<<"<PUnstructuredGrid GhostLevel=\"0\">\n";
		presults<<"<PPointData>\n";
		presults<<"<PDataArray type=\"Float64\" Name=\"energy\"/>\n";
		presults<<"<PDataArray type=\"Float64\" Name=\"total_energy\"/>\n";
		presults<<"<PDataArray type=\"Float64\" Name=\"spin\"/>\n";
		presults<<"</PPointData>\n";
		presults<<"<PPoints>\n";
        presults<<"<PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n";
		presults<<"</PPoints>\n";

		
		//std::cout<<vtu_source<<std::endl;

		presults<<"<Piece Source=\""<<vtu_source<<"\"/>\n";
		presults<<"</PUnstructuredGrid>\n";
		presults<<"</VTKFile>";	
		presults.close();
	
}


void MembraneSystem::vtk_rect_writer(int step){
	
	double ener=0.0, tener=0.0, sp=0.0;
	
	int dum=step;
	int n_dec=0;
	//find out how many digit is in time
	while(dum!=0){
		dum /=10; 
		n_dec++;
		}
	
	std::stringstream file_name;
	std::stringstream file_name2;

	//FIX: if data folder does not exist, it does nothing 
	//FIXED in parameter reader class
	
	file_name<<directory<<"/results-";
	file_name2<<directory<<"/results-";
	std::string vtu_source = "results-";
	//add enough 0's to the file name
	if (n_dec == 0) n_dec++; //otherwise after adding t=0 to file name, we have an extra digit
	while(n_dec != 8){
		file_name<<"0";
		file_name2<<"0";
		vtu_source +="0";
		n_dec++;
		}
	file_name<<step;
	file_name2<<step;
	vtu_source += to_string(step);
	file_name<<".vtu";
	file_name2<<".pvtu";
	vtu_source += ".vtu";
	//std::cout<<file_name.str()<<std::endl;

	//get the string, then make it char
	std::ofstream results(file_name.str().c_str(),std::ios::out);
	std::ofstream presults(file_name2.str().c_str(),std::ios::out);
	
	//write the header
	results<<"<?xml version=\"1.0\"?>\n";
	results<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
	results<<"\t<UnstructuredGrid>\n";

	//write the number of points
	results<<"\t<Piece NumberOfPoints=\""<<n_nodes<<"\" NumberOfCells=\""<<n_rectangles<<"\">\n";
	
	//write some point data
	results<<"\t\t<PointData>\n";
		
		//ENERGY
		results<< "<DataArray type=\"Float64\" Name=\""<<"energy"<<"\" format=\"ascii\" RangeMin=\""<<-10000<<"\" RangeMax=\""<<10000<<"\">\n";
		for(int mem_iter = 0; mem_iter<n_membranes; mem_iter++){
			int msize = membranes[mem_iter].get_n_nodes();
			for(int part_iter = 0; part_iter<msize; part_iter++){
				ener = membranes[mem_iter].get_node(part_iter).energy;
			results<<"\t\t\t"<<ener<<std::endl;					
			}
		}		
		results<<"</DataArray>\n";
		
		//TOTAL ENERGY
		results<< "<DataArray type=\"Float64\" Name=\""<<"total_energy"<<"\" format=\"ascii\" RangeMin=\""<<-10000<<"\" RangeMax=\""<<10000<<"\">\n";
		for(int mem_iter = 0; mem_iter<n_membranes; mem_iter++){
			int msize = membranes[mem_iter].get_n_nodes();
			for(int part_iter = 0; part_iter<msize; part_iter++){
				tener = membranes[mem_iter].get_node(part_iter).total_energy;
			results<<"\t\t\t"<<tener<<std::endl;					
			}
		}		
		results<<"</DataArray>\n";

		//SPIN
		results<< "<DataArray type=\"Float64\" Name=\""<<"spin"<<"\" format=\"ascii\" RangeMin=\""<<-10<<"\" RangeMax=\""<<10<<"\">\n";
		for(int mem_iter = 0; mem_iter<n_membranes; mem_iter++){
			int msize = membranes[mem_iter].get_n_nodes();
			for(int part_iter = 0; part_iter<msize; part_iter++){
				sp = membranes[mem_iter].get_node(part_iter).spin;
			results<<"\t\t\t"<<sp<<std::endl;					
			}
		}		
		results<<"</DataArray>\n";
		
	results<<"\t\t</PointData>\n";
	
	//Cell data
	results<<"\t\t<CellData>\n";
		//ASYMMETRY
		double as;
		results<< "<DataArray type=\"Float64\" Name=\""<<"asymmetry"<<"\" format=\"ascii\" RangeMin=\""<<-10<<"\" RangeMax=\""<<10<<"\">\n";
		for(int mem_iter = 0; mem_iter<n_membranes; mem_iter++){
			int esize = membranes[mem_iter].get_n_rectangles();
			for(int eiter = 0; eiter<esize; eiter++){
				as = membranes[mem_iter].get_rectangle(eiter).asymmetry;
			results<<"\t\t\t"<<as<<std::endl;					
			}
		}		
		results<<"</DataArray>\n";
	results<<"</CellData>\n";	
	
	results<<"\t\t<Points>\n";
	
	//Position data are coming
	results<<"\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"-100.0\" RangeMax=\"100.0\">\n";
	for(int mem_iter = 0; mem_iter<n_membranes; mem_iter++){
		int msize = membranes[mem_iter].get_n_nodes();
		for(int part_iter = 0; part_iter<msize; part_iter++){
			RealVectorValue p1 = membranes[mem_iter].get_node(part_iter).get_xyz();
			results<<"\t\t\t"<<p1(0)<<" "<<p1(1)<<" "<<p1(2)<<std::endl;
			//std::cout<<"vtk: "<<p1<<std::endl;					
		}
	}
	results<<"\t\t\t</DataArray>\n";
	results<<"\t\t</Points>\n";
	
	//Now Cell data are coming, for regular particles nothing interesting
	//FUTUR: may be interesting for bounded particles
	results<<"\t\t<Cells>\n";
	//connectivity
	results<<"\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\""<<n_rectangles<<"\">\n";
	
	int sofar=0;
	for(int mem_iter = 0 ; mem_iter < n_membranes; ++mem_iter){
		int esize = membranes[mem_iter].get_n_rectangles();
		for(int e_iter=0; e_iter < esize; e_iter++){
			Rectangle el = membranes[mem_iter].get_rectangle(e_iter);
			results<<"\t\t\t"<<el.get(0)+sofar<<" "<<el.get(1)+sofar<<" "<<el.get(2)+sofar<<" "<<el.get(3)+sofar<<std::endl;
		}
		sofar += membranes[mem_iter].get_n_nodes(); //shift by the sum of the n_nodes of prev membranes
	}
    results<<"\t\t\t</DataArray>\n";
    
    //offsets
    results<<"\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" RangeMin=\"1\" RangeMax=\""<<n_nodes<<"\">\n";
	for(int iter =0 ; iter < n_rectangles; ++iter){ 
		results<<"\t\t\t"<<4*(iter+1)<<std::endl;
	}  
	results<<"\t\t\t</DataArray>\n";
	//types: 1 for regular particles
	// 5 for triangles
	//9 for quad
	results<<"\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\" RangeMin=\"1\" RangeMax=\"9\">\n";
	for(int iter =0 ; iter < n_rectangles; ++iter){ 
		results<<"\t\t\t"<<9<<std::endl;
	}  
    results<<"\t\t\t</DataArray>\n";
	results<<"\t\t</Cells>\n";

//end of file
    results<<"\t</Piece>\n";
	results<<"\t</UnstructuredGrid>\n";
	results<<"</VTKFile>";
	results.close();
	

//WRITE THE PVTU file
		presults<<"<?xml version=\"1.0\"?>\n";
		presults<<"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
		presults<<"<PUnstructuredGrid GhostLevel=\"0\">\n";
		presults<<"<PPointData>\n";
		presults<<"<PDataArray type=\"Float64\" Name=\"energy\"/>\n";
		presults<<"<PDataArray type=\"Float64\" Name=\"total_energy\"/>\n";
		presults<<"<PDataArray type=\"Float64\" Name=\"spin\"/>\n";
		presults<<"</PPointData>\n";
		presults<<"<PCellData>\n";
		presults<<"<PDataArray type=\"Float64\" Name=\"asymmetry\"/>\n";
		presults<<"</PCellData>\n";
		presults<<"<PPoints>\n";
        presults<<"<PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n";
		presults<<"</PPoints>\n";

		
		//std::cout<<vtu_source<<std::endl;

		presults<<"<Piece Source=\""<<vtu_source<<"\"/>\n";
		presults<<"</PUnstructuredGrid>\n";
		presults<<"</VTKFile>";	
		presults.close();
	
}

void MembraneSystem::construct(){
	//since off file does not include bond data
	// we take care of it
	sample.find_bonds(alpha_1, alpha_2, alpha_3, d_bond, d_bond1, d_bond2, has_seam);
	sample.check_orientations(mesh_type);
	//dimer stucture and cooperativity
	if(has_seam) sample.alternate_seam(coop_coef); //insert strength of cooperativity
	else sample.alternate(coop_coef);
	
	sample.two_diagonals();
	sample.check_orientations(mesh_type); //recheck again for new bonds
	
	membranes[0] = sample;
	
	n_nodes = membranes[0].get_n_nodes();
	n_triangles = membranes[0].get_n_triangles();
	n_bonds =membranes[0].get_n_bonds();
	/*//For TAO
	ix_kds = new PetscInt[4*n_nodes];
	y_kds = new PetscScalar[4*n_nodes];
	for (int nit = 0; nit<n_nodes; nit++){
			ix_kds[4*nit + 0] = 4*nit + 0;
			ix_kds[4*nit + 1] = 4*nit + 1;
			ix_kds[4*nit + 2] = 4*nit + 2;		
			ix_kds[4*nit + 3] = 4*nit + 3;		
	}
	*/
	
}
void MembraneSystem::construct_from_rect(){
	
	sample.construct(coop_coef);
	
	membranes[0] = sample;
	
	n_nodes = membranes[0].get_n_nodes();
	n_rectangles = membranes[0].get_n_rectangles();
	n_triangles = membranes[0].get_n_triangles();
	
}

void MembraneSystem::clear_forces(){
	
	Eextern = 0.0;
	Estretch = 0.0;
	Ebending = 0.0;
	Ebend4 = 0.0;
	Ecoop = 0.0;
	Etotal = 0.0;
	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over different membranes
		Membrane& m1 = membranes[mit];
		int nsize = m1.get_n_nodes();
		
	
		for (int nit = 0; nit<nsize; nit++){ //iterate over triangles of a given membrane
			Node& n1 = m1.get_node(nit);

			n1.force(0) = 0.0;
			//bond_energy += 0.5*(xpos-d1)*(xpos-d1);
			n1.force(1) = 0.0;
			//bond_energy += 0.5*(ypos-d1)*(ypos-d1);
			n1.force(2) = 0.0;
			//bond_energy += 0.5*(zpos-d1)*(zpos-d1);
			n1.energy = 0.0;
			n1.total_energy = 0.0;
			
			n1.fspin = 0.0;
	
		}
		//std::cout<<"bond energy: "<<bond_energy<<std::endl;
	}
}

void MembraneSystem::add_noise(){
	double strength = noise_strength;
	RealVectorValue disp;
	
	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over different membranes
		Membrane& m1 = membranes[mit];
		int nsize = m1.get_n_nodes();
		
		for (int nit = 0; nit<nsize; nit++){ //iterate over triangles of a given membrane
			Node& n1 = m1.get_node(nit);
			
			disp(0) = strength*(r1->doub()-0.5);
			disp(1) = strength*(r1->doub()-0.5);
			disp(2) = strength*(r1->doub()-0.5);
			//std::cout<<"disp1: "<<disp<<std::endl;
			//std::cout<<"disp2: "<<n1.get_xyz() <<std::endl;
			n1.get_xyz() += disp;
			//std::cout<<"disp3: "<<n1.get_xyz() <<std::endl;

		}
		
	}
	std::cout<<"Noise added "<<std::endl;
}
void MembraneSystem::initial_conditions(){
	double strength = 1.0;
	RealVectorValue disp;
	
	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over different membranes
		Membrane& m1 = membranes[mit];
		int nsize = m1.get_n_nodes();
		
		for (int nit = 0; nit<nsize; nit++){ //iterate over triangles of a given membrane
			Node& n1 = m1.get_node(nit);
			
			disp(0) = 0.0;//strength*(r1->doub()-0.5);
			disp(1) = strength*sin(PI*n1.get_xyz()(2)/20.0); 
			disp(2) = 0.0;// strength*(r1->doub()-0.5);
			//std::cout<<"disp1: "<<disp<<std::endl;
			//std::cout<<"disp2: "<<n1.get_xyz() <<std::endl;
			n1.get_xyz() += disp;
			//std::cout<<"disp3: "<<n1.get_xyz() <<std::endl;
			
			
			//spin values
			//n1.spin = r1->doub();
			RealVectorValue xyz = n1.get_xyz0();
			if(xyz(0)>0)
				n1.spin = 1.0;
			else
				n1.spin = 0.0;
		
			//6 protofilaments
			//n1.spin = 0.5*xyz(0)+0.5;
			//8 protofilaments
			//n1.spin = 0.5*(xyz(0)/1.30656)+0.5;
			
		}
		
	
		
	}
	std::cout<<"initial conditions applied "<<std::endl;
}
void MembraneSystem::initial_torsion(){
	
	//double theta = 0.1;
//	double theta = 2.0*PI*4.0/3400.0;
	double theta = 0.0;
//	std::cout<<"PI: "<<PI<<" "<<theta<<std::endl;
	double costh, sinth;
	double xcomp, ycomp, zcomp;
	RealVectorValue disp,current, new_xyz;
	
	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over different membranes
		Membrane& m1 = membranes[mit];
		int nsize = m1.get_n_nodes();
		
		for (int nit = 0; nit<nsize; nit++){ //iterate over triangles of a given membrane
			Node& n1 = m1.get_node(nit);
			
			current = n1.get_xyz();
			xcomp = n1.get_xyz()(0);
			ycomp = n1.get_xyz()(1);
			zcomp = n1.get_xyz()(2);
			
			costh = cos(theta*zcomp);
			sinth = sin(theta*zcomp);
			
			new_xyz(0) = costh*xcomp - sinth*ycomp;
			new_xyz(1) = sinth*xcomp + costh*ycomp;
			new_xyz(2) = zcomp;

			n1.get_xyz() = new_xyz;	
		}
		
	}
	std::cout<<"torsion applied "<<std::endl;
}


void MembraneSystem::reset_bond_lengths(){
	double ds = 0.0;
	
	double sin_radial =  alpha_3; //0.464723;
	double l_diagonal = sqrt(2);
	
	//trick for bond type 3
	this->calculate_rectangle_normals();
	this->calculate_angles();
	
	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over different membranes
		Membrane& m1 = membranes[mit];

		int bsize = m1.get_n_bonds();
		for(int bit =0; bit<bsize; bit++){
			Bond& b1 = m1.get_bond(bit);

			int nd1 = b1.get_node(0);
			int nd2 = b1.get_node(1);
			Node& p1 = m1.get_node(nd1);
			Node& p2 = m1.get_node(nd2);				
			
			double dd= (p2.get_xyz()-p1.get_xyz()).size();
				
			if(b1.bond_type == 1){
				b1.length = dd;
				b1.sin_angle = alpha_1;
				b1.current_sin1 = 0.0; b1.current_sin2 = 0.0;
				if(Srigidity1>0.0){ b1.sCoef = Srigidity1; b1.strOn = true;} //otherwise they are
				if(Brigidity1>0.0){ b1.bCoef = Brigidity1; b1.benOn = true;} //initialized as FALSE
			} 
			if(b1.bond_type == 2){
				b1.length = dd; 
				//b1.sin_angle = -0.588;
				b1.sin_angle = alpha_2;
				b1.current_sin1 = 0.0; b1.current_sin2 = 0.0;				
				if(Srigidity2>0.0){ b1.sCoef = Srigidity2; b1.strOn = true;}
				if(Brigidity2>0.0){ b1.bCoef = Brigidity2; b1.benOn = true;}
			} 
			if(b1.bond_type == 3){
				b1.length = dd;
				b1.sin_angle = b1.current_sin1;
				//b1.sin_angle = -0.5;
				b1.current_sin1 = sin_radial; b1.current_sin2 = sin_radial;	
				if(Srigidity3>0.0){ b1.sCoef = Srigidity3; b1.strOn = true;}
				if(Brigidity3>0.0){ b1.bCoef = Brigidity3; b1.benOn = true;}
			} 
			if(b1.bond_type == 4){
				//b1.length = l_diagonal;
				b1.length = dd;
				b1.sin_angle = alpha_4;
				b1.current_sin1 = 0.0; b1.current_sin2 = 0.0;				
				if(Srigidity4>0.0){ b1.sCoef = Srigidity4; b1.strOn = true;}
				if(Brigidity4>0.0){ b1.bCoef = Brigidity4; b1.benOn = true;}
			} 
			if(b1.bond_type == 5){
				//b1.length = l_diagonal;
				b1.length = dd;
				b1.sin_angle = alpha_5;
				b1.current_sin1 = 0.0;b1.current_sin2 = 0.0;
				if(Srigidity5>0.0){ b1.sCoef = Srigidity5; b1.strOn = true;}
				if(Brigidity5>0.0){ b1.bCoef = Brigidity5; b1.benOn = true;}
			} 

				
			//CHECK BOND TYPES 1 and 2 
			if(b1.bond_type == 1){
				p1.spin += 10;
				p2.spin += 10;}
			if(b1.bond_type == 2){
				p1.spin += 100;
				p2.spin += 100;}
			
			//ANGLES	only set for type 2 and 3
/*			if(!b1.is_boundary()){
				

			//	if(b1.bond_type != 4){ 
			if(true){
				int nd3 = b1.get_lnode(0);
				int nd4 = b1.get_lnode(1);
				Node& p3 = m1.get_node(nd3);
				Node& p4 = m1.get_node(nd4);

				RealVectorValue u1 = p3.get_xyz() - p1.get_xyz();
				RealVectorValue u2 = p2.get_xyz() - p1.get_xyz();
				double u2size = u2.size();
				RealVectorValue u2n = u2/u2size;				
				RealVectorValue u3 = p4.get_xyz() - p1.get_xyz();

				RealVectorValue n_alpha = u1.cross(u2);
				double na_size = n_alpha.size();
				n_alpha /= na_size;
				RealVectorValue n_beta = u2.cross(u3);
				double nb_size = n_beta.size();
				n_beta /= nb_size;

				
				RealVectorValue cross_ab = n_alpha.cross(n_beta);
			//	double cross_size = cross_ab.size();
				
				
				b1.sin_angle = cross_ab*u2n;
				b1.current_sin = cross_ab*u2n;	
			//	if(b1.bond_type == 3)
				//		std::cout<<"angle"<<b1.sin_angle<<std::endl;
				
				}			
			}*/
		
		}
	//	m1.print();
	}
}

void MembraneSystem::set_temperature(double ttemp){
	this->temperature = ttemp;
	noise_coef = sqrt(2.0*boltzmann*temperature*gamma*dt);
	
	diffusion_coef = diffusion/(boltzmann*temperature);
	//brownian_coef = sqrt(2.0*diffusion*dt);	
	//std::cout<<"temperature is set to: "<<temperature<<std::endl;

}


void MembraneSystem::internal_forces(){
	Estretch = 0.0; Ebending = 0.0; Ebend4 = 0.0; Ecoop = 0.0;
	Einter = 0.0; Etotal = 0.0;
	
	double en = 0.0;
	double en_coop = 0.0;
	RealVectorValue fs;
	
	//double l1 = 0.0, l2 = 0.0, sp12 = 0.0; 
	//for stretching
	double ds = 0.0;
	RealVectorValue distance;

	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over different membranes
		Membrane& m1 = membranes[mit];
		int bsize = m1.get_n_bonds();
		
		for (int bit = 0; bit<bsize; bit++){ //iterate over bonds
			Bond& b1 = m1.get_bond(bit);
			
			//get a reference to the nodes of the bond
			int nd1 = b1.get_node(0); int nd2 = b1.get_node(1);
			Node& p1 = m1.get_node(nd1); Node& p2 = m1.get_node(nd2);
		
		if(stretching_on && b1.strOn){ //STRETCHING BEGIN		
			distance = -b1.u2n;
			ds = b1.u2size;
			/*distance = p1.get_xyz() - p2.get_xyz();
			ds = distance.size();
			distance/=ds;*/
			double sCoef = b1.sCoef;		
			double d1 = b1.get_length();
			
			fs = sCoef*(ds-d1)*distance;
			p1.force += fs;
			p2.force += -fs;
			
			en = 0.5*sCoef*(ds-d1)*(ds-d1);		
			Estretch += en;
			p1.total_energy += en;
			p2.total_energy += en;
		
		}//STRETCHING END
		
		if(bending_on && b1.benOn){//BENDING BEGIN
			//check if the bond is on boundary
			//if(!b1.is_boundary() && b1.bond_type == 1){
			if(!b1.is_boundary()){
				double angle0 = b1.get_target_angle();	
				Rectangle& r1 = m1.get_rectangle(b1.rectangles[0]);
				Rectangle& r2 = m1.get_rectangle(b1.rectangles[1]);
				int nd3, nd4, nd3s, nd4s;
				double diff;
				RealVectorValue u1,u3;
				/*
				 * APPLY HARMONIC POTENTIAL if we have bond-type 2 or 3
				 */
				//if(b1.bond_type != 1){ 
			//if(b1.bond_type != 2){ 
		
			//if(b1.bond_type==1){ 				
			if(b1.bond_type==1 || (b1.bond_type==2 && harmonic)){ 								
				// reference to secondary nodes
				nd3 = r1.get(1); nd4 = r2.get(3);
				Node& p3 = m1.get_node(nd3); Node& p4 = m1.get_node(nd4);
				
				u1 = p3.get_xyz()-p1.get_xyz();					
				u3 = p4.get_xyz()-p1.get_xyz();					
	
				diff = b1.current_sin1 - angle0;
					
				b1.bending_tensors(diff, r1.normals[1], r2.normals[0], u1, u3);
					
				p1.force += b1.bforce1;				
				p2.force += b1.bforce2;
				p3.force += b1.bforce3;	
				p4.force += b1.bforce4;
				en = 0.5*b1.bCoef*diff*diff;
				p1.total_energy += en;
				p2.total_energy += en;
				p3.total_energy += en;
				p4.total_energy += en;
				Ebending += en;
				//
				nd3s = r1.get(0); nd4s = r2.get(2);
				Node& p3s = m1.get_node(nd3s); Node& p4s = m1.get_node(nd4s);
				
				u1 = p3s.get_xyz()-p1.get_xyz();					
				u3 = p4s.get_xyz()-p1.get_xyz();					
	
				diff = b1.current_sin2 - angle0;
					
				b1.bending_tensors(diff, r1.normals[2], r2.normals[3], u1, u3);
					
				p1.force += b1.bforce1;				
				p2.force += b1.bforce2;
				p3s.force += b1.bforce3;	
				p4s.force += b1.bforce4;
				en = 0.5*b1.bCoef*diff*diff;
				p1.total_energy += en;
				p2.total_energy += en;
				p3s.total_energy += en;
				p4s.total_energy += en;
				Ebending += en;
				//	
			}
			if(b1.bond_type==3){ 				
				// reference to secondary nodes
				nd3 = r2.get(1); nd4 = r1.get(3);
				Node& p3 = m1.get_node(nd3); Node& p4 = m1.get_node(nd4);
				
				u1 = p3.get_xyz()-p1.get_xyz();					
				u3 = p4.get_xyz()-p1.get_xyz();					
	
				diff = b1.current_sin1 - angle0;
					
				b1.bending_tensors(diff, r2.normals[0], r1.normals[1], u1, u3);
					
				p1.force += b1.bforce1;				
				p2.force += b1.bforce2;
				p3.force += b1.bforce3;	
				p4.force += b1.bforce4;
				en = 0.5*b1.bCoef*diff*diff;
				p1.total_energy += en;
				p2.total_energy += en;
				p3.total_energy += en;
				p4.total_energy += en;
				Ebending += en;
				//
				nd3 = r2.get(2); nd4 = r1.get(0);
				Node& p3s = m1.get_node(nd3); Node& p4s = m1.get_node(nd4);
				
				u1 = p3s.get_xyz()-p1.get_xyz();					
				u3 = p4s.get_xyz()-p1.get_xyz();					
				//std::cout<<"Nodes: "<<nd1<<" "<<nd2<<std::endl;
				diff = b1.current_sin2 - angle0;
					
				b1.bending_tensors(diff, r2.normals[2], r1.normals[3], u1, u3);
					
				p1.force += b1.bforce1;				
				p2.force += b1.bforce2;
				p3s.force += b1.bforce3;	
				p4s.force += b1.bforce4;
				en = 0.5*b1.bCoef*diff*diff;
				p1.total_energy += en;
				p2.total_energy += en;
				p3s.total_energy += en;
				p4s.total_energy += en;
				Ebending += en;
				//	
			}
			if(b1.bond_type==4){ 				
				//// reference to secondary nodes
				nd3 = r1.get(2); nd4 = r1.get(0);
				Node& p3 = m1.get_node(nd3); Node& p4 = m1.get_node(nd4);
				
				u1 = p3.get_xyz()-p1.get_xyz();					
				u3 = p4.get_xyz()-p1.get_xyz();					
	
				diff = b1.current_sin1 - angle0;
					
				b1.bending_tensors(diff, r1.normals[1], r1.normals[0], u1, u3);
					
				p1.force += b1.bforce1;				
				p2.force += b1.bforce2;
				p3.force += b1.bforce3;	
				p4.force += b1.bforce4;
				en = 0.5*b1.bCoef*diff*diff;
				p1.total_energy += en;
				p2.total_energy += en;
				p3.total_energy += en;
				p4.total_energy += en;
				Ebending += en;
			}
			if(b1.bond_type==5){ 				
				// reference to secondary nodes
				nd3 = r1.get(1); nd4 = r1.get(3);
				Node& p3 = m1.get_node(nd3); Node& p4 = m1.get_node(nd4);
				
				u1 = p3.get_xyz()-p1.get_xyz();					
				u3 = p4.get_xyz()-p1.get_xyz();					
	
				diff = b1.current_sin1 - angle0;
					
				b1.bending_tensors(diff, r1.normals[3], r1.normals[2], u1, u3);
					
				p1.force += b1.bforce1;				
				p2.force += b1.bforce2;
				p3.force += b1.bforce3;	
				p4.force += b1.bforce4;
				en = 0.5*b1.bCoef*diff*diff;
				p1.total_energy += en;
				p2.total_energy += en;
				p3.total_energy += en;
				p4.total_energy += en;
				Ebending += en;
			}
			
			/*
			 * APPLY QUARTIC POTENTIAL if the bond-type is 2
			 */			
			if(b1.bond_type==2 && !harmonic){
			//if(false){
				
				double pCprime = pC+b1.mu_coop;
				double pX, coangle, pD, pCoef;
///////////				
				coangle = b1.coop_sin1; 
				pX = b1.current_sin1;
				pD = -2.0*b1.mu_coop*coangle;
				
				pCoef= 4.0*pA* pX*pX*pX + 3.0*pB* pX*pX + 2.0*pCprime* pX + pD;

				nd3 = r1.get(1); nd4 = r2.get(3);
				Node& p3 = m1.get_node(nd3); Node& p4 = m1.get_node(nd4);

				u1 = p3.get_xyz()-p1.get_xyz();					
				u3 = p4.get_xyz()-p1.get_xyz();					

				b1.bending_tensors_for_anharmonic(pX, pCoef, r1.normals[1], r2.normals[0], u1, u3);
	
				p1.force += b1.bforce1;				
				p2.force += b1.bforce2;
				p3.force += b1.bforce3;	
				p4.force += b1.bforce4;

				//FIX THE ENERGY
				en = pA*pX*pX*pX*pX + pB*pX*pX*pX + pC*pX*pX;
				//cooperativity energy
				en_coop = b1.mu_coop*(pX*pX - 2.0*pX*b1.coop_sin1 + b1.coop_residual1);

				p1.energy += en;
				p2.energy += en;

				p1.total_energy += en + en_coop;
				p2.total_energy += en + en_coop;
				p3.total_energy += en + en_coop;
				p4.total_energy += en + en_coop;
				
				//Ebending += en;	
				Ebend4 += en;	
				Ecoop += en_coop;
				
/////////////////////////				
				coangle = b1.coop_sin2; 
				pX = b1.current_sin2;
				pD = -2.0*b1.mu_coop*coangle;
				
				pCoef = 0.0;
				pCoef= 4.0*pA* pX*pX*pX + 3.0*pB* pX*pX + 2.0*pCprime* pX + pD;

				nd3 = r1.get(0); nd4 = r2.get(2);
				Node& p3s = m1.get_node(nd3); Node& p4s = m1.get_node(nd4);

				u1 = p3s.get_xyz()-p1.get_xyz();					
				u3 = p4s.get_xyz()-p1.get_xyz();				

				b1.bending_tensors_for_anharmonic(pX, pCoef, r1.normals[2], r2.normals[3], u1, u3);
	
				p1.force += b1.bforce1;				
				p2.force += b1.bforce2;
				p3s.force += b1.bforce3;	
				p4s.force += b1.bforce4;

				//FIX THE ENERGY
				en = pA*pX*pX*pX*pX + pB*pX*pX*pX + pC*pX*pX;
				//cooperativity energy
				en_coop = b1.mu_coop*(pX*pX - 2.0*pX*b1.coop_sin2 + b1.coop_residual2);

				p1.energy += en;
				p2.energy += en;

				p1.total_energy += en + en_coop;
				p2.total_energy += en + en_coop;
				p3s.total_energy += en + en_coop;
				p4s.total_energy += en + en_coop;
				
				//Ebending += en;	
				Ebend4 += en;	
				Ecoop += en_coop;
////////////								
			}
/////////////////////////////
			}
			
			
		}//BENDING END


		
		}
		//std::cout<<"bond energy: "<<bond_energy<<std::endl;
	}
	Etotal += Estretch + Ebending+ Ebend4 + Einter + Ecoop;// + Elandau; 
}

void MembraneSystem::update_coop_value(double fac){
	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over different membranes
		Membrane& m1 = membranes[mit];
		int bsize = m1.get_n_bonds();
		
		for (int bit = 0; bit<bsize; bit++){ //iterate over bonds
			Bond& b1 = m1.get_bond(bit);
			if(b1.bond_type==2){
				b1.mu_coop *= fac;
				//std::cout<<b1.mu_coop<<std::endl;
			}
		}
	}
	//std::cout<<"cooperativity values update"<<std::endl;		
}
	
void MembraneSystem::calculate_rectangle_normals(){
	
	Membrane& m1 = membranes[0];
	int nd0=-1, nd1=-1, nd2=-1, nd3=-1;
	RealVectorValue normal;
	double nor_size = 0.0;
	
	for (int rit = 0; rit<n_rectangles; rit++){ //iterate over elements
		//std::cout<<"element no: "<<rit<<std::endl;
		Rectangle& r1 = m1.get_rectangle(rit);
		nd0 = r1.get(0); nd1 = r1.get(1); nd2 = r1.get(2); nd3 = r1.get(3);
		//std::cout<<nd0<<" "<<nd1<<" "<<nd2<<" "<<nd3<<" "<<std::endl;
		
		Node& p0 = m1.get_node(nd0);
		Node& p1 = m1.get_node(nd1);
		Node& p2 = m1.get_node(nd2);
		Node& p3 = m1.get_node(nd3);
		
		RealVectorValue e0 = p1.get_xyz() - p0.get_xyz();
		RealVectorValue e1 = p2.get_xyz() - p1.get_xyz();
		RealVectorValue e2 = p3.get_xyz() - p2.get_xyz();
		RealVectorValue e3 = p0.get_xyz() - p3.get_xyz();
		//RealVectorValue d1 = p3.xyz() - p1.xyz();
		//RealVectorValue d2 = p0.xyz() - p2.xyz();
		
		normal = e3.cross(e0);
		nor_size = normal.size();
		normal /= nor_size;
		r1.normals[0] = normal;

		normal = e1.cross(e2);
		nor_size = normal.size();
		normal /= nor_size;
		r1.normals[1] = normal;

		normal = e2.cross(e3);
		nor_size = normal.size();
		normal /= nor_size;
		r1.normals[2] = normal;

		normal = e0.cross(e1);
		nor_size = normal.size();
		normal /= nor_size;
		r1.normals[3] = normal;
		
	}
	
}

		
void MembraneSystem::calculate_normals_and_angles(){

	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over different membranes
		Membrane& m1 = membranes[mit];
		int bsize = m1.get_n_bonds();
		
		for (int bit = 0; bit<bsize; bit++){ //iterate over bonds
			Bond& b1 = m1.get_bond(bit);
			
			//get a reference to the nodes of the bond
			int nd1 = b1.get_node(0);
			int nd2 = b1.get_node(1);
			Node& p1 = m1.get_node(nd1);
			Node& p2 = m1.get_node(nd2);
		
		
		if(bending_on && b1.benOn){//BENDING BEGIN
		
			if(!b1.is_boundary()){
				
				//now we also need reference to secondary nodes
				int nd3 = b1.get_lnode(0);
				int nd4 = b1.get_lnode(1);
			
				Node& p3 = m1.get_node(nd3);
				Node& p4 = m1.get_node(nd4);
				
				//std::cout<<nd1<<" "<<nd2<<" "<<nd3<<" "<<nd4<<" "<<b1.is_longitudinal()<<"\n";
				
				b1.u1 = p3.get_xyz() - p1.get_xyz();
				//u2 is the shared edge
				b1.u2 = p2.get_xyz() - p1.get_xyz();
				b1.u2size = b1.u2.size();
				b1.u2n = b1.u2/b1.u2size; //normalized U2
				 
				b1.u3 = p4.get_xyz() - p1.get_xyz();

				//triangle normals
				b1.n_alpha = b1.u1.cross(b1.u2);
				b1.na_size = b1.n_alpha.size();
				b1.n_alpha /= b1.na_size;
				b1.n_beta = b1.u2.cross(b1.u3);
				b1.nb_size = b1.n_beta.size();
				b1.n_beta /= b1.nb_size;

				
				b1.cross_ab = b1.n_alpha.cross(b1.n_beta);

			/*
			 * APPLY QUARTIC POTENTIAL if the bond-type is 1 (NOW IT IS 4)
			 */			
		//	if(b1.bond_type==2){
			if(true){	
				//double angle0 = b1.get_target_angle();
				//double angle0 = b1.coop_sin;
				double pX = b1.cross_ab*b1.u2n; 
				b1.current_sin = pX;
							
}}}}}}
		
void MembraneSystem::calculate_angles(){
	RealVectorValue kross;
	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over different membranes
		Membrane& m1 = membranes[mit];
		int bsize = m1.get_n_bonds();
		
		for (int bit = 0; bit<bsize; bit++){ //iterate over bonds
			Bond& b1 = m1.get_bond(bit);
			//get a reference to the nodes of the bond
			int nd1 = b1.get_node(0);
			int nd2 = b1.get_node(1);
			Node& p1 = m1.get_node(nd1);
			Node& p2 = m1.get_node(nd2);

			b1.u2 = p2.get_xyz() - p1.get_xyz();
			b1.u2size = b1.u2.size();
			b1.u2n = b1.u2/b1.u2size; //normalized U2 		

		if(!b1.is_boundary()){

			Rectangle& r1 = m1.get_rectangle(b1.rectangles[0]);
			Rectangle& r2 = m1.get_rectangle(b1.rectangles[1]);				

			if(b1.bond_type==1 || b1.bond_type==2){	
				kross = r1.normals[1].cross(r2.normals[0]);
				b1.current_sin1 = kross*b1.u2n;
			//	std::cout<<b1.current_sin1<<std::endl;
				
				kross = r1.normals[2].cross(r2.normals[3]);
				b1.current_sin2 = kross*b1.u2n;
				//std::cout<<b1.current_sin2<<std::endl<<std::endl;				
			}
			if(b1.bond_type==3){
				kross = r2.normals[0].cross(r1.normals[1]);
				b1.current_sin1 = kross*b1.u2n;
				//std::cout<<std::setprecision(9)<<b1.current_sin1<<std::endl;	

				kross = r2.normals[2].cross(r1.normals[3]);
				b1.current_sin2 = kross*b1.u2n;
				//std::cout<<std::setprecision(50)<<b1.current_sin2<<std::endl<<std::endl;				
			}
			if(b1.bond_type==4){
				kross = r1.normals[1].cross(r1.normals[0]);
				b1.current_sin1 = kross*b1.u2n;
				//std::cout<<b1.current_sin1<<std::endl<<std::endl;				
			}
			if(b1.bond_type==5){
				kross = r1.normals[3].cross(r1.normals[2]);
				b1.current_sin1 = kross*b1.u2n;
				//std::cout<<b1.current_sin1<<std::endl<<std::endl;				
			}	
		}
}}}		


void MembraneSystem::landau_forces(){
	Elandau = 0.0;

	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over membranes
		Membrane& m1 = membranes[mit];
		int psize = m1.get_n_nodes();
		
		for(int prit = 0; prit<psize; ++prit){ //iterate over nodes

			if(spin_on){
				Node& p1 = m1.get_node(prit);	
				//f(x) = 16(x^2)(x-1)^2
				double p1s = p1.spin; 
		
				p1.fspin += 32.0*deltaE*p1s*(p1s-1.0)*(2.0*p1s-1.0); 
				Elandau += deltaE *(p1s*p1s)*(p1s-1.0)*(p1s-1.0);
				
			}
		}		
	}
	Etotal += Elandau;	
	
}

//Collects the neigbors' current angles
void MembraneSystem::cooperativity(){
	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over membranes
		Membrane& m1 = membranes[mit];	
		//int bsize = m1.get_n_bonds();
		
		//for (int bit = 0; bit<bsize; bit++){ //iterate over bonds
			//Bond& b1 = m1.get_bond(bit);
			//if(!b1.is_boundary() && b1.bond_type==4){
				//std::cout<<b1.
			//}	
		int b1_idx, bn_idx, bp_idx;
		
		int number_of_dimers = m1.get_n_coop();
		for(int pit = 0; pit<13; pit++){
			for(int dit = 0; dit<number_of_dimers; dit++){		
				b1_idx = m1.get_coop_bond2(pit, dit);
				Bond& b1 = m1.get_bond(b1_idx);
				//b1.coop_sin = -0.2; //test
				if(dit ==0){
					bn_idx = m1.get_coop_bond2(pit, dit+1);
					Bond& b_next = m1.get_bond(bn_idx);
					b1.coop_sin1 = b_next.current_sin1;
					b1.coop_sin2 = b_next.current_sin2;
					b1.coop_residual1 = b_next.current_sin1*b_next.current_sin1;
					b1.coop_residual2 = b_next.current_sin2*b_next.current_sin2;
				}	
				else if(dit ==number_of_dimers-1){
					bp_idx = m1.get_coop_bond2(pit, dit-1);					
					Bond& b_prev = m1.get_bond(bp_idx);
					b1.coop_sin1 = b_prev.current_sin1;
					b1.coop_sin2 = b_prev.current_sin2;
					b1.coop_residual1 = b_prev.current_sin1*b_prev.current_sin2;
					b1.coop_residual2 = b_prev.current_sin2*b_prev.current_sin2;
				}	
				else{
					bp_idx = m1.get_coop_bond2(pit, dit-1);	
					Bond& b_prev = m1.get_bond(bp_idx);
					
					bn_idx = m1.get_coop_bond2(pit, dit+1);
					Bond& b_next = m1.get_bond(bn_idx);
					
					b1.coop_sin1 = 0.5*(b_prev.current_sin1 + b_next.current_sin1);
					b1.coop_sin2 = 0.5*(b_prev.current_sin2 + b_next.current_sin2);
					b1.coop_residual1 = 0.5*(b_next.current_sin1*b_next.current_sin1 + b_prev.current_sin1*b_prev.current_sin1);
					b1.coop_residual2 = 0.5*(b_next.current_sin2*b_next.current_sin2 + b_prev.current_sin2*b_prev.current_sin2);
				}					
			}
		}
		//Iterate over protofilaments and then over dimers of each protofilament
/*		int number_of_dimers = m1.get_n_coop();
		for(int pit = 0; pit<13; pit++){
			for(int dit = 0; dit<number_of_dimers; dit++){
				std::cout<<"do\n";
				Bond* b1 = (m1.get_coop_bond(pit, dit));
				b1->coop_sin = -0.2;
				std::cout<<b1->coop_sin<<"\n";

				//if(dit ==0){
					//Bond* b_next = m1.get_coop_bond(pit, dit+1);
					//b1.coop_sin = b_next->current_sin;
				//}	
				//else if(dit ==number_of_dimers-1){
					//Bond* b_prev = m1.get_coop_bond(pit, dit-1);
					//b1.coop_sin = b_prev->current_sin;
				//}	
				//else{
					//Bond* b_prev = m1.get_coop_bond(pit, dit-1);
					//Bond* b_next = m1.get_coop_bond(pit, dit+1);
					//b1.coop_sin = 0.5*(b_prev->current_sin + b_next->current_sin);
				//}	
			//	std::cout<<"Cooperativity: "<<b1.bond_index<<" "<<pit<<" "<<dit<<"\n";
				
			}
		}*/
		
	}	
}

void MembraneSystem::external_forces(){
	Eextern = 0.0;
	
	if(nodal_on){
		//this->nodal_forces();
		this->deflection_force();
		Eextern += Enodal;
	}
	if(container_on){
		this->container_forces();
		Eextern += Econtainer;
	}
	if(moment_on){
		this->nodal_moment();
		//Eextern += Emoment;
	}
	
	Etotal += Eextern;
	
}

void MembraneSystem::nodal_forces(){
	Enodal = 0.0;
	double ene = 0.0;
	for(int mit=0 ; mit < this->n_membranes; ++mit){
		Membrane& m1 = membranes[mit];
		int bsize = m1.get_n_bonds();
		
		for (int bit = 0; bit<bsize; bit++){ //iterate over bonds
			Bond& b1 = m1.get_bond(bit);
		
		//	if(b1.is_boundary()){
			if(true){
				int nd1 = b1.get_node(0);
				int nd2 = b1.get_node(1);
				Node& p1 = m1.get_node(nd1);
				Node& p2 = m1.get_node(nd2);
				
				
				RealVectorValue x1 = p1.get_xyz();
				RealVectorValue x10 = p1.get_xyz0();
				RealVectorValue x2 = p2.get_xyz();
				RealVectorValue x20 = p2.get_xyz0();
				//RealVectorValue r1;
				double dr;
			/*	if(x10(2)<boundary1){				
					//std::cout<<"external points: "<<bit<<" "<<nd1<<" "<<nd2<<std::endl;
						
					r1 = x1-x10;
					dr = r1.size();
					
					if(dr!=0.0){  //1st node 
						r1 /= dr;
					
						p1.force += mu_nodal*dr*r1;
					
						ene = 0.5*mu_nodal*dr*dr;
						p1.total_energy += ene;
						Enodal += ene;
					}
					
					r1 = x2-x20;
					dr = r1.size();
					if(dr!=0.0){	//2nd node  
						r1 /= dr;
										
						p2.force += mu_nodal*dr*r1;
					
						ene = 0.5*mu_nodal*dr*dr;
						p2.total_energy += ene;
						Enodal += ene;
					}
				}*/
				//other end
				//if(x10(2)>37.5){				
				if(x10(2)>boundary2){				
		//		if(true){				
					//double zz = 19.0;
					
					p1.force(1) += -force_strength;
					p2.force(1) += -force_strength;
				/*	
					Node& p3 = m1.get_node(nd1-8);
					Node& p4 = m1.get_node(nd2-8);
					p3.force(2) += -force_strength;
					p3.force(2) += -force_strength;
			*/
			//RADIAL
			/*		p1.force(0) += -force_strength*x10(0);
					p1.force(1) += -force_strength*x10(1);
					p2.force(0) += -force_strength*x20(0);
					p2.force(1) += -force_strength*x20(1);
				*/	
					
					
				/*	p1.force(2) += mu_nodal*(x1(2)-zz);
					p1.energy += 0.5*mu_nodal*(x1(2)-zz)*(x1(2)-zz);
					Enodal += 0.5*mu_nodal*(x1(2)-zz)*(x1(2)-zz);
					
					p2.force(2) += mu_nodal*(x2(2)-zz);
					p2.energy += 0.5*mu_nodal*(x2(2)-zz)*(x2(2)-zz);
					Enodal += 0.5*mu_nodal*(x2(2)-zz)*(x2(2)-zz);
				*/
				
				
				}
			}
		}		
	}	
}


void MembraneSystem::nodal_moment(){
	//Enodal = 0.0;
	//double ene = 0.0;
	
	double radius=0.0;
	double costh = 0.0, sinth = 0.0;
	RealVectorValue moment;
	RealVectorValue cntr;
	
	for(int mit=0 ; mit < this->n_membranes; ++mit){
		Membrane& m1 = membranes[mit];
		int psize = m1.get_n_nodes();
		
		
		//Find ring center
		for (int pit = 0; pit<psize; pit++){ //iterate over bonds
			Node& p1 = m1.get_node(pit);
			
			RealVectorValue x1 = p1.get_xyz();
			RealVectorValue x10 = p1.get_xyz0();
						
			if(x10(2)>moment_boundary_down && x10(2)<moment_boundary_up){				
				cntr += x1;
			}
		}
		cntr /= 13.0;
		
		for (int pit = 0; pit<psize; pit++){ //iterate over bonds
			Node& p1 = m1.get_node(pit);
			
			RealVectorValue x1 = p1.get_xyz();
			RealVectorValue x10 = p1.get_xyz0();
						
			if(x10(2)>moment_boundary_down && x10(2)<moment_boundary_up){				

				//p1.force(1) += -force_strength;
				costh = x1(0)-cntr(0); sinth = x1(1)-cntr(1);
				radius = sqrt(costh*costh + sinth*sinth);
				costh /= radius; sinth /= radius;
				moment(0) = -sinth;
				moment(1) = costh;
				moment(2) = 0.0;
				p1.force += moment_strength*moment; 
				
			}
		}
	}
}

void MembraneSystem::update_parameters(int nit){
	if(force_is_gradual)
		force_strength = nit*update_coef*max_force_strength;
	else
		force_strength = max_force_strength;
	
	if(moment_is_gradual)
		moment_strength = nit*update_coef*max_moment_strength;
	else
		moment_strength = max_moment_strength;
	
	//if(nit%printeach == 0)	
	    //std::cout<<"mm: "<<moment_strength<<"\n";
}

void MembraneSystem::deflection_force(){
	//std::cout<<"force strength is: "<<force_strength<<std::endl;
	Enodal = 0.0;
	double ene1 = 0.0;
	
	for(int mit=0 ; mit < this->n_membranes; ++mit){
		int psize = membranes[mit].get_n_nodes();
		
		for(int prit = 0; prit<psize; ++prit){
			//get a reference to the particle
			Node& p1 = membranes[mit].get_node(prit);
				
				RealVectorValue x1 = p1.get_xyz();
				RealVectorValue x10 = p1.get_xyz0();
				
				if(x10(2)>boundary2){				
					p1.force(1) += -force_strength;
					ene1 = (x1(1)-x10(1))*(-force_strength); 
					Enodal += ene1;
				}	
}}}

void MembraneSystem::container_forces(){
	Econtainer = 0.0;

	//std::cout<<"rad "<<cradius<<std::endl;
	//std::cout<<"mu "<<mu_cont<<std::endl;
	RealVectorValue radial;
	double en, rsize;
	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over membranes
		Membrane& m1 = membranes[mit];
		int psize = m1.get_n_nodes();
		
		for(int prit = 0; prit<psize; ++prit){ //iterate over nodes

			Node& p1 = m1.get_node(prit);	
			
			radial = p1.get_xyz();
			rsize = radial.size();
			if (rsize>cradius){
				p1.force += mu_cont*(rsize-cradius)*radial;
				
				en = 0.5*(rsize-cradius)*(rsize-cradius);
				p1.energy+=en;
				Econtainer += en;
				
			}
			
		}		
	}	
	
}


void MembraneSystem::damping(){
	
	RealVectorValue noise;
	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over membranes
		Membrane& m1 = membranes[mit];
		int psize = m1.get_n_nodes();
		
		for(int prit = 0; prit<psize; ++prit){ //iterate over nodes

			Node& p1 = m1.get_node(prit);	
				
				noise(0) = r1->gauss(noise_coef);
				noise(1) = r1->gauss(noise_coef);
				noise(2) = r1->gauss(noise_coef);

				p1.force += gamma*mass*p1.velocity - noise;

				//p1.fspin += gamma*p1.vel_spin;
			
		}		
	}
	//Etotal += Eextern;
}


void MembraneSystem::brownian_update(int nit){
	
	RealVectorValue new_position;
	RealVectorValue noise;
	
	for(int mit=0 ; mit < this->n_membranes; ++mit){
		int psize = membranes[mit].get_n_nodes();
		
		for(int prit = 0; prit<psize; ++prit){
			//get a reference to the particle
			Node& p1 = membranes[mit].get_node(prit);
			
			//generate noise
				//noise(0) = brownian_coef*r1->gauss(1.0);
				//noise(1) = brownian_coef*r1->gauss(1.0);
				//noise(2) = brownian_coef*r1->gauss(1.0);
				noise(0) = r1->gauss( brownian_coef);
				noise(1) = r1->gauss( brownian_coef);
				noise(2) = r1->gauss( brownian_coef);
			//std::cout<<noise<<std::endl;				
			//actually force is the gradient(= -force)
			//so we take -force to get actual force
			
			new_position = p1.position - p1.force*diffusion_coef*dt + noise;
			//if(prit == 39 && nit%100 ==0) std::cout<<p1.force<<std::endl;
			 	
			p1.p_position = p1.position;
			p1.position = new_position;
		
			p1.p_force = -p1.force;
			p1.force.zero();			
		}
	}	
}

void MembraneSystem::constrained_brownian_update(int nit){
	
	RealVectorValue new_position;
	RealVectorValue noise;
	
	for(int mit=0 ; mit < this->n_membranes; ++mit){
		int psize = membranes[mit].get_n_nodes();
		
		for(int prit = 0; prit<psize; ++prit){
			//get a reference to the particle
			Node& p1 = membranes[mit].get_node(prit);
			
			RealVectorValue x10 = p1.get_xyz0();
			
			//if(x10(2)>boundary1 || x10(2)<boundary0){	
			if(x10(2)>boundary1){
			//if(x10(2)<boundary1){
				noise(0) = r1->gauss( brownian_coef);
				noise(1) = r1->gauss( brownian_coef);
				noise(2) = r1->gauss( brownian_coef);
			
				//actually force is the gradient(= -force)
				//so we take -force to get actual force
			
				new_position = p1.position - p1.force*diffusion_coef*dt + noise;
				
			 	
				p1.p_position = p1.position;
				p1.position = new_position;
			}
				//even i do not update the position of some particles, i keep forces
				p1.p_force = -p1.force;
				p1.force.zero();
					
		}
	}	
}


void MembraneSystem::velocity_verlet_position(int nit){


	RealVectorValue noise;
	RealVectorValue new_position;
	double newspin;

	for(int mit=0 ; mit < this->n_membranes; ++mit){
		int psize = membranes[mit].get_n_nodes();
		
		for(int prit = 0; prit<psize; ++prit){
			//get a reference to the particle
			Node& p1 = membranes[mit].get_node(prit);

			RealVectorValue x10 = p1.get_xyz0();
			
			//if(x10(2)>boundary1){
			if(x10(2)>boundary1 || x10(2)<boundary0){
			
				noise(0) = r1->gauss(noise_coef);
				noise(1) = r1->gauss(noise_coef);
				noise(2) = r1->gauss(noise_coef);							
				//actually force is the gradient(= -force)
				//so we take -force to get actual force
				
				p1.velocity_half = p1.velocity + (-p1.force - gamma*p1.mass*p1.velocity + noise)*dt/(2.0*p1.mass);
				new_position = p1.position + p1.velocity_half*dt;
		 	
				p1.p_position = p1.position;
				p1.position = new_position;
			}
			
			p1.p_force = -p1.force;
			p1.force.zero();
			
			

		}
	}	
}

void MembraneSystem::velocity_verlet_velocity(int nit){
	RealVectorValue noise;
	RealVectorValue new_velocity;
	double new_v_spin;
	
	kinetic_energy=0.0;
	for(int mit=0 ; mit < this->n_membranes; ++mit){
		int psize = membranes[mit].get_n_nodes();
		for(int prit = 0; prit<psize; ++prit){
			//get a reference to the particle
			Node& p1 = membranes[mit].get_node(prit);

			RealVectorValue x10 = p1.get_xyz0();
			
			//if(x10(2)>boundary1){
			if(x10(2)>boundary1 || x10(2)<boundary0){

				noise(0) = r1->gauss(noise_coef);
				noise(1) = r1->gauss(noise_coef);
				noise(2) = r1->gauss(noise_coef);				
				//actually force is the gradient(= -force)
				//so we take -force to get actual force						
				//new_velocity = p1.velocity + (p1.p_force - p1.force)*(dt/(2.0*p1.mass));
				new_velocity = ( p1.velocity_half + (-p1.force+noise)*(dt/(2.0*p1.mass)) )/(1.0 + dt*gamma/2.0); 
	
				p1.velocity = new_velocity;
				kinetic_energy += 0.5*p1.mass* new_velocity*new_velocity;
			
			}
		}
	}
	//kinetic_energy *= (mass/n_particles);
	//potential_energy /= n_particles;
}

void MembraneSystem::write_measurements(int nit, int every) {
		
		//measurements

		measurements<<std::scientific<<nit*dt<<"\t"<<kinetic_energy
				<<"\t"<<Estretch<<"\t"<<Ebending<<"\t"<<Ebend4<<"\t"
					<<Eextern<<"\t"<<Ecoop<<"\t"<<Etotal+kinetic_energy<<std::endl;
			
}

void MembraneSystem::write_tube_data(int step){
	int dum=step;
	int n_dec=0;
	//find out how many digit is in time
	while(dum!=0){
		dum /=10; 
		n_dec++;
		}
		
	std::stringstream file_name;
	file_name<<directory<<"/tube-";
	if (n_dec == 0) n_dec++; //otherwise after adding t=0 to file name, we have an extra digit
	while(n_dec != 8){
		file_name<<"0";
		n_dec++;
		}
	file_name<<step;	
	file_name<<".dat";	
	std::ofstream tubefile(file_name.str().c_str(),std::ios::out);
	//tubefile<<"##those are actually the results of the previous step\n";
	
	int b1_idx;

	//double limit_angle = -0.05869; //barrier: -3.36 degree;
	int state = 0;
	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over membranes
		Membrane& m1 = membranes[mit];	

		int number_of_dimers = m1.get_n_coop();
		//iterate over the height
		for(int dit = 0; dit<number_of_dimers; dit++){
			//iterate over the section
			for(int pit = 0; pit<13; pit++){
				//tubefile<<dit<<"-"<<pit<<" ";				
				b1_idx = m1.get_coop_bond2(pit, dit);
				Bond& b1 = m1.get_bond(b1_idx);
				if((0.5*(b1.current_sin1+b1.current_sin2))<limit_angle) state = 1; else state =0;
				int nd1 = b1.get_node(0);int nd2 = b1.get_node(1);Node& p1 = m1.get_node(nd1);Node& p2 = m1.get_node(nd2);
				//tubefile<<asin(b1.current_sin)*180.0/PI<<"("<<nd1<<"-"<<nd2<<")\t";
				//tubefile<<state<<"("<<nd1<<"-"<<nd2<<")\t";
				tubefile<<state<<"\t";
			}
			tubefile<<"\n";
		}
	}
}

void MembraneSystem::write_bond_data(int step){
	int dum=step;
	int n_dec=0;
	//find out how many digit is in time
	while(dum!=0){
		dum /=10; 
		n_dec++;
		}
		
	std::stringstream file_name;
	file_name<<directory<<"/tbonds-";
	if (n_dec == 0) n_dec++; //otherwise after adding t=0 to file name, we have an extra digit
	while(n_dec != 8){
		file_name<<"0";
		n_dec++;
		}
	file_name<<step;	
	file_name<<".dat";	
	std::ofstream bondfile(file_name.str().c_str(),std::ios::out);
	
	bondfile<<"bond_index \t type \t is_boundary \t target_length \t current_length \t target_angle \t current_angle\n";
	
	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over membranes
		Membrane& m1 = membranes[mit];	
		int bsize = m1.get_n_bonds();
		for(int bit =0; bit<bsize; bit++){
			Bond& b1 = m1.get_bond(bit);
		
			int nd1 = b1.get_node(0);
			int nd2 = b1.get_node(1);
			Node& p1 = m1.get_node(nd1);
			Node& p2 = m1.get_node(nd2);
			
			RealVectorValue	distance = p1.get_xyz() - p2.get_xyz();
			double ds = distance.size();
			double d1 = b1.get_length();
			
			double angle0 = b1.get_target_angle();
			double angle1 = b1.current_sin1;	
			double angle2 = b1.current_sin2;	

			bondfile<<b1.bond_index<<"\t"<<b1.bond_type<<" "<<b1.is_boundary()<<" ";
			bondfile<<d1<<" "<<ds<<" ";
			bondfile<<angle0<<" "<<angle1<<" "<<angle2<<"\n";
			
		}

	}
}

void MembraneSystem::collect_asymmetry_from_diagonals(){
	double as=0.0;
	
	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over membranes
		Membrane& m1 = membranes[mit];	
		int bsize = m1.get_n_bonds();		
		int rsize = m1.get_n_rectangles();		
		
		//set asymmetries to zero
		for(int rit=0; rit<rsize; rit++){
			Rectangle& r1 = m1.get_rectangle(rit);
			r1.asymmetry = 0.0;
		}
		
		//iterate_over_diagonals
		for(int bit =0; bit<bsize; bit++){ 
			Bond& b1 = m1.get_bond(bit);
			if(b1.bond_type==4){
				as = b1.u2size;
				//std::cout<<"bond 4: "<<b1.bond_index<<" "<<as<<" "<<b1.rectangles[0]<<std::endl;
				Rectangle& rect = m1.get_rectangle(b1.rectangles[0]);
				//std::cout<<rect.asymmetry<<std::endl;
				rect.asymmetry += as;
				//std::cout<<rect.asymmetry<<std::endl;				
			}
			if(b1.bond_type==5){
				as = b1.u2size;
				Rectangle& rect = m1.get_rectangle(b1.rectangles[0]);
				rect.asymmetry -= as;
			}
		}
		
		
		
		}

}


void MembraneSystem::test(){
	//std::cout<<&membranes[0].get_node(1).get_xyz()<<std::endl;
	//membranes[0].get_face(1).print();
	std::cout<<"reference to MembraneSystem functions is working\n";
	}


/* //COMMENTED OUT TAO
void MembraneSystem::push_solution(){
	//PetscInt ix_kds[3*n_nodes];
	//PetscScalar y_kds[3*n_nodes];
	
	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over different membranes
		Membrane& m1 = membranes[mit];
		int nsize = m1.get_n_nodes();
		
		for (int nit = 0; nit<nsize; nit++){ //iterate over triangles of a given membrane
			Node& n1 = m1.get_node(nit);
			
			y_kds[4*nit + 0] = n1.get_xyz()(0);
			y_kds[4*nit + 1] = n1.get_xyz()(1);
			y_kds[4*nit + 2] = n1.get_xyz()(2);
			y_kds[4*nit + 3] = n1.spin;
		}
	}	
	
	//for (int iii  = 0; iii<3*n_nodes; iii++)
		//std::cout<<"index: "<<ix_kds[iii]<<"\n";
	
	
	VecSetValues(x,4*n_nodes,ix_kds,y_kds,INSERT_VALUES);
}

void MembraneSystem::pull_solution(){
	//PetscInt ix_kds[3*n_nodes];
	//PetscScalar y_kds[3*n_nodes];
	
	VecGetValues(x,4*n_nodes,ix_kds,y_kds);

	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over different membranes
		Membrane& m1 = membranes[mit];
		int nsize = m1.get_n_nodes();
		
		for (int nit = 0; nit<nsize; nit++){ //iterate over triangles of a given membrane
			Node& n1 = m1.get_node(nit);
			
			n1.get_xyz()(0) = y_kds[4*nit + 0];
			n1.get_xyz()(1) = y_kds[4*nit + 1];
			n1.get_xyz()(2) = y_kds[4*nit + 2];
			n1.spin = y_kds[4*nit + 3];
		}
	}	

	//for (int iii  = 0; iii<3*n_nodes; iii++)
		//std::cout<<"index: "<<y_kds[iii]<<"\n";	
}

void MembraneSystem::pull_solution_sol(PetscScalar *sl[]){
	//PetscInt ix_kds[3*n_nodes];
	//PetscScalar y_kds[3*n_nodes];
	
	//VecGetValues(x,3*n_nodes,ix_kds,y_kds);

	for(int mit=0 ; mit < this->n_membranes; ++mit){ //iterate over different membranes
		Membrane& m1 = membranes[mit];
		int nsize = m1.get_n_nodes();
		
		for (int nit = 0; nit<nsize; nit++){ //iterate over triangles of a given membrane
			Node& n1 = m1.get_node(nit);
			
			n1.get_xyz()(0) = (*sl)[4*nit+0]; 
			n1.get_xyz()(1) = (*sl)[4*nit+1]; 
			n1.get_xyz()(2) = (*sl)[4*nit+2]; 
			n1.spin = (*sl)[4*nit+3]; 
		}
	}	

	//for (int iii  = 0; iii<3*n_nodes; iii++)
		//std::cout<<"index: "<<y_kds[iii]<<"\n";	
}

void MembraneSystem::function_gradient(PetscScalar *sol[], double *ff, PetscScalar *grad[]){
	
	//this->pull_solution(); //transfer current solution to the mesh data
	this->pull_solution_sol(sol); //transfer current solution to the mesh data
	
	this->clear_forces();

	this->internal_forces();
	this->landau_forces();
	this->external_forces();
	
	PetscInt i;

	for (i=0; i<n_nodes; i++){
	    Membrane& m1 = membranes[0];
	    Node& n1 = m1.get_node(i);
    
		(*grad)[4*i + 0] = n1.force(0);   //GRADIENT = FORCE
		(*grad)[4*i + 1] = n1.force(1);
		(*grad)[4*i + 2] = n1.force(2);
		(*grad)[4*i + 3] = n1.fspin;
     }
	
	//assign FUNCTION value
	   
	//std::cout<<"01_Bond energy: "<<bond_energy<<std::endl;
	//std::cout<<"02_ff: "<<*ff<<std::endl;
	*ff = Etotal;
	
	//output
	info = TaoGetSolutionStatus(tao, &tao_iter, &fdummy, &gnormm, &cnormm, &xdifff, &reason);
	if(tao_iter>tao_iter0){	
		this->vtk_writer(vtk_iter);//std::cout<<"ITER: "<<tao_iter<<" "<<vtk_iter<<std::endl;	
		this->write_measurements(vtk_iter, 0);
		vtk_iter++;
	}
	tao_iter0= tao_iter;
}

//take the pointers of argc and argv
void MembraneSystem::tao_construct(int* argc, char ***argv){
	char help[] ="no help, for now\n";
	  // Initialize TAO and PETSc 
	PetscInitialize(argc,argv,(char *)0,help);
	TaoInitialize(argc,argv,(char *)0,help);
	
	zero=0.0;
	info = MPI_Comm_size(PETSC_COMM_WORLD,&size);// CHKERRQ(info);
	info = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);// CHKERRQ(info);

	if (size >1) {
		if (rank == 0)
			PetscPrintf(PETSC_COMM_SELF,"This example is intended for single processor use!\n");
			//SETERRQ(1,"Incorrect number of processors");
	}	
	
	// Initialize problem parameters 
	user.n = 2; user.alpha = 99.0;
	user.sys = this;
	
	// Check for command line arguments to override defaults 
  info = PetscOptionsGetInt(PETSC_NULL,"-n",&user.n,&flg);  //CHKERRQ(info);
  info = PetscOptionsGetReal(PETSC_NULL,"-alpha",&user.alpha,&flg); //CHKERRQ(info)
	
	
	  // Allocate vectors for the solution and gradient 
 // info = VecCreateSeq(PETSC_COMM_SELF,user.n,&x);// CHKERRQ(info);
  info = VecCreateSeq(PETSC_COMM_SELF,4*n_nodes,&x);// CHKERRQ(info);

  // Allocate storage space for Hessian matrix; 
   //  Hessian information is optional -- unless a Newton method is selected
  
  info = MatCreateSeqBAIJ(PETSC_COMM_SELF,2,user.n,user.n,1,PETSC_NULL,&H); //CHKERRQ(info);
  info = MatSetOption(H,MAT_SYMMETRIC,PETSC_TRUE); //CHKERRQ(info);

  // Create TAO solver with desired solution method 
	info = TaoCreate(PETSC_COMM_SELF,"tao_lmvm",&tao); //CHKERRQ(info);
//  info = TaoCreate(PETSC_COMM_SELF,"tao_cg",&tao); //CHKERRQ(info);
	info = TaoSetMaximumIterates(tao, 100000);
	info = TaoSetMaximumFunctionEvaluations(tao, 100000);

  info = TaoApplicationCreate(PETSC_COMM_SELF,&taoapp); //CHKERRQ(info

	
	std::cout<<"tao constructed"<<std::endl;

}



void MembraneSystem::tao_solve(){

  // Set solution vec and an initial guess 
	info = VecSet(x, zero); //CHKERRQ(info); //VecView(x,	PETSC_VIEWER_STDOUT_WORLD );
	
	//Set solution vec according to mesh data
  	push_solution(); //CHKERRQ(info); //VecView(x,	PETSC_VIEWER_STDOUT_WORLD );
	
  info = TaoAppSetInitialSolutionVec(taoapp,x); //CHKERRQ(info); 

  // Set routines for function, gradient, hessian evaluation 
  info = TaoAppSetObjectiveAndGradientRoutine(taoapp,FormFunctionGradient,(void *)&user);  //CHKERRQ(info);
  info = TaoAppSetHessianMat(taoapp,H,H); //CHKERRQ(info);
  info = TaoAppSetHessianRoutine(taoapp,FormHessian,(void *)&user); //CHKERRQ(info);


  // Check for TAO command line options 
  info = TaoSetOptions(taoapp,tao); //CHKERRQ(info);

  // SOLVE THE APPLICATION 
  info = TaoSolveApplication(taoapp,tao); //CHKERRQ(info);


  // Get termination information 
  info = TaoGetTerminationReason(tao,&reason); //CHKERRQ(info);
  std::cout<<"Tao terminate reason: "<<reason<<std::endl;
  if (reason <= 0)
    PetscPrintf(MPI_COMM_WORLD,"Try a different TAO method, adjust some parameters, or check the function evaluation routines\n");

	//output the final result
	info = TaoGetSolutionStatus(tao, &tao_iter, &fdummy, &gnormm, &cnormm, &xdifff, &reason);
	if(tao_iter>tao_iter0) {
		this->vtk_writer(vtk_iter); vtk_iter++;
		//std::cout<<"ITER: "<<tao_iter<<" "<<vtk_iter<<std::endl;	
	}


std::cout<<"tao solved"<<std::endl;

}

void MembraneSystem::tao_destroy(){
  // Free TAO data structures 
  info = TaoDestroy(tao); //CHKERRQ(info);
  info = TaoAppDestroy(taoapp); //CHKERRQ(info);

  // Free PETSc data structures 
  info = VecDestroy(x); //CHKERRQ(info);
  info = MatDestroy(H); //CHKERRQ(info);

  // Finalize TAO 
  TaoFinalize();
  PetscFinalize();	
}
*/
