#include <iostream>
#include "membrane_system.H"


#undef __FUNCT__
#define __FUNCT__ "main"
double find_big_doub(double p1, double p2){ return p1 > p2 ? p1 : p2; }

int main(int argc,char **argv){
	
	
	if(argc != 3){
		std::cout<<"Usage: "<<argv[0]<<" [parameter file] [output directory]\n";
		exit(0);
	}
	//we have just one membrane	
	MembraneSystem ms(1);
	
	//READ PARAMETERS
	//Reader parameters("parameters.in", "data"); //(parameter_file, folder to copy parameter_file) 
	Reader parameters(argv[1], argv[2]); //(parameter_file, folder to copy parameter_file) 
	ms.read_parameters(&parameters);
	
	//READ THE MESH
	//std::cout<<"Mesh_file"<<parameters.get_char("mesh_file")<<std::endl;
	ms.read_rectangular_off(parameters.get_char("mesh_file") );

	//INITIALIZE THE MESH and FIND BONDS
	ms.construct_from_rect();
	ms.test();

	ms.reset_bond_lengths(); //and angles
	
	//INITIAL CONDITIONS
	if(parameters.get_bool("read_from_mdf"))
		ms.input_solution( parameters.get_char("input_file") );	
	
	if(parameters.get_bool("add_noise"))
		ms.add_noise(); // to have non planar bonds
	if(parameters.get_bool("apply_ic"))
		ms.initial_conditions();

	//ms.initial_torsion();
	//ms.output_solution();
	
	//zeroth output
	ms.clear_forces();	
	ms.update_parameters(0);
	ms.calculate_rectangle_normals();
	ms.calculate_angles();
	ms.cooperativity();
	ms.internal_forces();
	//ms.landau_forces();
	ms.external_forces();
	ms.write_measurements(0,0);
    
    ms.collect_asymmetry_from_diagonals();
	ms.vtk_rect_writer(0);
	
	//ms.clear_forces();

	//ms.vtk_writer(1);
	//ms.output_solution();

	
	//initialize tao
	//ms.tao_construct(&argc, &argv);
//	let tao do the job
	//ms.tao_solve();
	//ms.tao_destroy();
	
	//SWITCH FOR BROWNIAN UPDATE TYPE
	bool is_brownian = parameters.get_bool("is_brownian");
	bool is_clamped = parameters.get_bool("is_clamped");
	int print_each = parameters.get_int("print_each");
	int pol_print_each = parameters.get_int("print_polarization_each");
	
	double coop_count = 1.0;
	int timespan  = parameters.get_int("number_of_steps");

	
	for(int it=1; it<timespan+1; ++it){
		
		ms.update_parameters(it);
		
		if(is_brownian)
			ms.constrained_brownian_update(it);
		else
			ms.velocity_verlet_position(it);
		//if(is_clamped) ms.constrained_brownian_update(it);
		//else ms.brownian_update(it);
		
		ms.clear_forces();
		ms.calculate_rectangle_normals();
		ms.calculate_angles();
		ms.cooperativity(); //collect the dimer angles from neighbors 		
		ms.internal_forces();
		//ms.landau_forces();
		ms.external_forces();

		if(!is_brownian){
			//ms.damping();
			ms.velocity_verlet_velocity(it);
		}
		
		//OUTPUTS
		if(it%pol_print_each ==0){
			ms.write_tube_data(it);
			//ms.write_bond_data(it);			
			}
		
		if(it%print_each == 0){
			ms.collect_asymmetry_from_diagonals();
			std::cout<<"writing the result for step "<<it<<std::endl;
			ms.vtk_rect_writer(it);
			ms.output_solution();
			ms.write_measurements(it, 0);
		}

	}

	
	//ms.tao_construct(&argc, &argv);
	//ms.tao_solve();
	//ms.tao_destroy();
	
	//OUTPUT THE RESULT MESH
	//ms.write_off();	
	//ms.output_solution();
	
	return 0;
}
	




