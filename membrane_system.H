#include "membrane.H"
#include <fstream>
#include <vector>
#include <iomanip> //to set output precision
//#include "tao.h"
#include "utility/ran.h"
#include "utility/reader.h"

template <class T>
inline std::string to_string (const T& t)
{
std::stringstream ss;
ss << t;
return ss.str();
}


class MembraneSystem;
/*
//CONTEXT FOR TAO
typedef struct {
	PetscInt n;          // dimension 
	PetscReal alpha;   // condition parameter 
	MembraneSystem *sys;
} AppCtx;
*/
class MembraneSystem{
	public:
		MembraneSystem();
		MembraneSystem(int n_memb);
		
		void read_off(const char* off_file);
		void read_rectangular_off(const char* off_file);
		void write_off();
		void vtk_writer(int step);
		void vtk_rect_writer(int step);
		void write_measurements(int nit, int every);
		void construct();
		void construct_from_rect();
		void test();
		void write_tube_data(int step);
		void write_bond_data(int step);
		void collect_asymmetry_from_diagonals();
		
		void read_parameters(Reader* params);
		void input_solution(const char* in_file);
		void output_solution();
		
		//tao
		/*
		void tao_construct(int* argc, char ***argv);
		void tao_solve();
		void tao_destroy();

		void push_solution();
		void pull_solution();
		void pull_solution_sol(PetscScalar *sl[]);
		
		void function_gradient(PetscScalar *sol[], double *ff, PetscScalar *grad[]);
		*/
		
		void update_parameters(int nit);
		
		//forces
		void internal_forces();
		void clear_forces();
		
		void landau_forces();
		
		void external_forces();//wrapper for all external forces
		void nodal_forces();
		void deflection_force();
		void nodal_moment();
		void container_forces();
		void set_radius(double rad){cradius = rad;};
		
		void cooperativity();
		void calculate_normals_and_angles();
		void update_coop_value(double fac);
		
		void calculate_rectangle_normals();
		void calculate_angles();		
		
		void damping();
		double gamma;
		double mass;
		double boltzmann;
		double noise_coef;
		double diffusion;
		double diffusion_coef;
		double brownian_coef;
	
		void set_temperature(double ttemp);
	
		void add_noise();
		void initial_conditions();
		void reset_bond_lengths();
		void initial_torsion();
		
	
		//integrator
		void velocity_verlet_position(int nit);
		void velocity_verlet_velocity(int nit);		
		
		void brownian_update(int nit);
		void constrained_brownian_update(int nit);
		
	protected:
		Membrane sample;
		Membrane* membranes;
		
		int n_membranes;
		int n_nodes;
		int n_triangles;
		int n_rectangles;
		int n_bonds;
	
		Ran* r1;
		
		double dt;


		double wallX, wallY;
		double kinetic_energy;
		double potential_energy;	
		double temperature;
		std::ofstream measurements; 
		
		double Estretch;
		double Ebending;
		double Ebend4;
		double Ecoop;
		double Etotal;
		
		bool stretching_on, bending_on;
		double Srigidity1,Srigidity2,Srigidity3,Srigidity4, Srigidity5;
		double Brigidity1, Brigidity2, Brigidity3, Brigidity4, Brigidity5;
		double d_bond, d_bond1, d_bond2; 
		double alpha_1, alpha_2, alpha_3, alpha_4, alpha_5;
		int mesh_type;
		bool has_seam;
		const char* mesh_file;
		
		//quartic potential parameters, heaviside value for polarization
		double pA, pB, pC, limit_angle;
		bool harmonic;
		
		//cooperativity
		double coop_coef;
		
		//spin
		double Elandau;
		double Einter;
		bool spin_on;
		double deltaE, JJ;
		
		
		double Eextern;
		//nodal forces
		bool nodal_on;
		double mu_nodal;
		double max_force_strength, force_strength, update_coef;
		bool force_is_gradual;
		double Enodal;
		double boundary0, boundary1, boundary2; //limits for external forces or clamping
		bool is_clamped;
		
		//moment
		double moment_boundary_up, moment_boundary_down;
		double max_moment_strength, moment_strength;
		bool moment_on, moment_is_gradual;
		
		int printeach;
		
		//container
		double mu_cont;
		bool container_on;
		double cradius;
		double Econtainer;
		
		double noise_strength;
		
		const char* directory;
	//to use in push and pull solution functions
	//PetscInt* ix_kds;
	//PetscScalar* y_kds;
	
	//int        info;                  /* used to check for functions returning nonzeros */
	//PetscScalar zero;
	//Vec        x;                     /* solution vector */
	//Mat        H;                     /* Hessian matrix */
	//TAO_SOLVER tao;                   /* TAO_SOLVER solver context */
	//TAO_APPLICATION taoapp;           /* TAO application context */
	//PetscTruth  flg;
	//int        size,rank;                  /* number of processes running */
	//TaoTerminateReason reason;
	//AppCtx     user;                  /* user-defined application context */
	
	int tao_iter, tao_iter0, vtk_iter;
	double fdummy, gnormm, cnormm, xdifff;		
		
};
