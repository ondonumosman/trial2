#include <iostream>
//#include "md_vector.h"
#include <math.h>
#include "utility/algebra_support.h"
#include <vector>
 
class Node{
	friend class Membrane;
	friend class MembraneSystem;
	
	public:
		Node(){	
			mass = 1.0;	
			energy = 0.0;
			fspin = 0.0;
			p_fspin = 0.0; 
			vel_spin = 0.0;
			total_energy=0.0;};

		
	RealVectorValue& get_xyz(){ return position;}
	RealVectorValue& get_xyz0(){ return position0;}
	//void print_xyz(){std::cout<<&position<<std::endl;}	
	RealVectorValue& get_force(){return force;}
	RealVectorValue& get_velocity(){return velocity;}
	
	void set_mass(double ms){mass = ms;} 
	
	protected:
		RealVectorValue position;
		RealVectorValue position0;
		RealVectorValue p_position;
		RealVectorValue velocity;
		RealVectorValue force;
		RealVectorValue p_force;
		RealVectorValue velocity_half;
		
		double spin;
		double p_spin;
		double fspin;
		double p_fspin;
		double vel_spin;
		
		double energy;
		double total_energy;
		
		double mass;
		unsigned int index;
	
};

class Bond_old{
	friend class Membrane;
	friend class MembraneSystem;
	
	public:
		Bond_old(){};
		
		Bond_old(unsigned int bi){
			this->bond_index = bi;
			elements[1] = -1; //-1 means empty for now
			lnodes[1] = -1; //-1 means empty for now
			bond_type = 0; //typeless
			coop_index[0]=-1; coop_index[1]=-1; 
			
			bCoef = 0.0; sCoef = 0.0;
			benOn = false; strOn = false;
		};
		
		int get_node(int gt){
			return nodes[gt];
		}
		int get_lnode(int gt){
			return lnodes[gt];
		}
		
		double get_target_angle(){ return sin_angle;}
		double get_length(){ return length;}
		
		//bool is_boundary(){
			//if(this->elements[1]==-1) return true;
			//else return false;
		//}
		//For new diagonals, i modified boundary check  
		//but THIS IS SO DIRTY. 
		bool is_boundary(){
			if(bond_type == 5) return false;
			else{
				if(this->elements[1]==-1) return true;
				else return false;
			}
		}

		
		void swap(){ //be careful with -1 indices
			//int n0 = nodes[0];
			//nodes[0] = nodes[1];
			//nodes[1] = n0;
			
			int n0 = lnodes[0];
			lnodes[0] = lnodes[1];
			lnodes[1] = n0;
			
			n0 = elements[0];
			elements[0] = elements[1];
			elements[1] = n0;
		}

	protected:
		unsigned int nodes[2];
		
		int lnodes[2];  //nodes of triangles that are not bond nodes
		int elements[2];
		
		double length, length1, length2;
		double sin_angle; //target angle for quadratic potential
		int bond_index;
		int index_in_vector;
		unsigned int bond_type;
		
		double bCoef, sCoef;
		bool benOn, strOn;
		
		//all vectors are here now
		RealVectorValue u1, u2, u3, u2n, n_alpha, n_beta, cross_ab;
		double na_size, nb_size, u2size;
		
		//COOP
		double current_sin; //to pass to the neighbors
		int coop_index[2];	//index: (protofilament, dimer)
		double mu_coop;	//cooperativity strength
		double coop_sin;	//constructed from current_sin of neighbors
		double coop_residual;	//constructed from current_sin of neighbors
		//Bond* next, prev;
};


class Bond{
	friend class Membrane;
	friend class MembraneSystem;
	
	public:
		Bond(){};
		
		Bond(unsigned int bi){
			this->bond_index = bi;
			rectangles[1] = -1; //-1 means empty for now
			//lnodes[1] = -1; //-1 means empty for now
			bond_type = 0; //typeless
			coop_index[0]=-1; coop_index[1]=-1; 
			
			bCoef = 0.0; sCoef = 0.0;
			benOn = false; strOn = false;
		};
		
		int get_node(int gt){
			return nodes[gt];
		}
		int get_lnode(int gt){
			return lnodes[gt];
		}		
		
		double get_target_angle(){ return sin_angle;}
		double get_length(){ return length;}
		
		//bool is_boundary(){
			//if(this->elements[1]==-1) return true;
			//else return false;
		//}
		//For new diagonals, i modified boundary check  
		//but THIS IS SO DIRTY. 
		bool is_boundary(){
				if(this->rectangles[0]==-1 || this->rectangles[1]==-1) return true;
				else return false;
		//	}
		}

		
		void swap(){ //be careful with -1 indices
			
			int n0 = rectangles[0];
			rectangles[0] = rectangles[1];
			rectangles[1] = n0;
		}
		
		void swap_nodes(){ //be careful with -1 indices
			int n0 = nodes[0];
			nodes[0] = nodes[1];
			nodes[1] = n0;
		}

	protected:
		int nodes[2];
		
		int lnodes[2];  //nodes of triangles that are not bond nodes
		int elements[2];
		int rectangles[2];
		
		double length, length1, length2;
		double sin_angle; //target angle for quadratic potential
		int bond_index;
		int index_in_vector;
		unsigned int bond_type;
		
		double bCoef, sCoef;
		bool benOn, strOn;
		
		//all vectors are here now
		RealVectorValue u1, u2, u3, u2n, n_alpha, n_beta, cross_ab;
		double na_size, nb_size, u2size;
		
		RealVectorValue E_del_alpha, E_del_beta,  E_del_u2;
		RealTensorValue alpha_del_1, alpha_del_2, beta_del_2, beta_del_3;
		void bending_tensors(double diff, RealVectorValue n_alpha, RealVectorValue n_beta,
				RealVectorValue u1, RealVectorValue u3);

		void bending_tensors_for_anharmonic(double pX, double pCoef, RealVectorValue n_alpha, RealVectorValue n_beta,
				RealVectorValue u1, RealVectorValue u3);
				
		RealVectorValue X_del_alpha, X_del_beta,  X_del_u2;
		
		RealVectorValue bforce1, bforce2, bforce3, bforce4;
		
		//COOP
		//to erase
		double current_sin, coop_sin, coop_residual;
		
		double current_sin1; //to pass to the neighbors
		double current_sin2; //to pass to the neighbors
		int coop_index[2];	//index: (protofilament, dimer)
		double mu_coop;	//cooperativity strength
		double coop_sin1;	//constructed from current_sin of neighbors
		double coop_sin2;	//constructed from current_sin of neighbors
		double coop_residual1;	//constructed from current_sin of neighbors
		double coop_residual2;	//constructed from current_sin of neighbors
		//Bond* next, prev;
};

class Triangle{
	public:
		Triangle(){
		};
		
		void set(int p1, int p2, int p3){
			nodes[0] = p1;
			nodes[1] = p2;
			nodes[2] = p3;
		}
		
		int get(int gt){
			return nodes[gt];
		}
		
		void print(){
			std::cout<<nodes[0]<<" "<<nodes[1]<<" "<<nodes[2]<<std::endl;
		}
	protected:
		unsigned int nodes[3];
};

class Rectangle{
	friend class Membrane;
	friend class MembraneSystem;	
	public:
		Rectangle(){
			asymmetry = 0.0;
		};
		
		void set(unsigned int id, int p1, int p2, int p3, int p4){
			nodes[0] = p1;
			nodes[1] = p2;
			nodes[2] = p3;
			nodes[3] = p4;
			idx = id;
		}
		
		int get(int gt){
			return nodes[gt];
		}
		
		void print(){
			std::cout<<nodes[0]<<" "<<nodes[1]<<" "<<nodes[2]<<" "<<nodes[3]<<std::endl;
		}
		
		
	protected:
		unsigned int nodes[4];
		RealVectorValue normals[4];
		unsigned int idx;
		double asymmetry;
};
 
class Membrane{
	public:
		Membrane(){	};
		Membrane(int n_nod, int n_rec){
			n_nodes = n_nod; n_rectangles = n_rec; n_triangles = 4*n_rec;
			nodes = new Node[n_nodes];
			rectangles = new Rectangle[n_rectangles];
			triangles = new Triangle[n_triangles];
		};
		
		void resize(int n_nod, int n_rec){
			n_nodes = n_nod; n_rectangles = n_rec; n_triangles = 4*n_rec;
			nodes = new Node[n_nodes];
			rectangles = new Rectangle[n_rectangles];
			triangles = new Triangle[n_triangles];
		};
		
		Node& get_node(int nn){return nodes[nn];}
		int get_n_nodes(){return n_nodes;}
		
		Triangle& get_triangle(int nf){return triangles[nf];}
		Rectangle& get_rectangle(int nf){return rectangles[nf];}
		int get_n_triangles(){return n_triangles;}
		int get_n_rectangles(){return n_rectangles;}
		
		void check_rectangle_order();
	
		Bond& get_bond(int nf){return bonds[nf];}
		int get_n_bonds(){return n_bonds;}
		
		void find_bonds(double a_1, double a_2, double a_3, double lngth, double lngth1, double lngth2, bool seam);
		void check_orientations(int mesh_type);
		void alternate(double cc);
		void alternate_seam(double cc);
		void two_diagonals();
		
		
		void construct(double cc);
		void initial_bonds();
		void add_diagonals();
		void print_rect();			
		
		void print();		
	
		
		Bond* get_coop_bond(int pp, int dd){return cooperativity[pp][dd];}
		int get_coop_bond2(int pp, int dd){return cooperativity2[pp][dd];}
		int get_n_coop(){return n_coop;}
	protected:
		Node* nodes;
		Triangle* triangles;
		Rectangle* rectangles;
		std::vector<Bond> bonds;
		std::vector<Bond>::iterator bond_it;
		
		int n_nodes, n_triangles, n_rectangles, n_bonds;
		Bond ***cooperativity; //2D array of pointers to bonds
		int **cooperativity2;
		int n_coop; //how many dimer along a protofilament

};
