#include "membrane.H"

int find_big(int p1, int p2){ return p1 > p2 ? p1 : p2; }
int find_small(int p1, int p2){ return p1 < p2 ? p1 : p2; } //do it with template

double signof(double px){ return (px > 0.0) - (px < 0.0); }

int index_gen(int p1, int p2){
	//order the node index
	int small = find_small(p1,p2);
	int big = find_big(p1,p2);
	
	//bond index for (5,3) is 30005
	return small*100000+big;
}

int protofilament_index(int np, RealVectorValue mid){
	//np: number of protofilaments in the tube
	//mid: the midpoint vector of a horizontal bond
	//angle is the polar angle of the midpoint of a horizontal bond
	
	double xcomp = mid(0), ycomp = -mid(1); //minus because, we move clockwise
	double smid = mid.size();
	mid /= smid;
	double angl = (180.0/PI)*atan(ycomp/xcomp); //CAREFUL, works only if x!=0
	if(ycomp<0.0){
		if(xcomp>0.0) angl += 360.0;
		else angl += 180.0;
	}
	if(ycomp>=0.0 && xcomp<0.0) angl +=180.0;
	
	double dth = 360.0/np;
	
	angl += dth*0.5 ; //+1 because integer cast rounds the number down
	double idx = angl/dth ;
	//std::cout<<"p index: "<<idx<<" "<< angl<<std::endl;//idx<<" "<<int(idx)<<std::endl; 
	return int(idx+0.01);
}

void Bond::bending_tensors(double diff, RealVectorValue n_alpha, RealVectorValue n_beta,
RealVectorValue u1, RealVectorValue u3){
	
	RealVectorValue cross_ab = n_alpha.cross(n_beta);
	double na_size = n_alpha.size();	
	double nb_size = n_beta.size();	
	
	E_del_alpha = this->bCoef*diff*(n_beta.cross(this->u2n));
	E_del_beta = this->bCoef*diff*(this->u2n.cross(n_alpha));
	E_del_u2 = this->bCoef*diff*(cross_ab - (cross_ab*this->u2n)*this->u2n)/this->u2size;
	
	alpha_del_1.zero();
	alpha_del_2.zero();
	beta_del_2.zero();
	beta_del_3.zero();

	for(int ii = 0; ii<3; ii++){
		for(int jj=0; jj<3; jj++){
			//std::cout<<"index "<<ii<<" "<<jj<<std::endl;
			for(int kk=0; kk<3; kk++){
				alpha_del_1(ii,jj) += this->u2(kk)*(le_ci[jj][ii][kk] 
						- n_alpha(jj)*n_alpha(0)*le_ci[0][ii][kk] 
						- n_alpha(jj)*n_alpha(1)*le_ci[1][ii][kk] 
						- n_alpha(jj)*n_alpha(2)*le_ci[2][ii][kk] ); 


				alpha_del_2(ii,jj) += u1(kk)*(le_ci[jj][kk][ii] 
						- n_alpha(jj)*n_alpha(0)*le_ci[0][kk][ii] 
						- n_alpha(jj)*n_alpha(1)*le_ci[1][kk][ii] 
						- n_alpha(jj)*n_alpha(2)*le_ci[2][kk][ii] ); 						
		
		
				beta_del_2(ii,jj) += u3(kk)*(le_ci[jj][ii][kk] 
						- n_beta(jj)*n_beta(0)*le_ci[0][ii][kk] 
						- n_beta(jj)*n_beta(1)*le_ci[1][ii][kk] 
						- n_beta(jj)*n_beta(2)*le_ci[2][ii][kk] ); 
				

				beta_del_3(ii,jj) += this->u2(kk)*(le_ci[jj][kk][ii] 
						- n_beta(jj)*n_beta(0)*le_ci[0][kk][ii] 
						- n_beta(jj)*n_beta(1)*le_ci[1][kk][ii] 
						- n_beta(jj)*n_beta(2)*le_ci[2][kk][ii] ); 									
				
			}
			alpha_del_1(ii,jj) /= na_size; 
			alpha_del_2(ii,jj) /= na_size; 
			beta_del_2(ii,jj) /= nb_size; 
			beta_del_3(ii,jj) /= nb_size; 
		}
	}	

	bforce3 = alpha_del_1*E_del_alpha;
						
	bforce4 = beta_del_3*E_del_beta;
	
	bforce1 = -alpha_del_1*E_del_alpha - alpha_del_2*E_del_alpha
							-beta_del_2*E_del_beta - beta_del_3*E_del_beta - E_del_u2;
				
	bforce2 = alpha_del_2*E_del_alpha + beta_del_2*E_del_beta + E_del_u2;		

}

void Bond::bending_tensors_for_anharmonic(double pX, double pCoef, RealVectorValue n_alpha, RealVectorValue n_beta,
RealVectorValue u1, RealVectorValue u3){
	
	RealVectorValue cross_ab = n_alpha.cross(n_beta);
	double na_size = n_alpha.size();	
	double nb_size = n_beta.size();	
	
	X_del_alpha = (n_beta.cross(this->u2n));
	X_del_beta = (this->u2n.cross(n_alpha));
	X_del_u2 = (cross_ab - pX*this->u2n)/this->u2size;	
	
	alpha_del_1.zero();
	alpha_del_2.zero();
	beta_del_2.zero();
	beta_del_3.zero();

	for(int ii = 0; ii<3; ii++){
		for(int jj=0; jj<3; jj++){
			//std::cout<<"index "<<ii<<" "<<jj<<std::endl;
			for(int kk=0; kk<3; kk++){
				alpha_del_1(ii,jj) += this->u2(kk)*(le_ci[jj][ii][kk] 
						- n_alpha(jj)*n_alpha(0)*le_ci[0][ii][kk] 
						- n_alpha(jj)*n_alpha(1)*le_ci[1][ii][kk] 
						- n_alpha(jj)*n_alpha(2)*le_ci[2][ii][kk] ); 


				alpha_del_2(ii,jj) += u1(kk)*(le_ci[jj][kk][ii] 
						- n_alpha(jj)*n_alpha(0)*le_ci[0][kk][ii] 
						- n_alpha(jj)*n_alpha(1)*le_ci[1][kk][ii] 
						- n_alpha(jj)*n_alpha(2)*le_ci[2][kk][ii] ); 						
		
		
				beta_del_2(ii,jj) += u3(kk)*(le_ci[jj][ii][kk] 
						- n_beta(jj)*n_beta(0)*le_ci[0][ii][kk] 
						- n_beta(jj)*n_beta(1)*le_ci[1][ii][kk] 
						- n_beta(jj)*n_beta(2)*le_ci[2][ii][kk] ); 
				

				beta_del_3(ii,jj) += this->u2(kk)*(le_ci[jj][kk][ii] 
						- n_beta(jj)*n_beta(0)*le_ci[0][kk][ii] 
						- n_beta(jj)*n_beta(1)*le_ci[1][kk][ii] 
						- n_beta(jj)*n_beta(2)*le_ci[2][kk][ii] ); 									
				
			}
			alpha_del_1(ii,jj) /= na_size; 
			alpha_del_2(ii,jj) /= na_size; 
			beta_del_2(ii,jj) /= nb_size; 
			beta_del_3(ii,jj) /= nb_size; 
		}
	}	

	bforce3 = pCoef*alpha_del_1*X_del_alpha;
						
	bforce4 = pCoef*beta_del_3*X_del_beta;
	
	bforce1 = pCoef*(-alpha_del_1*X_del_alpha - alpha_del_2*X_del_alpha
							-beta_del_2*X_del_beta - beta_del_3*X_del_beta -X_del_u2);
				
	bforce2 = pCoef*(alpha_del_2*X_del_alpha + beta_del_2*X_del_beta + X_del_u2);		

}

void Membrane::find_bonds(double a_1, double a_2, double a_3, double lngth, double lngth1, double lngth2, bool seam){
	double lower_limit = -0.001, upper_limit = 0.001;
	if(seam){lower_limit = 0.20; upper_limit = 0.25;}
	
	//the vector needed to compare if a bond
	//is lateral or longitudinal
	RealVectorValue xy(1.0,1.0,0.0);
	xy /= xy.size();
	RealVectorValue zz(0.0,0.0,1.0);
	int bond_counter = 0;
	for (int eit = 0; eit<n_triangles; eit++){ //iterate over elements
		//std::cout<<"element no: "<<eit<<std::endl;
		Triangle& e1 = triangles[eit];
		
		//std::cout<<e1.get(0)<<" "<<e1.get(1)<<" "<<e1.get(2)<<std::endl;
		
		for(int lit=0; lit<3; lit++){ //iterate over edges
				
			//node indices
			int nd1=e1.get(lit), nd2=e1.get((lit+1)%3); //get 0,1 then 1,2 then 2,0
			int nd3=e1.get((lit+2)%3); //get the other third node
			
			//generate corresponding bond index for this edge
			int dummy_idx = index_gen(nd1, nd2);
			bool newbond = true;
			
			//iterate over the existing bond list 
			for ( bond_it=bonds.begin() ; bond_it < bonds.end(); bond_it++ ){
					//if we have already this bond in the list
					// just add current element index to the element list					
					if (bond_it->bond_index == dummy_idx){
						//std::cout<<"inside\n";
						bond_it->elements[1] = eit;
						bond_it->lnodes[1] =nd3; //and add elements 3rd node
						newbond = false;
					}
					//else, do nothing
			}
				
			if(newbond){ //if we could not find it in the list
				Bond b1(dummy_idx);
				b1.elements[0] = eit; //add this triangle 
				b1.nodes[0] = find_small(nd1,nd2);
				b1.nodes[1] = find_big(nd1,nd2);
				b1.lnodes[0] = nd3;
				b1.length = lngth; b1.length1 = lngth1; b1.length2 = lngth2;
				
				//get a reference to the particle
				Node& p1 = nodes[nd1];
				Node& p2 = nodes[nd2];
				
				
				RealVectorValue orientation = p1.position - p2.position;
				orientation /= orientation.size();
				double sproduct = fabs(orientation*zz);
				//std::cout<<"debug: "<<nd1<<" "<<nd2<<" "<<sproduct<<std::endl;				
				
				//if(sproduct<0.001 and sproduct>-0.001) { //HORIZONTAL regular
				//if(sproduct<0.45 and sproduct>0.35) { //HORIZONTAL
				if(sproduct<upper_limit and sproduct>lower_limit) { //HORIZONTAL seam
					b1.bond_type = 1;
					b1.sin_angle = a_1;
				}
				else if(sproduct<1.001 and sproduct>0.999) { //VERTICAL
					b1.bond_type = 3;
					b1.sin_angle = a_3;
				}
				else{						//DIAGONAL
					b1.bond_type = 4;
					b1.sin_angle = a_2;
				}
				b1.index_in_vector = bond_counter;
				bond_counter++;
				bonds.push_back(b1);	
			}
		}
	}
	
	this->n_bonds = bonds.size();
	std::cout<<"bond done "<<bonds.size()<<"\n";
	//debug bonds
	
	//this->print();
} 

void Membrane::check_orientations(int mesh_type){
	//CHECK NORMAL ORIENTATIONS
	RealVectorValue compare;
	
	for ( bond_it=bonds.begin() ; bond_it < bonds.end(); bond_it++ ){
		
		Bond& b1 = *bond_it;
		int nd1 = b1.get_node(0);
		int nd2 = b1.get_node(1);
		//get a reference to primary nodes
		Node& p1 =this->get_node(nd1);
		Node& p2 = this->get_node(nd2);
		
		RealVectorValue u2 = p2.get_xyz() - p1.get_xyz();

		//CHECK INTERNAL BONDS, NOT BOUNDARY
		if(!b1.is_boundary()){
			
			//now we also need reference to secondary nodes
			int nd3 = b1.get_lnode(0);
			int nd4 = b1.get_lnode(1);
			Node& p3 = this->get_node(nd3);
			Node& p4 = this->get_node(nd4);

			RealVectorValue u1 = p3.get_xyz() - p1.get_xyz();
			RealVectorValue u3 = p4.get_xyz() - p1.get_xyz();

			RealVectorValue n_alpha = u1.cross(u2);
			double na_size = n_alpha.size();
			n_alpha /= na_size;
			RealVectorValue n_beta = u2.cross(u3);
			double nb_size = n_beta.size();
			n_beta /= nb_size;
	

			//RealVectorValue cross_ab = n_alpha.cross(n_beta);
			//double cross_size = cross_ab.size();
			
			compare.zero();
			if(mesh_type == 0) //for plane case
				compare(1) = 1.0;  //so it is (0,1,0)
			if(mesh_type ==1) //cylinder case
				compare = p1.get_xyz0();
			
		//	std::cout<<nd1<<" "<<nd2<<" n_alpha: "<<n_alpha<<std::endl;
		//	std::cout<<nd3<<" "<<nd4<<" n_beta: "<<n_beta<<std::endl;
		//	std::cout<<nd1<<" "<<nd2<<" comp: "<<n_alpha*compare<<std::endl;
		//	std::cout<<std::endl;
			
			if(compare*n_alpha>0){
				if(compare*n_beta<0){
					std::cout<<"\nCHECK ORIENTATIONS\n";
				}
				else{
				//	std::cout<<"swapping\n";
					b1.swap();
				}
			}
		}
		
		
	}
	//this->print();

}



void Membrane::print(){
	std::cout<<"BOND DEBUG START\n";
	for ( bond_it=bonds.begin() ; bond_it < bonds.end(); bond_it++ ){

		std::cout<<"bond vector index: "<<bond_it->index_in_vector<<std::endl;
		std::cout<<"bond index: "<<bond_it->bond_index<<std::endl;
		std::cout<<"bond nodes: "<<bond_it->nodes[0]<<" "<<bond_it->nodes[1]<<std::endl;
		std::cout<<"bond 2nd nodes: "<<bond_it->lnodes[0]<<" "<<bond_it->lnodes[1]<<std::endl;

			Bond& b1 = *bond_it;

			if(!b1.is_boundary()){
				
				int nd1 = b1.get_node(0);
				int nd2 = b1.get_node(1);
			
				//get a reference to the particle
				Node& p1 =this->get_node(nd1);
				Node& p2 = this->get_node(nd2);
				
				//now we also need reference to secondary nodes
				int nd3 = b1.get_lnode(0);
				int nd4 = b1.get_lnode(1);
			
				Node& p3 = this->get_node(nd3);
				Node& p4 = this->get_node(nd4);
				RealVectorValue u1 = p3.get_xyz() - p1.get_xyz();
				RealVectorValue u2 = p2.get_xyz() - p1.get_xyz();
				RealVectorValue u3 = p4.get_xyz() - p1.get_xyz();

				RealVectorValue n_alpha = u1.cross(u2);
				double na_size = n_alpha.size();
				n_alpha /= na_size;
				RealVectorValue n_beta = u2.cross(u3);
				double nb_size = n_beta.size();
				n_beta /= nb_size;
	
				RealVectorValue cross_ab = n_alpha.cross(n_beta);
				double cross_size = cross_ab.size();

				std::cout<<"angle: "<<cross_size<<std::endl;
				std::cout<<"orient a: "<<p1.get_xyz()*n_alpha<<std::endl;
				std::cout<<"orient b: "<<p1.get_xyz()*n_beta<<std::endl;

		}

		std::cout<<"bond elements: "<<bond_it->elements[0]<<" "<<bond_it->elements[1]<<std::endl;
		std::cout<<"bond type: "<<bond_it->bond_type<<std::endl;
		std::cout<<"bond is on boundary?: "<<bond_it->is_boundary()<<std::endl;
		std::cout<<"bond angle: "<<bond_it->sin_angle<<std::endl;
		std::cout<<"bond length: "<<bond_it->length<<std::endl;
		std::cout<<std::endl;
		
	}
	std::cout<<"BOND DEBUG END\n";	
}

void Membrane::alternate(double cc){
	std::cout<<"BOND DEBUG START\n";
	int counter = 0;
	
	for ( bond_it=bonds.begin() ; bond_it < bonds.end(); bond_it++ ){

	//if(!b1.is_boundary() && bond_it->bond_type == 1){
	if(bond_it->bond_type == 1){
		//std::cout<<"bond index: "<<bond_it->bond_index<<std::endl;
		//std::cout<<"bond index in vector: "<<bond_it->index_in_vector<<std::endl;		
		//std::cout<<"bond nodes: "<<bond_it->nodes[0]<<" "<<bond_it->nodes[1]<<std::endl;

		Bond& b1 = *bond_it;
					
		int nd1 = b1.get_node(0);
		int nd2 = b1.get_node(1);
			
		//get a reference to the particle
		Node& p1 =this->get_node(nd1);
		Node& p2 = this->get_node(nd2);
				
		RealVectorValue middle = 0.5*(p2.get_xyz() + p1.get_xyz());
		int indx = protofilament_index(13, middle) - 1;
		//double shift = 3.0/13.0;
		double zcomp = middle(2);
		//std::cout<<"z: "<<indx<<" "<<zcomp<<std::endl;
		int height_index = int(zcomp+0.01);
		if(height_index % 2){  //if it is odd
			b1.bond_type = 2;
			b1.coop_index[0] = indx;
			b1.coop_index[1] = (height_index-1)/2;
			//std::cout<<(height_index-1)/2<<std::endl;
			counter++;
		} 
		
		//std::cout<<"angle: "<<anglee<<std::endl;
//		std::cout<<"bond elements: "<<bond_it->elements[0]<<" "<<bond_it->elements[1]<<std::endl;
//		std::cout<<"bond type: "<<b1.bond_type<<std::endl;
//		std::cout<<"bond coop: "<<b1.coop_index[0]<<" "<<b1.coop_index[1]<<std::endl<<std::endl;
		
		//std::cout<<"bond is on boundary?: "<<bond_it->is_boundary()<<std::endl;
		//std::cout<<"bond angle: "<<bond_it->sin_angle<<std::endl;
		//std::cout<<"bond length: "<<bond_it->length<<std::endl;
		//std::cout<<std::endl;
		}
	}
	
	//FOR 13 PROTOFILAMENT
	n_coop = counter/13;
	cooperativity = new Bond**[13];
	cooperativity2 = new int*[13];
	for (int i = 0; i < 13; ++i) {
		cooperativity[i] = new Bond*[n_coop];
		cooperativity2[i] = new int[n_coop];
	}
	
	for ( bond_it=bonds.begin() ; bond_it < bonds.end(); bond_it++ ){
		Bond& b1 = *bond_it;
		if(b1.bond_type==2){
			
			cooperativity[b1.coop_index[0]][b1.coop_index[1]] = &b1;
			cooperativity2[b1.coop_index[0]][b1.coop_index[1]] = b1.index_in_vector;
			if(b1.coop_index[1] ==0 || b1.coop_index[1]==n_coop-1)
				b1.mu_coop = cc;
			else
				b1.mu_coop = 2.0*cc;
			//std::cout<<"Index: "<<b1.bond_index<<" "<<b1.coop_index[0]<<" "<<b1.coop_index[1]<<" "<<b1.mu_coop<<std::endl;
		}
	}
	
	//for ( bond_it=bonds.begin() ; bond_it < bonds.end(); bond_it++ ){
		//Bond& b1 = *bond_it;
		//if(b1.bond_type==2){
			//if(b1.coop_index[1] ==0
		//}
	//}		

	
	std::cout<<"BOND DEBUG END\n";	
}

void Membrane::alternate_seam(double cc){
	std::cout<<"BOND DEBUG START\n";
	int counter = 0;
	
	for ( bond_it=bonds.begin() ; bond_it < bonds.end(); bond_it++ ){

	//if(!b1.is_boundary() && bond_it->bond_type == 1){
	if(bond_it->bond_type == 1){
		//std::cout<<"bond index: "<<bond_it->bond_index<<std::endl;
//		std::cout<<"bond nodes: "<<bond_it->nodes[0]<<" "<<bond_it->nodes[1]<<std::endl;

		Bond& b1 = *bond_it;
					
		int nd1 = b1.get_node(0);
		int nd2 = b1.get_node(1);
			
		//get a reference to the particle
		Node& p1 =this->get_node(nd1);
		Node& p2 = this->get_node(nd2);
				
		RealVectorValue middle = 0.5*(p2.get_xyz() + p1.get_xyz());
		int indx = protofilament_index(13, middle) - 1;
		double shift = 3.0/13.0;
		double zcomp = middle(2) - 0.5*shift -indx*shift;
		int height_index = int(zcomp+0.01);
		if(height_index % 2){  //if it is odd
			b1.bond_type =2;
			b1.coop_index[0] = indx;
			b1.coop_index[1] = (height_index-1)/2;
			counter++;
		} 
		
		//std::cout<<"angle: "<<anglee<<std::endl;
		//std::cout<<"bond elements: "<<bond_it->elements[0]<<" "<<bond_it->elements[1]<<std::endl;
	//	std::cout<<"bond type: "<<bond_it->bond_type<<std::endl;
		//std::cout<<"bond is on boundary?: "<<bond_it->is_boundary()<<std::endl;
		//std::cout<<"bond angle: "<<bond_it->sin_angle<<std::endl;
		//std::cout<<"bond length: "<<bond_it->length<<std::endl;
		//std::cout<<std::endl;
		}
	}
	
	//FOR 13 PROTOFILAMENT
	n_coop = counter/13;
	cooperativity = new Bond**[13];
	cooperativity2 = new int*[13];
	for (int i = 0; i < 13; ++i) {
		cooperativity[i] = new Bond*[n_coop];
		cooperativity2[i] = new int[n_coop];
	}
	
	for ( bond_it=bonds.begin() ; bond_it < bonds.end(); bond_it++ ){
		Bond& b1 = *bond_it;
		if(b1.bond_type==2){
			
			cooperativity[b1.coop_index[0]][b1.coop_index[1]] = &b1;
			cooperativity2[b1.coop_index[0]][b1.coop_index[1]] = b1.index_in_vector;
			if(b1.coop_index[1] ==0 || b1.coop_index[1]==n_coop-1)
				b1.mu_coop = cc;
			else
				b1.mu_coop = 2.0*cc;
			//std::cout<<"Index: "<<b1.bond_index<<" "<<b1.coop_index[0]<<" "<<b1.coop_index[1]<<" "<<b1.mu_coop<<std::endl;
		}
	}
	
	//for ( bond_it=bonds.begin() ; bond_it < bonds.end(); bond_it++ ){
		//Bond& b1 = *bond_it;
		//if(b1.bond_type==4){
			//if(b1.coop_index[1] ==0
		//}
	//}		

	
	std::cout<<"BOND DEBUG END\n";	
}

void Membrane::two_diagonals(){
	int counter2 = 0; 
	int node1=-1, node2=-1, lnode1 = -1, lnode2 = -1; 
	int dummy_idx;
	std::vector<Bond> additional_bonds;
	
	for ( bond_it=bonds.begin() ; bond_it < bonds.end(); bond_it++ ){
		Bond& b1 = *bond_it;
		if(b1.bond_type==4) {
			//inverse nodes and lateral nodes
			node1 = b1.lnodes[0];
			node2 = b1.lnodes[1];
			lnode1 = b1.nodes[0];
			lnode2 = b1.nodes[1];
			dummy_idx = index_gen(node1, node2);
			
			Bond new_bond(dummy_idx);
			//triangles are not defined for these bonds
			new_bond.elements[0] = -1;
			new_bond.nodes[0] = find_small(node1,node2);
			new_bond.nodes[1] = find_big(node1,node2);
			
			//BE SURE
			new_bond.lnodes[0] = lnode1;
			new_bond.lnodes[1] = lnode2;
			new_bond.bond_type = 5;
			//we make use of previously calculated n_bonds in find_bonds()
			new_bond.index_in_vector = this->n_bonds + counter2;
			additional_bonds.push_back(new_bond);			
			counter2++;
		}
	}
	
	//Augment the bonds vector
	bonds.insert(bonds.end(),additional_bonds.begin(),additional_bonds.end());

	std::cout<<"2 diagonals: "<<counter2<< std::endl;
	std::cout<<"FOR NOW, the orientation of secondary \
	diagonals is not verified so AVOID to define a Bending forces\n";
	this->n_bonds = bonds.size();
	
	std::cout<<"bond done "<<bonds.size()<<"\n";
	//this->print();	
}



void Membrane::initial_bonds(){
	double lower_limit = -0.001, upper_limit = 0.001;
	//if(seam){lower_limit = 0.20; upper_limit = 0.25;}	
	RealVectorValue zz(0.0,0.0,1.0);

	int bond_counter = 0;
	for (int rit = 0; rit<n_rectangles; rit++){ //iterate over elements
		//std::cout<<"element no: "<<eit<<std::endl;
		Rectangle& r1 = rectangles[rit];
		for(int lit=0; lit<4; lit++){
			
			//node indices
			//std::cout<<lit<<" "<<(lit+1)%4<<std::endl;
			int nd1=r1.get(lit), nd2=r1.get((lit+1)%4); //get 0,1 then 1,2 then 2,3 then 3,0
			
			//generate corresponding bond index for this edge
			int dummy_idx = index_gen(nd1, nd2);
			bool newbond = true;

			//iterate over the existing bond list 
			for ( bond_it=bonds.begin() ; bond_it < bonds.end(); bond_it++ ){
					//if we have already this bond in the list
					// just add current element index to the element list					
					if (bond_it->bond_index == dummy_idx){
						//std::cout<<"inside\n";
						bond_it->rectangles[1] = rit;
					//	bond_it->lnodes[1] =nd3; //and add elements 3rd node
						newbond = false;
					}
					//else, do nothing
			}

			if(newbond){ //if we could not find it in the list
				Bond b1(dummy_idx);
				b1.rectangles[0] = rit; //add this rectangle							
				b1.nodes[0] = find_small(nd1,nd2);
				b1.nodes[1] = find_big(nd1,nd2);
				//get a reference to the particle
				
				Node& p1 = nodes[nd1];
				Node& p2 = nodes[nd2];				

				RealVectorValue orientation = p1.position - p2.position;
				orientation /= orientation.size();
				double sproduct = fabs(orientation*zz);
				
				if(sproduct<upper_limit and sproduct>lower_limit) { //HORIZONTAL seam
					b1.bond_type = 1;
				}
				else if(sproduct<1.001 and sproduct>0.999) { //VERTICAL
					b1.bond_type = 3;
				}
				b1.index_in_vector = bond_counter;
				bond_counter++;
				bonds.push_back(b1);
			}
							
		}
		
	}
	
	this->n_bonds = bonds.size();
	std::cout<<"bond done "<<bonds.size()<<"\n";
}

void Membrane::check_rectangle_order(){
	for ( bond_it=bonds.begin() ; bond_it < bonds.end(); bond_it++ ){
		Bond& b1 = *bond_it;
		if(b1.bond_type==3){
			//std::cout<<"Bond: "<<b1.bond_index<<" "<<b1.bond_type<<std::endl;
			//std::cout<<b1.rectangles[0]<<" "<<b1.rectangles[1]<<std::endl;
			int diff = fabs(b1.rectangles[0]-b1.rectangles[1]);
			if(diff>1) b1.swap();
			//std::cout<<b1.rectangles[0]<<" "<<b1.rectangles[1]<<std::endl;
		}
		if(b1.bond_type==1 && !b1.is_boundary()){
			//std::cout<<"Bond: "<<b1.bond_index<<" "<<b1.bond_type<<std::endl;
			//std::cout<<b1.nodes[0]<<" "<<b1.nodes[1]<<std::endl;
			int diff = fabs(b1.nodes[0]-b1.nodes[1]);			
			if(diff>1) b1.swap_nodes();
			//std::cout<<b1.nodes[0]<<" "<<b1.nodes[1]<<std::endl;
		}
	}
}

void Membrane::add_diagonals(){
	int counter2 = 0; 
	int node1=-1, node2=-1;//, lnode1 = -1, lnode2 = -1; 
	int dummy_idx;
	std::vector<Bond> additional_bonds;
	
	for (int rit = 0; rit<n_rectangles; rit++){ //iterate over elements
		//std::cout<<"element no: "<<eit<<std::endl;
		Rectangle& r1 = rectangles[rit];
		
////////TYPE 4
		node1=r1.get(1); node2=r1.get(3);
		//generate corresponding bond index for this edge
		dummy_idx = index_gen(node1, node2);
		Bond new_bond4(dummy_idx);
		//They are only in this rectangle
		new_bond4.rectangles[0] = rit;				
		new_bond4.rectangles[1] = rit;
		new_bond4.nodes[0] = find_small(node1,node2);
		new_bond4.nodes[1] = find_big(node1,node2);
		new_bond4.bond_type = 4;
		//we make use of previously calculated n_bonds in find_bonds()
		new_bond4.index_in_vector = this->n_bonds + counter2;
		additional_bonds.push_back(new_bond4);			
		counter2++;						
////////TYPE 4 END

////////TYPE 5
		node1=r1.get(0); node2=r1.get(2);
		//generate corresponding bond index for this edge
		dummy_idx = index_gen(node1, node2);
		Bond new_bond5(dummy_idx);
		//They are only in this rectangle
		new_bond5.rectangles[0] = rit;				
		new_bond5.rectangles[1] = rit;
		new_bond5.nodes[0] = find_small(node1,node2);
		new_bond5.nodes[1] = find_big(node1,node2);
		new_bond5.bond_type = 5;
		//we make use of previously calculated n_bonds in find_bonds()
		new_bond5.index_in_vector = this->n_bonds + counter2;
		additional_bonds.push_back(new_bond5);			
		counter2++;						
////////TYPE 5 END

	}

	//Augment the bonds vector
	bonds.insert(bonds.end(),additional_bonds.begin(),additional_bonds.end());

	std::cout<<"Z diagonals: "<<counter2<< std::endl;
	this->n_bonds = bonds.size();
	
	std::cout<<"bond done "<<bonds.size()<<"\n";	

}

void Membrane::construct(double cc){
	this->initial_bonds();
	
	this->check_rectangle_order(); //CAREFUL: WE DON'T HAVE TYPE 2 YET
	this->alternate(cc);
	
	this->add_diagonals();
	//this->print_rect();
}


void Membrane::print_rect(){
	//DEBUG
	for ( bond_it=bonds.begin() ; bond_it < bonds.end(); bond_it++ ){	
		Bond& b1 = *bond_it;
		//if(b1.bond_type==1 && !b1.is_boundary()){
		if(b1.bond_type==3){
			std::cout<<"bond vector index: "<<bond_it->index_in_vector<<std::endl;
			std::cout<<"bond index: "<<bond_it->bond_index<<std::endl;
			std::cout<<"bond nodes: "<<bond_it->nodes[0]<<" "<<bond_it->nodes[1]<<std::endl;		
			std::cout<<"rectangles: "<<b1.rectangles[0]<<" "<<b1.rectangles[1]<<std::endl;
			std::cout<<"rectangle nodes: "<<rectangles[b1.rectangles[0]].get(1)<<" "<<rectangles[b1.rectangles[0]].get(2)<<std::endl;
			std::cout<<"rectangle nodes: "<<rectangles[b1.rectangles[1]].get(0)<<" "<<rectangles[b1.rectangles[1]].get(3)<<std::endl;
			std::cout<<"boundary: "<<b1.is_boundary()<<std::endl;
			//std::cout<<"coop: "<<b1.coop_index[0]<<" "<<b1.coop_index[1]<<std::endl;
			
			std::cout<<std::endl;
		}
		
		
	}		
}
