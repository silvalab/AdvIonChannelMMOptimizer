#include <iostream>
#include "model.hpp"
#include <random>
#include <iomanip>
#include "helper.hpp"


int in_bounds(double val, double lower, double higher){
	
	if (val > lower && val < higher)
	return 1;
	
	return 0;
}

  

  
Model::Model(int id, int N, std::vector<int> edgelist,int root, const SimulationParameters& sim_params){
	this->id = id;
	
    G.set_topology(N,root,edgelist);
	
	rs.resize(P*(G.N)); 
	rs[0] = 0;
	rs[1] = 0;
    rk.resize(P*G.E);
    C.resize(G.N);
	F.resize(G.N);
	s.resize(G.N);
    C[root] = 1.0; F[root] = 1.0; 
	r_vec.resize(2*G.E*P);
	//disp(r_vec);
	
  
    nargs = 2*(P-1);
	
	args.resize(nargs);
	vars.resize(nargs);
	Q.resize(N*N);
    cost = 0;  	
	
}
  
  
Model::Model(int id, int N, std::vector<int> edgelist, int root,std::vector<double> rs, std::vector<double> rk, std::vector<double> args, const SimulationParameters& sim_parms){
	
	
	this->id = id;
	G.set_topology(N,root,edgelist);
	this->rs =rs;
    this->rk = rk;
	this->args = args;

	
	C.resize(G.N);
	F.resize(G.N);
    C[G.root] = 1.0; F[G.root] = 1.0; 
	s.resize(G.N);
	r_vec.resize(2*G.E*P);
	//disp(r_vec);
	nargs = 2*(P-1);
	args.resize(nargs);
	vars.resize(nargs);
	Q.resize(N*N);
    cost = 0;  	
	  

}	
  

Model::Model(){
	nargs = 0;
    cost = 0;
}

Model::Model(const Model& m){
	
	this->id = m.id;
	G = m.G; 
	
	this->rs = m.rs;
    this->rk = m.rk;
	this->C = m.C;
	this->F = m.F;
	this->s = m.s;
    C[G.root] = 1.0; F[G.root] = 1.0; 
	this->r_vec = m.r_vec;
	
	args = m.args;
	vars = m.vars;
    nargs = m.nargs;
	this->Q = m.Q;
    cost = m.cost;  	
	
	
	
}

  

std::ostream& operator<<(std::ostream& os, Model& m){
   int P=2;
   int N=m.G.N;
   //os << "Model ID:" << "\t" << m.id << std::endl;
   os << "Num_Nodes:" << "\t" << N << std::endl;
   int E = m.G.E;
   os << "Num_Edges:" << "\t" << E << std::endl;
   int root = m.G.root;
   os << "Root:" << "\t" << root << std::endl;
	
	std::vector<int> ic (2*(E*N),0);
	
	
	for(int i=0; i<m.G.E; i++) {
		int e1 = m.G.edges[i].V1; 
		int e2 = m.G.edges[i].V2;
		ic[2*m.G.E*(e1) + 2*i] = -1;
		ic[2*m.G.E*(e2) + 2*i] = 1;
		ic[2*m.G.E*(e1) + 2*i+1] = 1;
		ic[2*m.G.E*(e2) + 2*i+1] = -1;
	}
	os << "ic = [\n";
	for (int i=0; i<m.G.N; i++) {
		os << "[";
		for (int j=0; j<2*m.G.E; j++) {
			os << ic[i*2*m.G.E+j] << "\t";
		}
		os << "]\n";
	}
	os << (" ] \n"); 
	
	os << "rs =" << "[" << std::endl;
    for (int i=0; i<N; i++) {
		os << " ["; 
      for (int j=0; j<P; j++) {
        os << std::setprecision(18) << m.rs[i*P+j] << "\t";
      }
      os << " ] \n"; 
    }
   os << " ]" << std::endl;
    
	os << "rk = [" << std::endl;
    for (int i=0; i<E; i++) {
		os << " ["; 
      for (int j=0; j<P; j++) {
       os << std::setprecision(18) << m.rk[i*P+j] << "\t";
      }
      os << " ] \n";
    }
	os << " ]" << std::endl;
	
    os << "G = [ " << std::endl;
    for (int i=0; i<N; i++) {
      os << m.C[i] << "\t";
    }

  os << " ] \n" <<std::endl;

	os << "args = [" <<std::endl;
	for (int i = 0; i < 2; i++){
	   os << std::setprecision(18) << m.args[i] << "\t";


			}
	os << " ] " << std::endl;
	
	
   
   
	return os;
}

void Model::vfunc(double vm){
   
    vars[0] = 1.0;
    
    double b = args[1];
    double a = args[0];
    vars[1] = tanh((vm+a)/b);
   
 
}

void Model::initial_state(double vm){

   
	
    vfunc(vm);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, s.size(), P,1.0, rs.data(), P, vars.data(), 1, 0.0, s.data(), 1);
	
	for (int i=0; i<s.size(); i++) s[i] = exp(s[i]); 
		 
		
		double scale = 1.0/cblas_dasum(s.size(), s.data(), 1);
		cblas_dscal(s.size(), scale, s.data(), 1); 
		
  }



void Model::rate_vector(){
    
	int N=G.N;
	int E=G.E;
	double* p_rs = rs.data(); //returns pointer to the contiguous memory that the vector owns
	double* p_rk = rk.data(); //returns pointer to the contiguous memory that the vector owns
    std::vector<double> b(2*E*P,0);   //create storage for computation
	double* p_b = b.data();
	std::fill(r_vec.begin(), r_vec.end(), 0);
	double* p_rvec = r_vec.data();
	//**************see equation 16 in Teed et al. for the following computations******************
    for ( int i=0; i<E; i++ ) {
        int e1=G.edges[i].V1, e2=G.edges[i].V2;
			//vdSub(num_params, vector a, vector b, storage) =a-b
			vdSub(P, &p_rs[e2*P], &p_rs[e1*P], &p_b[i*P]); //equivalent  to the (I_e)^T*rs
    }
	std::copy(rk.begin(),rk.end(),&p_b[P*E]); //append on rk    =[(I_e)^T*rs; rk]
	
    for ( int i=0; i<E; i++ ) { //after creating [D;abs(D)] matrix has in Teed et al. and inverting, the following vector
								 //additions/subtractions are the equivalent to the action of [D ;abs(D)]^-1*[(I_e)^T*rs; rk]. 
								 //For D formation in Teed, use vdAdd then vdSub, 
								 //assign parameters with D matrix formation as in Menon vdSub then vdAdd.
        vdAdd(P, &p_b[P*(i+E)], &p_b[P*i], &p_rvec[P*(2*i+0)]);
        vdSub(P, &p_b[P*(i+E)], &p_b[P*i], &p_rvec[P*(2*i+1)]);
      }


      cblas_dscal(2*E*P, .5, p_rvec, 1);  //when adding the sum and difference equations each rate is the form of 2(r_ij)
	 
     
    
 }

void Model::transition_matrix(double vm){
	//completes equation 17 in Teed et al.
    int N=G.N, E=G.E;
	std::fill(Q.begin(), Q.end(), 0);
    vfunc(vm); //get voltage dependence
	rate_vector(); //calculate rate rector

	std::vector<double> e_vec (2*E,0);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 2*E, P, 1.0,
        r_vec.data(), P, vars.data(), 1, 0.0, e_vec.data(), 1); //apply voltage dependence to rates

    for (int i=0; i<2*E; i++) e_vec[i] = exp(e_vec[i]); //exponenate the rates
    

    for (int i=0; i<E; i++) { //populate the Q matrix with the rates. Q(i,i) is the negative sum of Q(:,i)
      int e1 = G.edges[i].V1; 
	  int e2 = G.edges[i].V2;
	  
	  Q[N*e1+e1] -= e_vec[2*i];  //subtract rate from  Q(e1,e1) 
	  Q[N*e2+e2] -= e_vec[2*i+1]; //subtract rate from  Q(e2,e2) 
	  Q[N*e2+e1] += e_vec[2*i]; 	//populate Q(e2,e1)
	  Q[N*e1+e2] += e_vec[2*i+1]; //populate Q(e1,e2)
    }
}
	  	  
double Model::mot(double vm){
	
	int E = G.E;
	std::vector<double> e_vec(2*E,0);
	vfunc(vm);	
	rate_vector();
	
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 2*E, P, 1.0,
        r_vec.data(), P, vars.data(), 1, 0.0, e_vec.data(), 1); //apply voltage dependence to rates

    for (int i=0; i<2*E; i++) e_vec[i] = exp(e_vec[i]);
	
	
	std::vector<double> leaving_rates;
	
	
	//e_vec[2*i]; 	//populate Q(e2,e1)
   // e_vec[2*i+1]; //populate Q(e1,e2)
	for(int i = 0; i < E; i++){
		int e1 = G.edges[i].V1; 
		int e2 = G.edges[i].V2;
		if(e1 == G.root){
			leaving_rates.push_back(e_vec[2*i+1]);
		}
		else if(e2 == G.root){
			
			leaving_rates.push_back(e_vec[2*i]);
			
		}
	}
	
	
	
	if(leaving_rates.size() == 0){
		
		std::cout << "Warning: no leaving edges from open state found!" << std::endl;
		
		return 0;
		
	}
	
	
	double scale = 0;
	for(int i = 0; i < leaving_rates.size(); i++){
		
		scale += leaving_rates[i];
		
		
	}

	
	return 1/scale;
	
	
}
	  
void Model::sobol(Math& math_params,SimulationParameters& sim_params){
	
	double rate_min = sim_params.rate_min;
	double rate_max = sim_params.rate_max;	
	double arg_min = sim_params.arg_min;
	double arg_max = sim_params.arg_max;
	int N = G.N;
	int E = G.E;
	

	std::vector<double> Sobol_nums = math_params.sobol();
	std::cout << "Sobol nums" << std::endl;
	for(int i = 0; i < (math_params.get_dim()); i++){
		
		std::cout << Sobol_nums[i] << "\t";
	}
	std::cout << std::endl;
	for(int i = 0; i < (((N-1)*P)); i++){
		rs[i+2] = (Sobol_nums[i]*(rate_max-rate_min))+ rate_min;
		
		
	}
	
	for(int i = ((N-1)*P); i < (((N-1)+E)*P); i++){
		rk[i-((N-1)*P)] = (Sobol_nums[i]*(rate_max-rate_min))+ rate_min;
		
		
	}
	
	for(int i = (((N-1)+E)*P); i < ((((N-1)+E)*P) + nargs); i++){
		args[i-(((N-1)+E)*P)] = (Sobol_nums[i]*(arg_max-arg_min))+ arg_min;
		
	}

	
	//memory dynamically allocated in modified_sobol.cc
}	
	  
Model Model::perturb(Math& math_params,const SimulationParameters& sim_params){
	
	Model n = Model(*this);
    int N = G.N;
	int E = G.E;
	

	
	
	std::vector<double> r (P*(N+E),0);
	std::vector<double> ra (nargs,0);
	
	math_params.norm_dist(r,sim_params.update_mu_rates,sim_params.update_std_rates);
    
	
	
    for(int i=0; i<(P*N); i++) {
		n.rs[i] += r[i];
		
    }
	 for(int i=(P*N); i<(P*(N+E)); i++) {
		n.rk[i-(P*N)] += r[i];
		
    }
	math_params.norm_dist(ra,sim_params.update_mu_args,sim_params.update_std_args);
  
	  
    for (int i = 0; i < nargs; i++) {
		if(in_bounds(n.args[i] + ra[i], sim_params.arg_min, sim_params.arg_max)){
			n.args[i] += ra[i];
		}
    } 
	
    return n;
	
}

void Model::print_Q(){
	for (int i = 0; i < (G.N); i++){	
		for (int j = 0; j < G.N; j++){
				std::cout << Q[(i*G.N)+j] << "\t";
		}
		std::cout << std::endl;
	}
}

void Model::printQexp(double* Qm){

	for (int i = 0; i < (G.N); i++){
		for (int j =0; j < (G.N); j++){
			std::cout << Qm[(i*G.N)+j] << "\t";
		}
	std::cout << std::endl;
	}
}
