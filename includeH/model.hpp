#ifndef MODEL_HPP_
#define MODEL_HPP_
#include "graph.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <stdio.h>
#include "MarkovChannel.hpp"
#include "math.hpp"

class Model {
	private:
	int P =2; //number of voltage dependent arguments, need to change to explore other vfuncs
	int id;
	int nargs;
	
	public:
    Model(int id, int N, std::vector<int> edgelist, int root, const SimulationParameters& sim_parms);
	Model(int id, int N, std::vector<int> edgelist, int root,std::vector<double> rs, std::vector<double> rk, std::vector<double> args, const SimulationParameters& sim_parms);
    Model();
	Model(const Model& m);
	

   
    Graph G;
	
	std::vector<double> rs;
	std::vector<double> rk;
	std::vector<double> vars;
    std::vector<int> C;
	std::vector<int> F;
    std::vector<double> r_vec;
	std::vector<double> s;
    std::vector<double> args;
	std::vector<double> Q;
    
	int n_states() {return G.N;};
    int n_edges() {return G.E;};
	int get_P() {return P;}
	int get_id() {return id;}
    double cost;
	void sobol(Math& math_params,SimulationParameters& sim_parms);
	std::vector<int> in_bounds_rates();
	double mot(double vm);
	void transition_matrix(double vm);
	void rate_vector();
	void initial_state(double vm);
	void vfunc(double vm);
	void print_Q(); 
	void printQexp(double*);
	Model perturb(Math& math_params, const SimulationParameters& sim_parms);
	
	friend std::ostream& operator<<(std::ostream& os, Model& m);
	
	
  };





#endif
