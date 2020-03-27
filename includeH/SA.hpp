#ifndef SA_HPP_
#define SA_HPP_

#include "model.hpp"
#include "math.hpp"
#include "MarkovChannel.hpp"
#include "graph.hpp"
#include "cost.hpp"
#include "setup.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <stdio.h>
#include <fstream>
#include <memory>
#include "restart.hpp"
#include <unistd.h>
#include <stdlib.h>
#include <sys/wait.h>



class SA{
	
	private:
	std::string path; //where optimization results will be stored
	std::string version;
	std::string subpath;
	int time;
	double generalization_loss();
	double progress(int);
	int detect_overfitting();
	int inc3(std::vector<double> v);
	void print_snapshot(std::ofstream& out,Setup s, int i);
	void print_display(std::ofstream& out, Model m, Setup s, Math math_params, int i,SimulationParameters& sim_params);
	double fmin_val;
	std::vector<double> cost_val_history; // vectors for tracking when to stop the optimization to prevent over fitting
	std::vector<double> cost_tr_history;
	std::vector<double> PQ_history;
	std::vector<Model> models;
	std::vector<double> f_val;
	std::vector<int> r;
	std::vector<double> T_j;
	Model fmin_model;
	int S3_modelfiles(int N, std::string S3_bucket_path)
	public:
	SA(std::string,int time);
	SA(std::string,int time,Restart_params& restart_params);
	void anneal(Math& math_params, Model& m, Setup& s, SimulationParameters& sim_params, std::ofstream& out);
	void anneal_restart(Math& math_params, Setup& s, SimulationParameters& sim_params, Restart_params& restart_params, std::ofstream& out);
};



#endif