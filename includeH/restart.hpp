#ifndef restart_HPP_
#define restart_HPP_

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include "helper.hpp"
class Restart_params{
private:
int iterations;
int times;
int chains;
int k_max;
int root;
int Sobol_indx;
int random_counter_params;
int random_counter_params_args;
int random_counter_uni_real;
std::vector<int> acceptance_times;
std::vector<double> acceptance_rates;
int times_since_exchange;
public:
int get_Sobol_indx() {return Sobol_indx;}
int get_random_counter_uni_real() {return random_counter_uni_real;}
int get_times() {return times;}
int get_iterations() {return iterations;}
int get_k_max() {return k_max;}
Restart_params(std::ifstream& restarttxtfile);
std::vector<std::vector<double>> model_rs;
std::vector<std::vector<double>> model_rk;
std::vector<std::vector<double>> model_args;	
std::vector<int> r;
std::vector<double> cost_val_history;
std::vector<double> cost_tr_history;
std::vector<double> PQ_history;
std::vector<double> fmin_rs;
std::vector<double> fmin_rk;
std::vector<double> fmin_args;	
};

#endif