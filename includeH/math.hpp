#ifndef MATH_HPP_
#define MATH_HPP_

#include <fstream>
#include <iostream>
#include <random>
#include <cassert>
#include <vector>
//#include "gsl/gsl_qrng.h"
#include <string>
#include "mkl.h"
#include <cmath>
#include "modified_sobol.hpp"
#include <math.h>
# define M_PI           3.14159265358979323846  /* pi */

class Math {
private:
	long long random_counter_uni_real;
	int sobol_indx;
	int dim;
	char* dir_num_file;
	std::mt19937 mt_uni_real;
public:
	Math() {};
	Math(int dim,char* dir_num_file);
	Math(long long counter, int sobol_indx, int dim, char* dir_num_file);
	void norm_dist(std::vector<double>& r,double mu, double sigma); 
	void uni_dist(std::vector<double>& r);
	double uni_dist();
	double norm_dist();
	std::vector<double> sobol();
	int get_sobol_indx() {return sobol_indx;}
	long long get_random_counter_uni_real() {return random_counter_uni_real;}
	int get_dim() {return dim;}
};

#endif