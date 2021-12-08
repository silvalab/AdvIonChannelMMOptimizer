#ifndef MATH2_HPP_
#define MATH2_HPP_

#include <fstream>
#include <iostream>
#include <random>
#include <cassert>
#include "gsl/gsl_qrng.h"
#include <string>
namespace Math2 {



void save_stream(std::string);
void import_stream(char const*);


void norm_dist_params(int N, double* r); 
void norm_dist_args(int N, double* r);
void uni_dist(int N, double* r);
double uni_dist();


gsl_qrng* init_Sobol(int dim);

gsl_qrng* init_Sobol(int dim,int Sobol_indx);

void Sobol(gsl_qrng * q, double* x, int& Sobol_indx);

void free_Sobol(gsl_qrng*q);



/* void norm_dist_restart(int N, double* r); 
void norm_dist_args_restart(int N, double* r);
void uni_dist_restart(int N, double* r);
double uni_dist_restart(); */


}

#endif