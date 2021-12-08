#ifndef COST_HPP_
#define COST_HPP_

#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <mkl.h>
#include "mkl_lapacke.h"
#include "model.hpp"
#include "expm.hpp"
#include "MarkovChannel.hpp"
#include <algorithm> 
//#include "helper.hpp"



double cost_main(Model& m, std::vector<ProtocolParameter> protos);

double cost_validation(Model& m, std::vector<ProtocolParameter> protos, std::vector<ProtocolParameter> valids);

std::vector<double> cost_validation_comp(Model& m, std::vector<ProtocolParameter> protos,std::vector<ProtocolParameter> valids);

std::vector<double> cost_comp(Model& m, std::vector<ProtocolParameter> protos);

double cost(Model& m, ProtocolParameter proto, int val_run, Expm& Qexp);

double cost(Model& m, ProtocolParameter proto, int val_run, std::vector<double>& data, std::vector<double>& model);

double calculate_tau(std::vector<double> x, double v1, double v2);

double calculate_tau(std::vector<double> x, double v1, double v2, double v3, double v4);

double calc_rcond(Model& m);

double calc_eigs(Model& m);

double model_penalty(Model& m);

std::vector<double> transpose(Model& m);

int signC(std::vector<double> x);

double peak(std::vector<double> x, bool mini);

int inbounds(double output, double data, double SE);

double print_eigenvalues( char* desc, MKL_INT n, double* wr, double* wi );
void print_eigenvectors( char* desc, MKL_INT n, double* wi, double* v,
      MKL_INT ldv );


#endif
