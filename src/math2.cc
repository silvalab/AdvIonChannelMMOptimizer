#include "math2.hpp"

namespace Math2 {

  std::mt19937 mt(0);
  std::mt19937 mt_restart;
  std::normal_distribution<double> params (0,1);
  std::normal_distribution<double> params_restart (0,1);
  std::normal_distribution<double> params_args(0,3);
  std::normal_distribution<double> params_args_restart(0,3);
  std::uniform_real_distribution<double> uni_real (0,1);
  std::uniform_real_distribution<double> uni_real_restart (0,1);

 
 
 
void save_stream(std::string fname){
	
	
	std::ofstream fout(fname + "seed.dat");
	fout << mt;
    fout.close();
	
	std::ofstream fout2(fname+"params.dat");
    fout2 << params;
    fout2.close();
	
	std::ofstream fout3(fname+"params_args.dat");
    fout3 << params_args;
    fout3.close();
	
	std::ofstream fout4(fname+"uni_real.dat");
    fout4 << uni_real;
    fout4.close();
}

void import_stream(char const* c){
	
	std::string path(c);
	std::ifstream fin(path+"seed.dat");
    fin >> mt;
    fin.close();
	
    std::ifstream fin2(path+"params.dat");
    fin2 >> params;
    fin2.close();
	
	std::ifstream fin3(path+"params_args.dat");
    fin3 >> params_args;
    fin3.close();
	
	std::ifstream fin4(path+"uni_real.dat");
    fin4 >> uni_real;
    fin4.close();
}
  
  
void norm_dist_params(int N, double* r){
	
	for (int i = 0; i < N; i++){
		r[i] = params(mt);
	}
	
} 

void norm_dist_args(int N, double* r){
	
	for (int i = 0; i < N; i++){
		r[i] = params_args(mt);
	}
}
 
void uni_dist(int N, double* r){
	
	for (int i = 0; i < N; i++){
		r[i] = uni_real(mt);
	}
} 

double uni_dist(){
	double  r[1];
	for (int i = 0; i < 1; i++){
		r[i] = uni_real(mt);
	}
	
	return *r;
} 

gsl_qrng* init_Sobol(int dim){
	 gsl_qrng * q = gsl_qrng_alloc (gsl_qrng_sobol, dim);
	  return q;
}
  
gsl_qrng* init_Sobol(int dim,int Sobol_indx){
   gsl_qrng * q = gsl_qrng_alloc (gsl_qrng_sobol, dim);
   double x[dim];
   for (int i = 0; i < Sobol_indx-1; i++){
	   gsl_qrng_get (q, x);	   
   }
   return q;
  
}

void Sobol(gsl_qrng* q, double* x, int& Sobol_indx){
	gsl_qrng_get (q, x);
	Sobol_indx++;
}
	  
void free_Sobol(gsl_qrng* q){
	gsl_qrng_free (q);
}
  
}