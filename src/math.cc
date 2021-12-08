#include "math.hpp"
#include <iomanip>  
#include "MarkovChannel.hpp"

  
Math::Math(int dim,char* dir_num_file){
	

	mt_uni_real.seed(0);
	this->random_counter_uni_real = 0;
	this->sobol_indx = 1;
	this->dim = dim;
	this->dir_num_file = dir_num_file;
	
	
}
Math::Math(long long counter, int sobol_indx, int dim,char* dir_num_file){
	 
	
	mt_uni_real.seed(0);

	 
	this->random_counter_uni_real = 0;
	
	  for (int i = 0; i < counter; i++){
		 //std::cout << "Random value:" << mt_uni_real() << std::endl;
		 mt_uni_real();
		 this->random_counter_uni_real++;
	}
	

	this->sobol_indx = sobol_indx;
	this->dim = dim;
	this->dir_num_file = dir_num_file;
 }
 

void Math::norm_dist(std::vector<double>& r,double mu, double sigma){
	
	int times = r.size()/2;
	int remainder = r.size()%2;
	
	for (int i = 0; i < times; i++){
		double u1 = uni_dist();
		double u2 = uni_dist();
		r[i*2] = (sigma*sqrt (-2 * log(u1)) * cos (2 * M_PI * u2))+ mu;
		r[(i*2)+1] = (sigma*sqrt (-2 * log(u1)) * sin (2 * M_PI * u2)) + mu;
		//random_counter_uni_real+=2;		
	}
	for (int i = 0; i < remainder; i++){
		
		double u1 = uni_dist();
		double u2 = uni_dist();
		r[r.size()-1] = (sigma*sqrt (-2 * log(u1)) * cos (2 * M_PI * log(u2)))+ mu;
		//random_counter_uni_real+=2;
	}
	
}  


 
void Math::uni_dist(std::vector<double>& r){
	
	for (int i = 0; i < r.size(); i++){
		r[i] = ((double)mt_uni_real())/(mt_uni_real.max());
	}
	random_counter_uni_real+=r.size();
} 


double Math::norm_dist(){
	double u1 = uni_dist();
	double u2 = uni_dist();
	double  r[1];
	r[0] = sqrt (-2 * log(u1)) * cos (2 * M_PI * log(u2));
	return *r;
}


double Math::uni_dist(){
	double  r[1];
	r[0] = ((double)mt_uni_real())/(mt_uni_real.max());
	random_counter_uni_real++;
	return *r;
} 


std::vector<double> Math::sobol(){
	
	
	std::vector<double> points = sobol_points(sobol_indx,dim,dir_num_file);
	sobol_indx++;
	
	return points;
}
	  
  
