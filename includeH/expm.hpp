#include <iostream>
#include <mkl.h>
#include <math.h>
#include <vector>
class Expm  {
		public:
		double *T2,*T4, *T6, *T8, *T10, *T12;
		double *F;
		std::vector<double*> Tpowers;	
        int m;
		int s;
		int n;
		double *ans;
	    double *Fcopy;

 
	void calculate(double*,double*);
	/* expm_params(double *,int& s, int& m,int n, vector<double*>& Tpowers);
	pade_approx(double* , vector<double*>&, int n , int m, double *); */
    Expm(int);
    ~Expm();



};