#include <iostream>
#include <array>
#include <vector>
#include "mkl.h"
#include "math.h"
#include "string.h"
#include "expm.hpp"

using namespace std;


double norm1(int n, double* a){
	
double colmax = 0;
double* result = new double[n*n];
cblas_dcopy(n*n,a,1,result,1);

for (int i = 0; i < n*n; i++){

	if( result[i] < 0){
		
		result[i] = result[i]*-1;
	}


}

double tempcolmax = 0;

for (int i = 0; i < n; i++){ // col index
	
	for(int j = 0; j < n*n; j+=n){
		tempcolmax += result[i+j];
	}
	
	
	
	if(tempcolmax >= colmax){
		
		colmax = tempcolmax;
		//std::cout << "colmax" << colmax << std::endl;
	}
	
	tempcolmax = 0;
} 	
	
	if(isinf(colmax)){
		
		colmax = 1.7e308;
		
	}
	
	
	delete [] result;
	return colmax;
	
	
}



void print(int n, double *a){
	
	for (int i = 0; i < n*n; i++){
		
		std::cout << a[i] << "\t";
		
		
	}
	std::cout<< std::endl;
}
void identity_matrix(int n,double* a){
	
for (int i = 0; i < n; i++){
	for (int j = 0; j < n; j++){
		
		if ( i ==j){	
		a[(i*n)+j] = 1;	
		}
		//std::cout << a[(i*N)+j] << "\t"; 
	}
	//std::cout << std::endl;
}		
}

vector<double>  get_pade_coefficients(int m){
//%get_pade_coefficients Coefficients of numerator P of Pade approximant
//%    C = get_pade_coefficients returns coefficients of numerator
//%    of [m/m] Pade approximant, where m = 3,5,7,9,13.
	//double c[]= {17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1};	

switch (m){
    case 3:
	{		
	double c []= {120.0, 60.0, 12.0, 1.0};
	std::vector<double> cvec (c, c + sizeof(c) / sizeof(double) );
	return cvec;
	}
    case 5:{
	double c [] = {30240, 15120, 3360, 420, 30, 1};
	std::vector<double> cvec (c, c + sizeof(c) / sizeof(double) );
	return cvec;
	}
     case 7:{
	double c [] = {17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1};
	std::vector<double> cvec (c, c + sizeof(c) / sizeof(double) );
	return cvec;
	}
    case 9:{
     double c[]= {(double)17643225600, (double)8821612800, 2075673600, 302702400, 30270240, 2162160, 110880, 3960, 90, 1};
	 std::vector<double> cvec (c, c + sizeof(c) / sizeof(double) );
	return cvec;
	 }
    case 13:{
     double c[] = {(double)64764752532480000, (double)32382376266240000, (double)7771770303897600, (double)1187353796428800, (double)129060195264000,  (double) 10559470521600, 
            (double)670442572800,   (double)33522128640,    (double)1323241920,
	(double)40840800,(double)960960,(double)16380,(double)182,(double)1};
	 std::vector<double> cvec (c, c + sizeof(c) / sizeof(double) );
	return cvec;
	 }
 }


}

void expm_params(double *T, int& s, int& m, int n, vector<double*>& Tpowers){ //need to eventually pass in Tpowers

	//Tpowers has {T^2, T^4, T^6, T^8, T^10}
	
	double coeff[] = { 9.9206e-06, 9.9413e-11, 2.2282e-16, 1.6908e-22, 8.8300e-36};
	
	s=0;
	
	double theta [] = {1.495585217958292e-002,  2.539398330063230e-001, 9.504178996162932e-001,  
    
    2.097847961257068e+000,5.371920351148152e+000};
	
	
	

	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, T, n, T, n, 0.0, Tpowers[0], n);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, Tpowers[0], n, Tpowers[0], n, 0.0, Tpowers[1], n);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, Tpowers[1], n, Tpowers[0], n, 0.0, Tpowers[2], n);
	
	//cout << "Tpowers[0]" << endl;
	//print(n,Tpowers[0]);
	//cout << "Tpowers[1]" << endl;
	//print(n,Tpowers[1]);
	//std::cout << "Tpowers[2]" << std::endl;
	//print(n,Tpowers[2]);
	//double d4 = LAPACKE_slange(CblasRowMajor,'1', n, n, Tpowers[1], n);
	double d4 = norm1(n,Tpowers[1]);
	//std::cout << "d4 by lapake" << "\t" << d4 << endl;
	//const int size = n*n;
	//int index = cblas_isamax(size, Tpowers[1], 1);
	//double d42 = Tpowers[1][index];
	//std::cout << "d4 isamax" << "\t" << d42 << endl;
	//double d6 = LAPACKE_slange(CblasRowMajor,'1', n, n, Tpowers[2], n);
	double d6 = norm1(n,Tpowers[2]);
	//std::cout << "Tpowers[2]" << std::endl;
	//print(n,Tpowers[2]);
	//std::cout << "d6 before pow" << "\t" << d6 <<std::endl;
	d4 = pow(d4,0.25);
	//d42 = pow(d42,0.25);
//std::cout << "d4" << "\t" << d4 << std::endl;
	//std::cout << "d4_2" << d42 << std::endl;
	d6 = pow(d6, 0.1667);
	//std::cout << "d6 after 1/6 pow" << "\t" << d6 << std::endl;
	double eta1 = max(d4,d6);
	//std::cout << "eta1" << "\t" << eta1 << std::endl;
	
	
	if (eta1 <= theta[0]){
		
		m = 3;
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, Tpowers[1], n, Tpowers[1], n, 0.0, Tpowers[3], n); //T8
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, Tpowers[3], n, Tpowers[0], n, 0.0, Tpowers[4], n);//T10
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, Tpowers[4], n, Tpowers[0], n, 0.0, Tpowers[5], n); //T12
		return;
	}
	if (eta1 <= theta[1]){
		
		m=5;
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, Tpowers[1], n, Tpowers[1], n, 0.0, Tpowers[3], n); //T8
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, Tpowers[3], n, Tpowers[0], n, 0.0, Tpowers[4], n);//T10
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, Tpowers[4], n, Tpowers[0], n, 0.0, Tpowers[5], n); //T12
		return;
	}
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, Tpowers[1], n, Tpowers[1], n, 0.0, Tpowers[3], n); 
	//cout << "Tpowers[3] == T^4*T^4" << endl;
	//print(n,Tpowers[3]);
	//double d8 = LAPACKE_slange(CblasRowMajor,'1', n, n, Tpowers[3], n);
	double d8 = norm1(n,Tpowers[3]);
	//std::cout << "d8" << "\t" << d8 << std::endl;
	d8 = pow(d8,0.125);
	//std::cout << "d8" << d8 << std::endl;
	double eta3 = max(d6,d8);
	
	if (eta3 <= theta[2]){
		
		m=7;
		
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, Tpowers[3], n, Tpowers[0], n, 0.0, Tpowers[4], n);//T10
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, Tpowers[4], n, Tpowers[0], n, 0.0, Tpowers[5], n); //T12
		return;
	}
	
	if (eta3 <= theta[3]){
		
		m = 9;
		
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, Tpowers[3], n, Tpowers[0], n, 0.0, Tpowers[4], n);//T10
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, Tpowers[4], n, Tpowers[0], n, 0.0, Tpowers[5], n); //T12
		return;
		
		
	}
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, Tpowers[3], n, Tpowers[0], n, 0.0, Tpowers[4], n);//T10
	//double d10 = LAPACKE_slange(CblasRowMajor,'1', n, n, Tpowers[4], n);
	double d10 = norm1(n,Tpowers[4]);
	d10 = pow(d10,0.1);
	
	double eta4 = max(d8,d10);
	double eta5 = min(eta3,eta4);
	//std::cout << "eta5" << eta5 << std::endl;
	s = max(ceil(log2(eta5/theta[4])),0.0);
	// ell function should call here if needed
	m = 13;
	
		
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, Tpowers[4], n, Tpowers[0], n, 0.0, Tpowers[5], n); //T12
	
	
}

void pade_approx(double *T, vector<double*>& Tpowers, int n, int m, double* F){
	double *U = (double*) malloc(n*n*sizeof(double));
	double *V = (double*) malloc(n*n*sizeof(double));
	vector<double> c = get_pade_coefficients(m);
	/* std::cout << "c" << std::endl;
	for(int i = 0; i < c.size(); i++){
		
		
		
	std::cout << c[i] << "\t";
	
	 }
	
	std::cout << std::endl; */
	memset(U,0.0, n*n*sizeof(double));
	memset(V,0.0,n*n*sizeof(double));
	identity_matrix(n,U);
	identity_matrix(n,V);
	
	//std::cout << "m is" << m << std::endl;
	
	cblas_dscal(n*n, c[1], U, 1);
	cblas_dscal(n*n, c[0], V, 1);
	
	
	 switch(m){
		
		case 13: {
			
			double *c13T6 = (double*) malloc(n*n*sizeof(double));
			double *c11T4 = (double*) malloc(n*n*sizeof(double));
			double *c9T2 = (double*) malloc(n*n*sizeof(double));
			
			
			double *sum = (double*) malloc(n*n*sizeof(double));
			cblas_dcopy(n*n,Tpowers[2], 1, c13T6,1);
			cblas_dcopy(n*n,Tpowers[1], 1, c11T4,1);
			cblas_dcopy(n*n,Tpowers[0], 1, c9T2,1);
			
			cblas_dscal(n*n, c[13], c13T6, 1);
			cblas_dscal(n*n, c[11], c11T4, 1);
			cblas_dscal(n*n, c[9], c9T2, 1);
			
			/* std::cout << "C13T6" << std::endl;
			print(n,c13T6);
			std::cout << "C11T4" << std::endl;
			print(n,c11T4);
			std::cout << "C9T2" << std::endl;
			print(n,c9T2); */
			
			cblas_daxpy(n*n,1,c13T6,1,c11T4,1);
			cblas_daxpy(n*n,1,c9T2,1,c11T4,1);
			
			cblas_dcopy(n*n,c11T4, 1,sum,1);
			//std::cout << "sum" << std::endl;
			//print(n,sum);
			double *ans2 = (double*) malloc(n*n*sizeof(double));
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, Tpowers[2], n, sum, n, 0.0, ans2, n); 
			//T6*(c13*T6+C11*T4+C9*T2)
			
			//std::cout << "ans2" << std::endl;
			//print(n,ans2);
			
			double *c7T6 = (double*) malloc(n*n*sizeof(double));
			double *c5T4 = (double*) malloc(n*n*sizeof(double));
			double *c3T2 = (double*) malloc(n*n*sizeof(double));
			
			cblas_dcopy(n*n,Tpowers[2], 1, c7T6,1);
			cblas_dcopy(n*n,Tpowers[1], 1, c5T4,1);
			cblas_dcopy(n*n,Tpowers[0], 1, c3T2,1);
			
			cblas_dscal(n*n, c[7], c7T6, 1);
			cblas_dscal(n*n, c[5], c5T4, 1);
			cblas_dscal(n*n, c[3], c3T2, 1);
			
			// std::cout << "c7T6" << std::endl;
			// print(n,c7T6);
			// std::cout << "c5T4" << std::endl;
			// print(n,c5T4);
			// std::cout << "c3T2" << std::endl;
			// print(n,c3T2);
			
			cblas_daxpy(n*n,1,c7T6,1,ans2,1);
			cblas_daxpy(n*n,1,c5T4,1,ans2,1);
			cblas_daxpy(n*n,1,c3T2,1,ans2,1);
			cblas_daxpy(n*n,1,U,1,ans2,1);
			
			
			//find U
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, T, n, ans2, n, 0.0, U, n); 
			//std::cout << "U" << std::endl;
			//print(n,U);
			
			
			
			
			
			
			
			//find V
			double *c13T6_2 = (double*) malloc(n*n*sizeof(double));
			double *c11T4_2 = (double*) malloc(n*n*sizeof(double));
			double *c9T2_2 = (double*) malloc(n*n*sizeof(double));
			double *sum2 = (double*) malloc(n*n*sizeof(double));
			double *ans3 = (double*) malloc(n*n*sizeof(double));
			cblas_dcopy(n*n,Tpowers[2], 1, c13T6_2,1);
			cblas_dcopy(n*n,Tpowers[1], 1, c11T4_2,1);
			cblas_dcopy(n*n,Tpowers[0], 1, c9T2_2,1);
			
			cblas_dscal(n*n, c[12], c13T6_2, 1);
			cblas_dscal(n*n, c[10], c11T4_2, 1);
			cblas_dscal(n*n, c[8], c9T2_2, 1);
			
			cblas_daxpy(n*n,1,c13T6_2,1,c11T4_2,1);
			cblas_daxpy(n*n,1,c9T2_2,1,c11T4_2,1);
			
			cblas_dcopy(n*n,c11T4_2, 1,sum2,1);
			/* std::cout << "sum2" << std::endl;
			print(n,sum2); */
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, Tpowers[2], n, sum2, n, 0.0, ans3, n); 
			
			
			double *c7T6_2 = (double*) malloc(n*n*sizeof(double));
			double *c5T4_2 = (double*) malloc(n*n*sizeof(double));
			double *c3T2_2 = (double*) malloc(n*n*sizeof(double));
			
			
			
			cblas_dcopy(n*n,Tpowers[2], 1, c7T6_2,1);
			cblas_dcopy(n*n,Tpowers[1], 1, c5T4_2,1);
			cblas_dcopy(n*n,Tpowers[0], 1, c3T2_2,1);
			
			cblas_dscal(n*n, c[6], c7T6_2, 1);
			cblas_dscal(n*n, c[4], c5T4_2, 1);
			cblas_dscal(n*n, c[2], c3T2_2, 1);
			
			cblas_daxpy(n*n,1,c7T6_2,1,ans3,1);
			cblas_daxpy(n*n,1,c5T4_2,1,ans3,1);
			cblas_daxpy(n*n,1,c3T2_2,1,ans3,1);
			cblas_daxpy(n*n,1,V,1,ans3,1);
			
			cblas_dcopy(n*n,ans3, 1, V,1); //find V
			
			/* std::cout << "U" << std::endl;
			print(n,U);
			std::cout << "V" << std::endl;
			print(n,V); */
			
			free(sum2);
			free(sum);
			free(ans2);
			free(ans3);
			free(c13T6);
			free(c11T4);
			free(c9T2);
			free(c7T6);
			free(c5T4);
			free(c3T2);
			free(c13T6_2);
			free(c11T4_2);
			free(c9T2_2);
			free(c7T6_2); 
			free(c5T4_2); 
			free(c3T2_2);
			break;
			
		}
		
		case 9 :{
				int count = 6;	
				double *ans = (double*) malloc(n*n*sizeof(double));
			for (int j = m; j >= 3; j-=2){
				/* std::cout << "j" << j << std::endl;
				
				std::cout << "count" << count << std::endl;
				std::cout << "cj1" << c[j] << std::endl;
				std::cout << "cj" << c[j-1] << std::endl;
			std::cout << "Tpowers[m-count}" << std::endl;
			print(n,Tpowers[j-count]); */
				cblas_daxpy(n*n,c[j],Tpowers[j-count],1,U,1);
				cblas_daxpy(n*n,c[j-1],Tpowers[j-count],1,V,1);
				/* std::cout << "U" << U << std::endl;
				print(n,U);
				std::cout << "V" << V << std::endl;
				print(n,V); */
				count--;
			}
			
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, T, n, U, n, 0.0, ans, n);
			cblas_dcopy(n*n,ans,1,U,1);
			free(ans);
			break;
			}
			
			
			case 7 :{
					int count = 5;
					double *ans = (double*) malloc(n*n*sizeof(double));
			for (int j = m; j >= 3; j-=2){
				/* std::cout << "j" << j << std::endl;
				
				std::cout << "cj1" << c[j] << std::endl;
				std::cout << "cj" << c[j-1] << std::endl;
			std::cout << "Tpowers[m-count}" << std::endl;
			print(n,Tpowers[j-count]); */
				cblas_daxpy(n*n,c[j],Tpowers[j-count],1,U,1);
				cblas_daxpy(n*n,c[j-1],Tpowers[j-count],1,V,1);
				/* std::cout << "U" << U << std::endl;
				print(n,U);
				std::cout << "V" << V << std::endl;
				print(n,V); */
				count--;
			}
			
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, T, n, U, n, 0.0, ans, n);
			cblas_dcopy(n*n,ans,1,U,1);
			free(ans);
			break;
			}
			
			
			case 5 :{
					int count = 4;
					double *ans = (double*) malloc(n*n*sizeof(double));
			for (int j = m; j >= 3; j-=2){
				/* std::cout << "j" << j << std::endl;
				
				std::cout << "cj1" << c[j] << std::endl;
				std::cout << "cj" << c[j-1] << std::endl;
			std::cout << "Tpowers[m-count}" << std::endl;
			print(n,Tpowers[j-count]); */
				cblas_daxpy(n*n,c[j],Tpowers[j-count],1,U,1);
				cblas_daxpy(n*n,c[j-1],Tpowers[j-count],1,V,1);
				/* std::cout << "U" << U << std::endl;
				print(n,U);
				std::cout << "V" << V << std::endl;
				print(n,V); */
				count--;
			}
			
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, T, n, U, n, 0.0, ans, n);
			cblas_dcopy(n*n,ans,1,U,1);
			free(ans);
			
			break;
			}
			
			case 3 :{
					int count = 3;
					double *ans = (double*) malloc(n*n*sizeof(double));
			for (int j = m; j >= 3; j-=2){
				/* std::cout << "j" << j << std::endl;
				
				std::cout << "cj1" << c[j] << std::endl;
				std::cout << "cj" << c[j-1] << std::endl;
			std::cout << "Tpowers[m-count}" << std::endl;
			print(n,Tpowers[j-count]); */
				cblas_daxpy(n*n,c[j],Tpowers[j-count],1,U,1);
				cblas_daxpy(n*n,c[j-1],Tpowers[j-count],1,V,1);
				/* std::cout << "U" << U << std::endl;
				print(n,U);
				std::cout << "V" << V << std::endl;
				print(n,V); */
				count--;
			}
			
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, T, n, U, n, 0.0, ans, n);
			cblas_dcopy(n*n,ans,1,U,1);
			free(ans);
			break;
			}
		}
		//std::cout << "U outside switch" << std::endl;
		//print(n,U);
		//std::cout<< "V outside switch" << std::endl;
		//print(n,V);
		
		double *U2 = (double*) malloc(n*n*sizeof(double));
		double *ans = (double*) malloc(n*n*sizeof(double));
		double *ans2 = (double*) malloc(n*n*sizeof(double));
		memset(ans2,0.0,n*n*sizeof(double));
		cblas_dcopy(n*n,U,1,U2,1);
		cblas_dscal(n*n, 2, U2, 1); //2U
		cblas_daxpy(n*n,-1,U,1,V,1);
		
		/* std::cout << "V-U" << std::endl;
		print(n,V); */
		
		
		/* std::cout << "2U" << std::endl;
		print(n,U2); */
		
		int ipiv[n];
		int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, V, n, ipiv );
		int info2 = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, V, n, ipiv); 
		
		/* std::cout << "V-Uinverse" << std::endl;
		print(n,V); */
		
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, V, n, U2, n, 0.0, ans2, n);
		
		/* std::cout << "inv(V-U)*U2" << std::endl;
		print(n,ans2); */
		memset(ans,0.0,n*n*sizeof(double));
		identity_matrix(n,ans);
		
		cblas_daxpy(n*n,1,ans2,1,ans,1);
		
		// std::cout << "I should be F" << std::endl;
		// print(n,ans);
		// std::cout << "F" << F << std::endl;
		
		
		cblas_dcopy(n*n,ans, 1, F,1);
		
		//std::cout << "F" << std::endl;
		//print(n,F);
		
		free(U2);
		free(ans);
		free(ans2);
		free(U);
		free(V);
}

void scale(int n, double * a, int factor){
	
	for (int i = 0; i < n*n; i++){
			
		a[i] = (a[i])/(pow(2,factor));	
			
		}
}


void Expm::calculate(double * matrix,double* result){
	
	//std::cout << "m" << m << std::endl;
	//std::cout << "s" << s << std::endl;
	//std::cout << "n" << n << std::endl;
	expm_params(matrix,s,m,n,Tpowers);

		/* std::cout << "m" << m << std::endl;
	std::cout << "s" << s << std::endl;
	std::cout << "n" << n << std::endl; */

	if (s != 0){
		
		
		scale(n,matrix,s);
		for (int j = 0; j < Tpowers.size(); j++){
			int factor = s*((2*j)+2);
			scale(n,Tpowers[j],factor);
			
			
		}
		
	}
	
	//std::cout << "FbeforePadeapproxcall" << F << std::endl;
	pade_approx(matrix, Tpowers, n, m,F);
	//std::cout << "FafterPadeapproxcall"  << std::endl;
	//print(n,F);
	cblas_dcopy(n*n,F, 1, Fcopy,1);
	for (int k = 0; k < s; k++){
		
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, F, n, F, n, 0.0, ans, n);
	cblas_dcopy(n*n,ans, 1, F,1);	
		
		}
	// std::cout << "Ffinal" << std::endl;
	// print(n,F);
	cblas_dcopy(n*n,F, 1, result,1);
}

   Expm::Expm(int size)  {
	// std::cout << "size" << size << std::endl;
	 n=size; 
	//std::cout << "n" << n << std::endl;
	T2 = (double*) malloc(n*n*sizeof(double));
	//memset(T2,0.0,n*n*sizeof(double));
	T4 = (double*) malloc(n*n*sizeof(double));
	//memset(T4,0.0,n*n*sizeof(double));
	T6 = (double*) malloc(n*n*sizeof(double));
	//memset(T6,0.0,n*n*sizeof(double));
	T8 = (double*) malloc(n*n*sizeof(double));	
	//memset(T8,0.0,n*n*sizeof(double));
	T10 = (double*) malloc(n*n*sizeof(double));	
	//memset(T10,0.0,n*n*sizeof(double));
	T12 = (double*) malloc(n*n*sizeof(double));
	//memset(T12,0.0,n*n*sizeof(double));
	F = (double*) malloc(n*n*sizeof(double));
	
    Tpowers = {T2,T4,T6,T8,T10,T12};		
		/* for (int i = 0; i < Tpowers.size(); i++){
			
			std::cout << Tpowers[i] << std::endl;
			
		} */
	
	Fcopy = (double*) malloc(n*n*sizeof(double));   
	ans = (double*) malloc(n*n*sizeof(double));   
}


Expm::~Expm()  {
	free(Fcopy);
	free(ans);
	 free(T2);
	 free(T4);
	 free(T6);
	 free(T8);
	 free(T10);
	 free(T12);
	 free(F);
}
/* int main() {
	
	
	
	double test [] = { -15.1655 ,   0.1618 ,   0.5230  ,       0  ,       0  ,       0    ,     0 ,   0.5423,
    0.0612,   -0.1618 ,        0  ,       0   ,      0   ,      0 ,        0  ,       0,
   14.7491 ,        0,   -0.7010,    2.4224,         0,         0,         0,    2.3167,
         0,         0 ,   0.1243,   -4.5455,   17.4880,         0,    0.0369,         0,
         0,         0,         0,    1.2422,  -21.3085,    1.9626,         0,         0,
         0 ,        0 ,        0,         0,    3.8205,  -12.0081,    0.0582,         0,
         0 ,        0 ,        0,    0.8808 ,        0,   10.0455,   -0.3664,   14.2824,
    0.3552 ,        0 ,   0.0538,         0,         0 ,        0 ,   0.2712 , -17.1413
	};

	
	int n = 8;
	
	Expm* testexp = new Expm(8);
	double * result = (double*) malloc(n*n*sizeof(double));   
	
	testexp->calculate(test,result);
	
	free(testexp);
free(result);
	
	return 0;
} */