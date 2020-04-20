#include "modified_sobol.hpp"
using namespace std;
// ----- SOBOL POINTS GENERATOR BASED ON GRAYCODE ORDER -----------------
// INPUT: 
//   N         number of points  (cannot be greater than 2^32)
//   D         dimension  (make sure that the data file contains enough data!!)      
//   dir_file  the input file containing direction numbers
//
// OUTPUT:
//   A 2-dimensional array POINTS, where
//     
//     POINTS[i][j] = the jth component of the ith point,
//   
//   with i indexed from 0 to N-1 and j indexed from 0 to D-1
//
// ----------------------------------------------------------------------

std::vector<double> sobol_points(unsigned N, unsigned D, char *dir_file)
{
  ifstream infile(dir_file,ios::in);
  if (!infile) {
    cout << "Input file containing direction numbers cannot be found!\n";
    exit(1);
  }
  char buffer[1000];
  infile.getline(buffer,1000,'\n');
  
  // L = max number of bits needed 
  unsigned L = (unsigned)ceil(log((double)N)/log(2.0)); 

  // C[i] = index from the right of the first zero bit of i
  //unsigned *C = new unsigned [N];
  std::vector<unsigned int> C(N,0);
  C[0] = 1;
  for (unsigned i=1;i<=N-1;i++) {
    C[i] = 1;
    unsigned value = i;
    while (value & 1) {
      value >>= 1;
      C[i]++;
    }
  }
  
  // POINTS[i][j] = the jth component of the ith point
  //                with i indexed from 0 to N-1 and j indexed from 0 to D-1
  std::vector<double> POINTS(D,0); 
  //for (unsigned j=0;j<=D-1;j++) POINTS[j] = 0; 

  // ----- Compute the first dimension -----
  
  // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
  //unsigned *V = new unsigned [L+1]; 
  std::vector<unsigned int> V(L+1,0);
  for (unsigned i=1;i<=L;i++) V[i] = 1 << (32-i); // all m's = 1

  // Evalulate X[0] to X[N-1], scaled by pow(2,32)
  //unsigned *X = new unsigned [N];
  std::vector<unsigned int> X(N,0);
  X[0] = 0;
  for (unsigned i=1;i<=N-1;i++) {
    X[i] = X[i-1] ^ V[C[i-1]];
    if(i == N-1)POINTS[0] = (double)X[i]/pow(2.0,32); // *** the actual points
  }
  
  X.clear();
  V.clear();
  
  
  // ----- Compute the remaining dimensions -----
  for (unsigned j=1;j<=D-1;j++) {
    
    // Read in parameters from file 
    unsigned d, s;
    unsigned a;
    infile >> d >> s >> a;
    //unsigned *m = new unsigned [s+1];
	std::vector<unsigned int> m(s+1,0);
    for (unsigned i=1;i<=s;i++) infile >> m[i];
    
    // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
    //unsigned *V = new unsigned [L+1];
	V.resize(L+1,0);
    if (L <= s) {
      for (unsigned i=1;i<=L;i++) V[i] = m[i] << (32-i); 
    }
    else {
      for (unsigned i=1;i<=s;i++) V[i] = m[i] << (32-i); 
      for (unsigned i=s+1;i<=L;i++) {
	V[i] = V[i-s] ^ (V[i-s] >> s); 
	for (unsigned k=1;k<=s-1;k++) 
	  V[i] ^= (((a >> (s-1-k)) & 1) * V[i-k]); 
      }
    }
    
    // Evalulate X[0] to X[N-1], scaled by pow(2,32)
    //unsigned *X = new unsigned [N];
	X.resize(N,0);
    X[0] = 0;
    for (unsigned i=1;i<=N-1;i++) {
      X[i] = X[i-1] ^ V[C[i-1]];
	  if(i == N-1) POINTS[j] = (double)X[i]/pow(2.0,32); // *** the actual points
      //        ^ j for dimension (j+1)
   }
    
    // Clean up
   
    //delete [] V;
    //delete [] X;
  }

  
  return POINTS;
}


/* int main(int argc, char **argv)
{
  if (argc != 4) {
    cout << endl << "input format: sobol N D FILENAME" << endl << endl;
    cout << "The program prints the first N sobol points in D dimensions." << endl;
    cout << "The points are generated in graycode order." << endl;
    cout << "The primitive polynomials and initial direction numbers are" << endl
	 << "given by the input file FILENAME." << endl << endl;
    return 0;
  }

  int N = atoi(argv[1]);
  int D = atoi(argv[2]);
  

  // display points
  cout << setprecision(20);
  //cout << setiosflags(ios::scientific) << setprecision(10);
  for (unsigned i=0;i<=D-1;i++) {
     cout << P[i] << " " ;
  }
  cout << endl;
} */
