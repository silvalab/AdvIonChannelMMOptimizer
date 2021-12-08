#include "cost.hpp"
#include <sstream>
#include "helper.hpp"


/////////////////////////////////////////////////////////////////////////////

//Helper Functions for step_exp

int signC(std::vector<double> x){

	
	int idxmax = cblas_idamax(x.size(), x.data(), 1);
	double resmax = x[idxmax];

	int idxmin = cblas_idamin(x.size(), x.data(), 1);
	double resmin = x[idxmin];


	if (resmax > 0 && resmin > 0) {
		return 0;
	}
	else if (resmax > 0 && resmin < 0) {
		return 1;
	}
	else {
		return 0;
	}
}


double peak(std::vector<double> x, bool mini=false){ 
  int idx=0;
  if (mini) {
    idx = cblas_idamin(x.size(), x.data(), 1);
  } 
  else {
	idx = cblas_idamax(x.size(), x.data(), 1);
  }
  
  double res = x[idx]; 
  return res;
}

double calculate_tau(std::vector<double> x, double v1, double v2){
	int idx = cblas_idamax(x.size(), x.data(), 1);
	cblas_dscal(x.size(), 1/x[idx], x.data(), 1);
	const double eps=1e-4; 
	double res;
	std::vector<double> diff_1 (x.size(),0);
	std::vector<double> diff_2 (x.size(),0);
 	if (v1 > v2){
	int count = 0;
        for(int i=idx; i < x.size(); i++){
			diff_1[i-idx] = x[i]-v1;
			diff_2[i-idx] = x[i]-v2;
			count++;
		}
        int tau1 = cblas_idamin(count, diff_1.data(), 1);
        int tau2 = cblas_idamin(count, diff_2.data(), 1);
        
		res = abs(tau2-tau1);
	}
	
	else { // 
	int count = 0;
		for(int i=0; i <= idx ; i++){
			diff_1[i] = x[i]-v1;
			diff_2[i] = x[i]-v2;
			count++;
		}
		int tau1 = cblas_idamin(count,diff_1.data(), 1);
	
		int tau2 = cblas_idamin(count,diff_2.data(), 1);
        
		res = abs(tau2-tau1);
	}

	
	return res;
}



int inbounds(double output, double data, double SE){
	double lower = data-SE;
	double higher = data + SE;
	if(output >= lower && output <= higher){
		return 1;
	}
	else{
		return 0;
	}
}




void step_exp(Model& m, std::vector<ProtocolStep> steps, std::vector<double>& out, Expm& Qexp){ 

	int N = m.n_states(); 
	double t;
	double current;
	std::vector<double> y0 (N,0);
	std::vector<double> y (N,0);
	y = m.s;
	

	for ( int i=0; i<steps.size(); i++ ) {
		y0 = y;
		double dt = steps[i].get_dt();
		double vm=steps[i].get_vm(); 
		m.transition_matrix(vm);
		if ( steps[i].get_stype() == NONE ) {
	
			
			std::vector<double> Qpointer (N*N,0);
			cblas_dscal(m.G.N*m.G.N, dt, m.Q.data(), 1);
			Qexp.calculate(m.Q.data(),Qpointer.data());

			cblas_dgemv(CblasRowMajor, CblasNoTrans, m.G.N, m.G.N, 1.0, Qpointer.data(), m.G.N, y0.data(), 1, 0.0, y.data(), 1);
			current = y[m.G.root]*(vm-54.4);
			
		}	
								
		else {
			int n_steps = ceil(dt / steps[i].get_stepsize());
			
			t = steps[i].get_stepsize();
			
			std::vector<double> times(n_steps,0);
			std::vector<double> vals(n_steps,0);
			vals[0] = y0[m.G.root];
			times[0] = 0;
			std::vector<double> Qpointer (N*N,0);
			cblas_dscal(m.G.N*m.G.N, t, m.Q.data(), 1);
			Qexp.calculate(m.Q.data(),Qpointer.data());
			
	 
			for ( int j=1; j<n_steps; j++ ) {
				cblas_dgemv(CblasRowMajor, CblasNoTrans, m.G.N, m.G.N, 1.0, Qpointer.data(), m.G.N, y0.data(), 1, 0.0, y.data(), 1);
				vals[j] = y[m.G.root];
				times[j] = j*t; 
				
				y0 = y;
			}
			
			if (steps[i].get_stype() == TAUD){
		
				for (int j =0; j<n_steps; j++){
					vals[j] = vals[j]*(vm-54.4); 
					
				}
				
				double taucons = 0.50; //current we want to look for
				for (int j = 0; j < n_steps; j++) {
					vals[j] = ((vals[j] - (taucons*current))*-1); 
					
				}
		
			
				if ( (vals[0] > vals[n_steps/2]) && (vals[0] != vals[n_steps-1])){ //difference vector decreasing

					int idx = cblas_idamin(n_steps, vals.data(), 1); //find min of the differences 
					int idx2;
	

					

					if ( signC(vals)) { //if the sign changes, then difference vector contains desired value

						for(int i=0; i < n_steps; i++){
							if ( vals[i] < 0){
								idx2 = i;
								break;				
							}
						}
						//std::cout << "option2" << std::endl;
						double timer = times[idx2];			
						out.push_back(timer);
					}
					else { //return the maximum solving time
						double time = times[n_steps-1];
					//	std::cout << "option3" << std::endl;
						out.push_back(time); 
						//std::cout << "maximum solving time stagnates" << std::endl;
					}
				}
				else { //if the difference vector keeps increasing or is the same, return the maximum time
					double time1 = times[n_steps-1];
				//	std::cout << "option4" << std::endl;
					out.push_back(time1); 
					//std::cout << "diff vector keeps increasing or same" << std::endl;
				}

			}
	
			else if (steps[i].get_stype() == LATEINA){
			
				for (int j = 0; j < n_steps; j++) {
					vals[j] = vals[j]*(vm-54.4); 
				}
				double scale = 0;
				for (int i = 0; i < n_steps; i++) {
					scale = std::min(scale, vals[i]);
				}
				for (int i = 0; i < n_steps; i++) {
					vals[i] /= scale;
				}
				out.push_back(vals[n_steps-1]);
			}
	
			else if ( steps[i].get_stype() == TAU ){
				double res;
				if(steps[i].get_extra_args().size() == 2){
					double arg1=steps[i].get_extra_args()[0], arg2=steps[i].get_extra_args()[1];
					res = calculate_tau(vals, arg1, arg2);
				}
				
					out.push_back(res*steps[i].get_stepsize());
			}

			else if ( steps[i].get_stype() == PEAK ) {
				double res = peak(vals);
				out.push_back(res);
			}

			else if ( steps[i].get_stype() == MIN ) {
				double res = peak(vals, true);
				out.push_back(res);
			}

			else if ( steps[i].get_stype() == TRACE ) {
				//std::cout << "trace " << std::endl;
				for ( int j=0; j<n_steps; j++ ) {
					out.push_back(vals[j]);
				}
			}
			else if ( steps[i].get_stype() == MAXPO ) {
				double res = peak(vals);
				out.push_back(res);
			}
			
			
		}

	}

}
/////////////////////////////////////////////////////////////////////////////////////////
//calculate cost per model and per voltage protocol

double cost(Model& m, ProtocolParameter proto, int valid_run,  Expm& Qexp){
	std::vector<double> output;
	int N = m.G.N;
	double err = 0;
	
	if(valid_run){
		for ( int i=0; i<proto.sweeps_full; i++ ) { //populate output vector
	  
			m.initial_state(proto.v0);
			step_exp(m, proto.full_steps[i],output, Qexp); 
			
		}

		int idx = 0;
		int n = output.size();
	
		double dx = 0;
		if ( proto.normalize) {
			double scale = 0;
			for (int i = 0; i < output.size(); i++) {
				scale = std::max(scale, output[i]);
			}
			for (int i = 0; i < output.size(); i++) {
				output[i] /= scale;
			}
		}
		int size = (proto.n_traces_full-proto.n_traces);
		//need to only push out the error of the last data points appended for validation
		if(proto.has_validation_points && proto.full_steps[0][0].get_stype() == TRACE){
		
			for (int i = 0; i < output.size(); i++) {
					
					if(contains(proto.trace_val_indxs,i)){
						if(!inbounds(output[i],proto.full_data[i],proto.full_SE[i])){
							dx = (output[i] - proto.full_data[i])/(proto.full_data[i]);
							err += (1.0 / size) * dx * dx;
						}//
					}
				}
		
		}
		else{
			for (int i = proto.n_traces; i < output.size(); i++) {
				if(!inbounds(output[i],proto.full_data[i],proto.full_SE[i])){
					dx = (output[i] - proto.full_data[i])/(proto.full_data[i]);
					err += (1.0 / size) * dx * dx;
				}//
			
			}
		
		}
		//disp(output);
	}
	else{
		for ( int i=0; i<proto.sweeps; i++ ) { //populate output vector
	  
			m.initial_state(proto.v0);
			step_exp(m, proto.steps[i], output, Qexp); 
			
		}

		int idx = 0;
		int n = output.size();
		double dx = 0;
		err = 0;
		if ( proto.normalize) {
			double scale = 0;
			for (int i = 0; i < output.size(); i++) {
				scale = std::max(scale, output[i]);
			}
			for (int i = 0; i < output.size(); i++) {
				output[i] /= scale;
			}
		}
		if(proto.has_validation_points && proto.full_steps[0][0].get_stype() == TRACE){		
				int size = (proto.n_traces-(proto.n_traces_full-proto.n_traces));
					for (int i = 0; i < output.size(); i++) {
						if(!contains(proto.trace_val_indxs,i)){
								if(!inbounds(output[i],proto.data[i],proto.SE[i])){
									dx = (output[i] - proto.data[i])/(proto.data[i]);
									err += (1.0 / size) * dx * dx;
								}//
						}
					}			
		}
		else{
			
			for (int i = 0; i < output.size(); i++) {
					if(!inbounds(output[i],proto.data[i],proto.SE[i])){
						dx = (output[i] - proto.data[i])/(proto.data[i]);
						err += (1.0 / output.size()) * dx * dx;
					}//
					
			}
		}
	}
	
	
	if (std::isnan(err))
			return 1e6;
	
  
  return err;  
  

}

double cost(Model& m, ProtocolParameter proto, int valid_run, std::vector<double>& data, std::vector<double>& model){
	

	std::vector<double> output;
	int N = m.G.N;
	Expm Qexp(m.G.N);
	
	
	if(valid_run){
		
		for ( int i=0; i<proto.sweeps_full; i++ ) {
			
			m.initial_state(proto.v0);
			step_exp(m,proto.full_steps[i], output, Qexp);
			
		}
		
		
		if (proto.normalize) {
			double scale = 0;
			for (int i = 0; i < output.size(); i++) {
				scale = std::max(scale, output[i]);
			}
			for (int i = 0; i < output.size(); i++) {
				output[i] /= scale;
			}
		}
		int size = (proto.n_traces_full-proto.n_traces);
		//std::cout << "size" << size << std::endl;
		data.resize(size,0);
		model.resize(size,0);
		if (proto.has_validation_points && proto.full_steps[0][0].get_stype() == TRACE){
			for(int i = 0; i < proto.trace_val_indxs.size(); i++){
				int indx = proto.trace_val_indxs[i];
				model[i] = output[indx];
			}
			
			for (int i = proto.n_traces; i < proto.n_traces_full; i++) {
				data[i-proto.n_traces] = proto.full_data[i];
				
			}
		}
		else{
			
			for (int i = proto.n_traces; i < output.size(); i++) {
				model[i-proto.n_traces] = output[i];
				data[i-proto.n_traces] = proto.full_data[i];
			}
			
		}
	
	double err = 0;
	double dx = 0;
	if(proto.has_validation_points && proto.full_steps[0][0].get_stype() == TRACE){
		
		for (int i = 0; i < output.size(); i++) {
				
				if(contains(proto.trace_val_indxs,i)){
					if(!inbounds(output[i],proto.full_data[i],proto.full_SE[i])){
						dx = (output[i] - proto.full_data[i])/(proto.full_data[i]);
						err += (1.0 / size) * dx * dx;
					}//
				}
			}
		
	}
	else{
		for (int i = proto.n_traces; i < output.size(); i++) {
			if(!inbounds(output[i],proto.full_data[i],proto.full_SE[i])){
				dx = (output[i] - proto.full_data[i])/(proto.full_data[i]);
				err += (1.0 / size) * dx * dx;
			}//
			
		}
		
	}
	//disp(output);
	if (std::isnan(err))
			return 1e6;
	return err;
		
	}
	
	else{
		for ( int i=0; i<proto.sweeps; i++ ) {
			
			m.initial_state(proto.v0);
			step_exp(m, proto.steps[i], output, Qexp);
			
		}
		
		data.resize(output.size());
		model.resize(output.size());
		
		
		if (proto.normalize) {
			double scale = 0;
			for (int i = 0; i < output.size(); i++) {
				scale = std::max(scale, output[i]);
			}
			for (int i = 0; i < output.size(); i++) {
				output[i] /= scale;
			}
		}


		for (int i = 0; i < output.size(); i++) {
			model[i] = output[i];
			data[i] = proto.data[i];
		}
		
		
	double err = 0;
	double dx = 0;
	if(proto.has_validation_points && proto.full_steps[0][0].get_stype() == TRACE){
		
		int size = (proto.n_traces-(proto.n_traces_full-proto.n_traces));
		for (int i = 0; i < output.size(); i++) {
			if(!contains(proto.trace_val_indxs,i)){
					if(!inbounds(output[i],proto.data[i],proto.SE[i])){
						dx = (output[i] - proto.data[i])/(proto.data[i]);
						err += (1.0 / size) * dx * dx;
					}//
			}
		}
		
	}
	else{
		
		for (int i = 0; i < output.size(); i++) {
				if(!inbounds(output[i],proto.data[i],proto.SE[i])){
					dx = (output[i] - proto.data[i])/(proto.data[i]);
					err += (1.0 / output.size()) * dx * dx;
				}//
				
		}
	}
	
//	disp(output);
	if (std::isnan(err))
			return 1e6;
	return err;
	}
	
}

std::vector<double> transpose(Model& m){
	
	int N = m.G.N;
	
	
	std::vector<double> col_major_Q(N*N,0);

	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			
			col_major_Q[j*N+i] = m.Q[i*N+j];
		}
		
	}
	
	return col_major_Q;
}

double calc_rcond(Model& m){
     int i, info, n, lda;
    double anorm, rcond;
	n = m.G.N;
	
	std::vector<double> colQ = transpose(m);
	   
	lda = n;	
    int iw[n];
    double w[4*n];
    /* Computes the norm of m.Q */
    anorm = dlange_("1", &n, &n,colQ.data(), &lda, w);
    /* Modifies m.Q in place with a LU decomposition */
    dgetrf_(&n, &n, colQ.data(), &lda, iw, &info);
    if (info < 0) fprintf(stderr, "LU failure with error %d\n", info);
    /* Computes the reciprocal norm */
    dgecon_("1", &n,colQ.data(), &lda, &anorm, &rcond, w, iw, &info);
    if (info < 0) fprintf(stderr, " con failure with error %d\n", info);
	double err = 0;
	//std::cout << "rcond" << "\t" << rcond << std::endl;
	if(rcond < 1e-19){ //1e-19
		//std::cout << "poorly conditioned" << std::endl;
		double mag = rcond*1e19; //1e19
	//	std::cout << "mag" << "\t" << mag << std::endl;
		err = log((1/mag))/log(10); //(/log(10))
		//std::cout << "err" << "\t" << err << std::endl;
	}
	

    return err; 
}

double model_penalty(Model& m){
  double rcond, anorm, penalty = 0;
  int N = m.G.N, E = m.G.E;


  int idmax, idmin;

  for (double vm = -80; vm <= 60; vm += 20) {
    m.transition_matrix(vm);
	//std::cout << "vm" << "\t" << vm << std::endl;
	//std::cout << "transition matrix" << std::endl;
	//m.print_Q();
	penalty+= calc_rcond(m);
	//std::cout << "penalty" << "\t" << penalty << std::endl;

  }
	
	penalty = penalty/8;
  
	//std::cout << "penalty" << penalty << std::endl;
  return penalty;
}

double calc_eigs(Model& m){
		int N = m.G.N;
		int LDA = N;
		int LDVL = N;
		int LDVR = N;
		
		 std::vector<double> colQ = transpose(m);
		 //m.print_Q();
	/* /* Locals */
        int n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR, info;
        /* Local arrays */
        double wr[N], wi[N], vl[LDVL*N], vr[LDVR*N];
       /*  std::vector<double> a({-1.01,  3.98,  3.30,  4.43,  7.31,
            0.86,  0.53,  8.26,  4.96, -6.43,
           -4.60, -7.04, -3.89, -7.66, -6.16,
            3.31,  5.29,  8.20, -7.33,  2.47,
           -4.81,  3.55, -1.51,  6.18,  5.58}); */
		   
        /* Executable statements */
       // printf( "LAPACKE_dgeev (column-major, high-level) Example Program Results\n" );
        /* Solve eigenproblem */
        info = LAPACKE_dgeev( LAPACK_COL_MAJOR, 'V', 'V', n, colQ.data(), lda, wr, wi,
                        vl, ldvl, vr, ldvr );
        /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                exit( 1 );
        }
        /* Print eigenvalues */
        double ratio = print_eigenvalues( "Eigenvalues", n, wr, wi );
        /* Print left eigenvectors */
        /* print_eigenvectors( "Left eigenvectors", n, wi, vl, ldvl );
        /* Print right eigenvectors */
        //print_eigenvectors( "Right eigenvectors", n, wi, vr, ldvr );
        
		
		//return ratio;
		if(ratio > 10){
			return ratio;
		}
		else{
			return 0;
		}
	
		
	
}


double print_eigenvalues( char* desc, MKL_INT n, double* wr, double* wi ) {
        MKL_INT j;
        //printf( "\n %s\n", desc );
		std::vector<double> abswr(n,0);
		std::vector<double> abswi(n,0);
   for( j = 0; j < n; j++ ) {
      if( wi[j] == (double)0.0 ) {
        // printf( " %6.2f", wr[j] );
		 abswr[j] = fabs(wr[j]);
      } else {
        // printf( " (%6.2f,%6.2f)", wr[j], wi[j] );
		 abswr[j] = sqrt(pow(wr[j],2) + pow(wi[j],2));
      }
   }
   //printf( "\n" );
   // disp(abswr);
	std::sort(abswr.begin(),abswr.end()); 
   //disp(abswr);
   int not_zero_indx = 0;
   for(int i = 0; i < n; i++){
		if(abswr[i] > 1e-5){
			not_zero_indx = i;
			break;
		}
	}
	//std::cout << not_zero_indx << std::endl;
   double ratio = abswr[n-1]/abswr[not_zero_indx];
   //ratio = log(ratio);
   //ratio = ratio/200;
   return ratio;
}

/* Auxiliary routine: printing eigenvectors */
void print_eigenvectors( char* desc, MKL_INT n, double* wi, double* v, MKL_INT ldv ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
   for( i = 0; i < n; i++ ) {
      j = 0;
      while( j < n ) {
         if( wi[j] == (double)0.0 ) {
            printf( " %6.2f", v[i+j*ldv] );
            j++;
         } else {
            printf( " (%6.2f,%6.2f)", v[i+j*ldv], v[i+(j+1)*ldv] );
            printf( " (%6.2f,%6.2f)", v[i+j*ldv], -v[i+(j+1)*ldv] );
            j += 2;
         }
      }
      printf( "\n" );
   }
}






////main cost function
double cost_main(Model& m, std::vector<ProtocolParameter> protos){
  m.cost = 0;
  //double val;
  Expm Qexp(m.G.N);
  for ( int i = 0; i < protos.size(); i++ ) {
	  //val =  (cost(m, protos[i],Qexp)*protos[i].weight);
    m.cost += (cost(m, protos[i], 0, Qexp)*protos[i].weight);
  }
  
 
  m.cost+=model_penalty(m);
  
  if(m.cost > 1e8) m.cost = 1e8;
  return  m.cost;
}

double cost_validation(Model& m, std::vector<ProtocolParameter> protos,std::vector<ProtocolParameter> valids){
  
  Expm Qexp(m.G.N);
  double val = 0;
  for(int i = 0; i < protos.size(); i++){
	  if(protos[i].has_validation_points){
		  
		  val +=  (cost(m, protos[i], 1, Qexp)*protos[i].weight);
		  
	  }
  }
  for ( int i = 0; i < valids.size(); i++ ) {
    val +=  (cost(m, valids[i], 0, Qexp)*valids[i].weight);
	//std::cout << val << std::endl;
  }
 
  return  val;
}

std::vector<double> cost_validation_comp(Model& m, std::vector<ProtocolParameter> protos,std::vector<ProtocolParameter> valids){
  std::vector<double> cost_val_vec;
  Expm Qexp(m.G.N);
  double val = 0;
  for(int i = 0; i < protos.size(); i++){
	  if(protos[i].has_validation_points){
		  
		  val =  (cost(m, protos[i], 1, Qexp)*protos[i].weight);
		   cost_val_vec.push_back(val);
	  }
	 
  }
  for ( int i = 0; i < valids.size(); i++ ) {
    val =  (cost(m, valids[i], 0, Qexp)*valids[i].weight);
	//std::cout << "val" << val << std::endl;
	cost_val_vec.push_back(val);
  }

  return  cost_val_vec;
}

std::vector<double> cost_comp(Model& m, std::vector<ProtocolParameter> protos){
  std::vector<double> cost_vec (protos.size());
  //std::cout << &m << std::endl;
  Expm Qexp(m.G.N);
  for ( int i = 0; i < protos.size(); i++ ) {
	  //std::cout << protos[i].get_name() << std::endl;
    cost_vec[i] =  (cost(m, protos[i], 0, Qexp)*protos[i].weight);
	//std::cout << "cost vector" << cost_vec[i] << std::endl;
  }
  return  cost_vec;
}

