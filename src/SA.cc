#include "SA.hpp"

SA::SA(std::string s,int time){ //set up paths for saving results;
	
	this->path = s;
	this->time = time;
	int size = path.size() - 2;
	while (path[size] != '/') size--;
	subpath = path.substr(0, size+1);
	size = subpath.size() - 2;
	while (subpath[size] != '/') size--;
	std::string sub2path = subpath.substr(0,size);
	size = sub2path.size() - 2;
	while (sub2path[size] != '/') size--;
	version = sub2path.substr(size+1); //version is used when saving to Amazon S3
}

SA::SA(std::string s,int time,Restart_params& restart_params){ ////set up paths for saving results after restart and reset overfitting tracking vectors
	
	this->path = s;
	this->time = time;
	this->cost_val_history = restart_params.cost_val_history;
	this->cost_tr_history = restart_params.cost_tr_history;
	this->PQ_history = restart_params.PQ_history;
	this->r = restart_params.r;
	int size = path.size() - 2;
	while (path[size] != '/') size--;
	subpath = path.substr(0, size+1);
	size = subpath.size() - 2;
	while (subpath[size] != '/') size--;
	std::string sub2path = subpath.substr(0,size);
	size = sub2path.size() - 2;
	while (sub2path[size] != '/') size--;
	version = sub2path.substr(size+1);
}

double SA::generalization_loss(){
	//how much the current validation error has increased since the overall minimum
	// see https://page.mi.fu-berlin.de/prechelt/Biblio/stop_tricks1997.pdf for the derivation
	double min_idx = min(cost_val_history);
	double current_idx = cost_val_history.size()-1;
	double GL = 100*((cost_val_history[current_idx]/cost_val_history[min_idx])-1);
	if( GL <= 1e-3) GL = 0;
	return GL;
}

double SA::progress(int k){ //int k is a display mod
	//how much progress has been made 
	// see https://page.mi.fu-berlin.de/prechelt/Biblio/stop_tricks1997.pdf for the derivation
	double p_k;
	
	if (cost_tr_history.size() < k) {
		return p_k = 1;
	}
	else{
		std::vector<double> v1(cost_tr_history.end()-k, cost_tr_history.end());
		double min_idx = min(v1);
		double sumv1 = sum(v1);
		p_k = 1000*((sumv1/(k*v1[min_idx]))-1);
		if( p_k <= 1) p_k = 1;
		return p_k;
	}
}


int SA::inc3(std::vector<double> v){ // needs to be size 4
	//we decided on 3 consecutive increases in the ratio of (PQ = progress/generalization) loss are needed for overfitting but feel free to change
	//for your problem!
	
	if(v[0] < v[1] && v[1] < v[2] && v[2] < v[3]){
	//std::cout << v[0] << "<" << v[1] << "<" << v[2] << "<" << v[3] << std::endl;		
		return 1;
		
	}
	else{
		return 0;
	}
	
	
}


int SA::detect_overfitting(){ 
	
	if (PQ_history.size() < 4) { //can't overfit, not enough history yet!
		
		return 0;
	}
	else{

		for(int i = 0; i < PQ_history.size()-3; i++){ //iterate over PQ history to pick out vectors of size 4 to look for 3 consecutive increases in PQ 
			std::vector<double> vtemp(PQ_history.begin()+i, PQ_history.begin()+i+3+1); 
			
			if(inc3(vtemp)) return 1; //if there have been three consecutive increases in PQ history, overfitting has occurred
		}
			return 0;
	}
}


int SA::S3_modelfiles(int N, std::string S3_bucket_path){ //call this function to save to your Amazon S3 bucket when running code in a containerized instance
	//the bucket path is passed in solver.txt and becomes part of the SimulatedParameters object
	pid_t pid;
	pid = fork();
	int status;
	//std::cout << "S3 bucket path" << S3_bucket_path  << std::endl;
	std::stringstream ss;
	ss << S3_bucket_path << "State" << N << "/" << version << "/";
	std::string path_dest = ss.str();
	//std::cout << "path_dest" << path_dest << std::endl;
	std::stringstream ss1;
	ss1 << "./State" << N << "/" << version << "/";
	std::string path_origin = ss1.str();
	//std::cout << "dir" << path_origin << std::endl;
	if(pid == 0){
		 execl("/usr/bin/aws", "aws", "s3","sync",path_origin.c_str(),path_dest.c_str(),(char*)NULL); //actual calling to S3
		_exit(127);
	}
	else if (pid > 0) {
		/* the parent process calls waitpid() on the child */
		if (waitpid(pid, &status, 0) > 0) {
			if (WIFEXITED(status) && !WEXITSTATUS(status)) {
				/* the program terminated normally and executed successfully */
				//std::cout << "everything good" << std::endl;
				return WEXITSTATUS(status);
			} 
			else if (WIFEXITED(status) && WEXITSTATUS(status)) {
				if (WEXITSTATUS(status) == 127) {
					/* execl() failed */
					std::cout << "execl failed" << std::endl;
					std::cout << "exit status" << (WEXITSTATUS(status)) << std::endl;
					return (WEXITSTATUS(status)); 
				} 
				else {
					/* the program terminated normally, but returned a non-zero status */
					std::cout << "the execl terminated normally, but S3 returned a non-zero status" << std::endl;
					std::cout << "exit status" << (WEXITSTATUS(status)) << std::endl;
					return (WEXITSTATUS(status)); 
					
				}
			} 
			else {
				/* the program didn't terminate normally */
				perror("the execl program didn't terminate normally");
			}
		} 
		else {
			/* waitpid() failed */
		}
	} 
	else {
		/* it was not possible to create child process, so print error message */
		perror("fork failed");
	}
	
	
	
	
}


void SA::print_display(std::ofstream& out, Model m, Setup s, Math math_params, int i,SimulationParameters& sim_params){
		//populate modelX.txt and the progress file iter_state.txt periodically throughout the optimization 
		
		const int display = sim_params.display;
		const int k_max = sim_params.k_max;
		out << i <<"\t" << fmin_val;
		//std::cout  << i << std::endl;
		auto costs_comps = cost_comp(fmin_model,s.protos);
		std::ostringstream ss;
		for(int k = 0; k < costs_comps.size(); k++){
			
			ss << costs_comps[k] << "\t";
		
		}
		std::string strVal= ss.str();
		out << "(" << strVal << ")" << std::endl;
	
		ss.str("");
		
		ss << subpath << "iter_state" << ".txt";
		std::string state_file = ss.str();
		std::ofstream state(state_file.c_str());
		
		
	
		cost_val_history.push_back(cost_validation(fmin_model,s.protos,s.valids));
		double GL = generalization_loss();
		//out << "GL" << GL << std::endl;
		//out << "validation error" << cost_val_history.back() << std::endl;
		double progressv = progress(display);
		//out << "progress" << progressv << std::endl;
		double PQ = GL/progressv;
		//out << "generalization loss /progress" << PQ << std::endl;
		PQ_history.push_back(PQ);
		state << "States:" << "\t" << m.G.N << std::endl;
		state << "Model:" << "\t" << m.get_id() << std::endl;
		state << "Time:" << "\t" << time << std::endl;
		state << "Iterations completed:" << "\t" << i << "\t" << "out of" << "\t" << k_max << std::endl;
		state << "Sobol Indx:" << "\t" << math_params.get_sobol_indx() << std::endl;
		state << "uni_Random Counter_uni_real:" << "\t" << math_params.get_random_counter_uni_real() << std::endl;
		//std::cout << "uni_Random Counter_uni_real:" << "\t" << math_params.get_random_counter_uni_real() << std::endl;
		state << "r vector:" << std::endl;
		print(r,state);
		//std::cout << " r vector" << std::endl;
		//disp(r);
		state << "Cost History" << std::endl;
		print(cost_tr_history,state);
		state << "Validate Vector" << std::endl;
		print(cost_val_history,state);
		state << "PQ History" << std::endl;
		print(PQ_history,state);
		state << "fmin_model" << std::endl;
		state << fmin_model << std::endl;
		
		state << fmin_val << std::endl;
		state << "Chains:" << "\t" << models.size() << std::endl;
		for(int k = 0; k < models.size(); k++){ //output models chain for restart
			state << models[k] << std::endl; 
			
		}
	
	
    S3_modelfiles(fmin_model.G.N,sim_params.AWS_S3_path);
	
}

void SA::print_snapshot(std::ofstream& out, Setup s, int i){

		//populate the progress files with model parameters and associated fits
		/* std::ostringstream iss; 
		std::string snapshot_file;
		iss << path << "iter_" << i << ".model";
		snapshot_file = iss.str();  */
	
		std::ostringstream mss; 
		std::string modelfit_file;
		mss << path << "iter_" << i << ".txt";
		modelfit_file = mss.str();

		 /* std::ofstream fss(snapshot_file.c_str());
		fss << (fmin_model) << std::endl; 
		fss.close();  */
		std::ofstream os(modelfit_file.c_str());
		os << (fmin_model) << std::endl;
				

		std::vector<double> modelval1;
		std::vector<double> data1;
		os << "TRAINING" << std::endl;
		os << "cost" << "\t" << fmin_val << std::endl;
		
		auto costs_comps = cost_comp(fmin_model,s.protos);
		for(int k = 0; k < costs_comps.size(); k++){
			
			os << costs_comps[k] << "\t";
		
		}
			os << std::endl;
			
		auto validation_cost_comp = cost_validation_comp(fmin_model,s.protos,s.valids);
		os << "Validation Cost"  << "\t" << sum(validation_cost_comp) << std::endl;
		for(int k = 0; k < validation_cost_comp.size(); k++){
		
			os << validation_cost_comp[k] << "\t";
	
		}
		os << std::endl;
			
		for (int j = 0; j < s.protos.size(); j++) {
			
			os << s.protos[j].get_name() << "\n";

			auto cost2 = cost((fmin_model), s.protos[j], 0, data1, modelval1);
			
			os << "Data" << "\t" << "Model" << "\n";
			for (int k = 0; k < data1.size(); k++) {

				os << data1[k] << "\t" << modelval1[k] << "\n" ;


			}
			data1.clear();
			modelval1.clear();
			if(s.protos[j].has_validation_points){
				os << "VALIDATION" << std::endl;
				auto cost2 = cost((fmin_model), s.protos[j], 1, data1, modelval1);
				
				os << "Data" << "\t" << "Model" << "\n";
				for (int k = 0; k < data1.size(); k++) {

					os << data1[k] << "\t" << modelval1[k] << "\n" ;


				}
			}

			data1.clear();
			modelval1.clear();
			os << std::endl;
		} 
		os << "VALIDATION" << std::endl;
		
		for (int j = 0; j < s.valids.size(); j++) {

			os << s.valids[j].get_name() << "\n";

			auto cost2 = cost((fmin_model), s.valids[j], 0, data1, modelval1);
			
			os << "Data" << "\t" << "Model" << "\n";
			for (int k = 0; k < data1.size(); k++) {

				os << data1[k] << "\t" << modelval1[k] << "\n" ;


			}
		data1.clear();
		modelval1.clear();
		os << std::endl;
		} 
		
		os.close(); 
			
		
	
}


void SA::anneal_restart(Math& math_params, Setup& s, SimulationParameters& sim_params, Restart_params& restart_params, std::ofstream& out){
	
	const int k_max = sim_params.k_max;
	const int step = sim_params.step;
	const int n_chains = sim_params.n_chains;
	const double T0 = sim_params.t0;
	const double gamma = sim_params.gamma;
	const int display = sim_params.display;
	const int snapshot = sim_params.snapshot;
	const std::string S3_bucket_path = sim_params.AWS_S3_path;
	const int N = s.get_N();

	f_val.resize(n_chains, 1e12);
	T_j.resize(n_chains);
	for(int i = 0; i < T_j.size(); i++){
		T_j[i] = T0 + (gamma*log(1+r[i]));
	}
	for (int i = 0; i < n_chains; i++){ //initialize the chain of models where the optimization stopped
		
		models.push_back(Model(s.get_model(),s.get_N(),s.edges[s.get_model()-1],s.roots[s.get_model()-1],restart_params.model_rs[i], restart_params.model_rk[i], restart_params.model_args[i],sim_params));
		f_val[i] = cost_main((models[i]),s.protos);
	}
	
	fmin_model = Model(s.get_model(),s.get_N(),s.edges[s.get_model()-1],s.roots[s.get_model()-1],restart_params.fmin_rs, restart_params.fmin_rk, restart_params.fmin_args, sim_params);
	std::cout << fmin_model << std::endl;
	fmin_val = cost_main(fmin_model,s.protos);
	 if(detect_overfitting()){
				
			std::cout << "overfitting detected" << std::endl;
			return;
				
	} 
	
	for (int i=restart_params.get_iterations()+1; i<= k_max; i++) {
	
		 if (i % step == 0) {
			disp(restart_params.r);
			disp(T_j);
		} 
		for ( int j=0; j<n_chains; j++ ) {
			
			Model candidate_params = models[j].perturb(math_params,sim_params);

			double fn2 = cost_main(candidate_params,s.protos); 
			double diff = (fn2 - f_val[j]);

			 
			
			if(diff > 1e-6){ // model is worse
				
				r[j]++;
				T_j[j] = T0 + (gamma*log(1+r[j]));
				
			}
			
			else if(diff <= 1e-6 && diff >= 0){
				//r[j] counter remains unchanged and temperature remains the same
				T_j[j] = T0 + (gamma*log(1+r[j]));
				
			}
			
			else{ //solution is better
				r[j] = 0; //reset counter of sequential worse costs
				
				T_j[j] = T0;

				
				
			}	

			 if (math_params.uni_dist() < exp(-(fn2 - f_val[j])/T_j[j]))  {
				//std::cout << "new model accepted" << std::endl;
				
				f_val[j] = fn2;
				models[j] = candidate_params;
			} 

		
		 
		}//n_chains
		
		for ( int j=0; j<n_chains; j++ ) {//if a cost in the chains is lower than current min, update!
			if ( f_val[j] < fmin_val ) {
				fmin_val = f_val[j];
				fmin_model = models[j];
			}
		}
		cost_tr_history.push_back(fmin_val);
		
		
		if (i % snapshot == 0 ) { //snapshot needs be multiple of display mods
				
			//S3_modelfiles(N, S3_bucket_path);
			print_snapshot(out,s,i);
			
		} 
		
		
		if (i % display == 0){
			
			print_display(out,fmin_model,s,math_params,i,sim_params);
			
			 if(detect_overfitting()){
				
				std::cout << "overfitting detected" << std::endl;
				break;
				
			} 			
		}
		
	}
	
	
}
void SA::anneal(Math& math_params, Model& m, Setup& s, SimulationParameters& sim_params, std::ofstream& out){
	
	const int k_max = sim_params.k_max;
	const int step = sim_params.step;
	const int n_chains = sim_params.n_chains;
	const double T0 = sim_params.t0;
	const double gamma = sim_params.gamma;
	const int display = sim_params.display;
	const int snapshot = sim_params.snapshot;
	const std::string S3_bucket_path = sim_params.AWS_S3_path;
	int N =  m.G.N;
	
	

	models.resize(n_chains,m);
	f_val.resize(n_chains, 1e12);
	r.resize(n_chains,0);
	T_j.resize(n_chains,T0);
	for (int i = 0; i < n_chains; i++){ //initialize a chain of new models

		Model n(m);
		f_val[i] = cost_main(n,s.protos);
	}
	
	
	
	fmin_model = models[0]; //setup our fmins
	fmin_val = f_val[0];
	
	
	for (int i=0; i<= k_max; i++){
		
		
		/* if (i % step == 0) {
			//std::cout << "i" << i << std::endl;		
			disp(r);
			disp(T_j);
		} */
		
		
		for ( int j=0; j<n_chains; j++ ) {
			
			Model candidate_params = models[j].perturb(math_params,sim_params);
			double fn2 = cost_main(candidate_params,s.protos); 
			double diff = (fn2 - f_val[j]);

			 
			
			if(diff > 1e-6){ // model is worse
				
				r[j]++;
				T_j[j] = T0 + (gamma*log(1+r[j]));
				//disp(r);
			}
			
			else if(diff <= 1e-6 && diff >= 0){
				
				//r[j] counter remains unchanged and temperature remains the same
				T_j[j] = T0 + (gamma*log(1+r[j]));
				//disp(r);
			}
			
			else{ //solution is better
				//std::cout << "solution better" << std::endl;
				r[j] = 0; //reset counter of sequential worse costs
				
				T_j[j] = T0;
				
				 //disp(r);
				
			}	

			 if (math_params.uni_dist() < exp(-(fn2 - f_val[j])/T_j[j]))  {
				
				f_val[j] = fn2;
				models[j] = candidate_params;
			} 

		
		 
		}//n_chains
		
		
		for ( int j=0; j<n_chains; j++ ) {//if a cost in the chains is lower than current min, update!
			if ( f_val[j] < fmin_val ) {
				fmin_val = f_val[j];
				fmin_model = models[j];
			}
		}
	
	
	
		cost_tr_history.push_back(fmin_val); //save cost history for calculating progress
		
		
		if (i % snapshot == 0 ) { //snapshot needs be multiple of display mods
				
			//S3_modelfiles(N, S3_bucket_path);
			print_snapshot(out,s,i);
			
		} //snapshot
		
		if (i % display == 0){
			std::cout << i << "\t" << fmin_val << std::endl;
			print_display(out,m,s,math_params,i,sim_params);
			
			 if(detect_overfitting()){
				
				std::cout << "overfitting detected" << std::endl;
				break;
				
			} 
			
			
			
		}//display
		
	
	}// i iter
}
