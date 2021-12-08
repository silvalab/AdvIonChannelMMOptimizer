#include "SA.hpp"
#include <sstream>
#include <omp.h>
#include <memory>
#include <ctime>
#include <thread>
#include "same.hpp"
#include "cluster.hpp"


template <class T>
void print_vec_vec_vec(std::vector<std::vector<T>> v){
	
	for (int i = 0; i < v.size(); i++){
		for (int j = 0; j < v[i].size(); j++){
			for(int k = 0; k < v[i][j].size(); k++){
				std::cout << v[i][j][k] << "\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}



template <class T>
void print_vec_vec(std::vector<std::vector<T>> v){
	
	for (int i = 0; i < v.size(); i++){
		
		for (int j = 0; j < v[i].size(); j++){
			
			std::cout << v[i][j] << "\t";
			
		}
		std::cout << std::endl;
	}
}

template <class T>
void print_vec(std::vector<T> v){
	
	for (int i = 0; i < v.size(); i++){
		
		std::cout << v[i] << "\t";
		
	}
	std::cout << std::endl;
}

template<class T>
double sum_vector(std::vector<T> costs){
	
	double sum = 0;
	for(int i = 0; i < costs.size(); i++){
		
		sum = sum + costs[i];
	}
	
	return sum;
}


template <class T>
int vec_min(std::vector<T> v){
	int min_idx = 0;
	T min = v[0];
	for(int i = 1; i < v.size(); i++){
		if( v[i] < min){
			min = v[i];
			min_idx = i;
		}
	}
	return min_idx;
}

int domination(Model* current_model, Model* new_model, std::vector<ProtocolParameter> protos){
	
	std::vector<double> costs_current = cost_comp(*current_model,protos);
	std::vector<double> costs_new = cost_comp(*new_model,protos);
	
	//current dominates new point: return 0
	//current and new are nondominating to each other: return 1
	// new dominates current: return 2
	
	double cost = sum_vector(costs_current);
	//std::cout << "cost current" << cost << std::endl;
	double cost1 = sum_vector(costs_new);
	//std::cout << "cost new" << cost1 << std::endl;
	int better = 0; 
	int worse = 0;
	int equal = 0;
	if( cost <= cost1){
		//int and better are from the perspective of current dominating
		//check for dominating or non 
		for(int i = 0; i < costs_current.size(); i++){
			if(costs_current[i] < costs_new[i]){
				better++;
			}
			else if(costs_current[i] > costs_new[i]){
				worse++;
				
			}
			else{
				equal++;
			}
		}
		
		if(better >= 1 && better+equal == costs_current.size()){

			return 0; //current dominates new
		}
		else{
			return 1; //current and new are nondominating
		}
	}
	
	else{ //cost1 < cost 
		//int and better are from the perspective of new dominating
		//check for dominating or non 
		for(int i = 0; i < costs_current.size(); i++){
			if(costs_new[i] < costs_current[i]){
				better++;
			}
			else if(costs_new[i] > costs_current[i]){
				worse++;
				
			}
			else{
				equal++;
			}
		}
		
		if(better >= 1 && better+equal == costs_current.size()){

			return 2; //current dominates new
		}
		else{
			return 1; //current and new are nondominating
		}
		
		
	}
	
	
	
	
}

int domination(Model* current_model, Model* new_model, std::vector<ProtocolParameter> protos, std::vector<double>& costs_current, std::vector<double>& costs_new){
	
	costs_current = cost_comp(*current_model,protos);
	costs_new = cost_comp(*new_model,protos);
	
	//current dominates new point: return 0
	//current and new are nondominating to each other: return 1
	// new dominates current: return 2
	
	double cost = sum_vector(costs_current);
	//std::cout << "cost current" << cost << std::endl;
	double cost1 = sum_vector(costs_new);
	//std::cout << "cost new" << cost1 << std::endl;
	int better = 0; 
	int worse = 0;
	int equal = 0;
	if( cost <= cost1){
		//int and better are from the perspective of current dominating
		//check for dominating or non 
		for(int i = 0; i < costs_current.size(); i++){
			if(costs_current[i] < costs_new[i]){
				better++;
			}
			else if(costs_current[i] > costs_new[i]){
				worse++;
				
			}
			else{
				equal++;
			}
		}
		
		if(better >= 1 && better+equal == costs_current.size()){

			return 0; //current dominates new
		}
		else{
			return 1; //current and new are nondominating
		}
	}
	
	else{ //cost1 < cost 
		//int and better are from the perspective of new dominating
		//check for dominating or non 
		for(int i = 0; i < costs_current.size(); i++){
			if(costs_new[i] < costs_current[i]){
				better++;
			}
			else if(costs_new[i] > costs_current[i]){
				worse++;
				
			}
			else{
				equal++;
			}
		}
		
		if(better >= 1 && better+equal == costs_current.size()){

			return 2; //current dominates new
		}
		else{
			return 1; //current and new are nondominating
		}
		
		
	}
	
}



double delta_dom(std::vector<double> costs_current,std::vector<double> costs_new,int proto_num){
	std::vector<double> Ranges;
	if(proto_num == 6) // needs to be updated manually for each list of voltage protocols!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		Ranges = {192,47,1.5,119.16, 41.04, 15.28};
	else if(proto_num == 5) 
		Ranges = {192,47,1.5, 119.16, 41.04};
	else if(proto_num == 4) 
		Ranges = {192,47,1.5, 119.16};
	else if(proto_num == 3) 
		Ranges = {192,47,1.5};
	else if(proto_num == 2) 
		Ranges = {192, 47};
	else if(proto_num == 1) 
		Ranges = {192};
	
	double dom =1;
	for (int i = 0; i < Ranges.size(); i++){
		
		float d = std::abs(costs_current[i]-costs_new[i]);
		d = d/Ranges[i];
		dom = dom*d;
	}
	//std::cout << dom << std::endl;
	return dom;
	
}

int contains(std::vector<int> k_vec, int x){
	
	//check to see if there is at least one 2 in k_vec
	for(int i = 0; i < k_vec.size(); i++){
		if(k_vec[i] == x) return 1;
		
	}
	return 0;
}

std::vector<int> find_indxs(std::vector<int> k_vec,int x){
	std::vector<int> indxs;
	for(int i = 0; i < k_vec.size(); i++){
		if(k_vec[i] == x)
		indxs.push_back(i);
	}
	return indxs;
}


int szof(std::vector<std::vector<Model*>> models){
	int sum = 0;
	
	for (int i = 0; i < models.size(); i++){
		
		sum += models[i].size();
	}
	return sum;
}


int min_cost_idx(std::vector<std::vector<double>> f_val){ //find which j has the best archive so far
	
	double min = f_val[0][0];
	int idx = 0;
	for (int i = 0; i < f_val.size(); i++){
		for (int j = 0; j < f_val[i].size(); j++){		
			if(f_val[i][j] < min){
				min = f_val[i][j];
				idx = i;
			}
		}
	}
	return  idx;
}

double SAmosa(std::string path,Model* m, std::vector<ProtocolParameter> protos, SimulationParameters sim_params, std::ofstream& out){


	const int k_max = sim_params.k_max;
	const int step = sim_params.step;
	const int n_chains = sim_params.n_chains;
	const double T0 = sim_params.t0;
	const double gamma = sim_params.gamma;
	const int display = sim_params.display;
	const int snapshot = sim_params.snapshot;
	const std::string snapshotdir = sim_params.snapshotdir;
	int N =  m->G.N;
	//omp_set_dynamic(0);

	//omp_set_num_threads(numThreads);

	int archive_size_initial = 15; //needs to be y*HL
	int archive_iter_initial = 25;
	int SL = 15;
	int HL = 10;
	
	
	
	std::vector<std::vector<Model*>> models (n_chains, std::vector<Model*> (archive_size_initial));
	std::vector<std::vector<double>> f_val (n_chains, std::vector<double> (archive_size_initial)); //f_val is not changing
	std::vector<std::vector<std::vector<double>>> model_vals(n_chains, std::vector<std::vector<double>> (archive_size_initial));
	
	//#pragma omp parallel for
	for (int i = 0; i < n_chains; i++){
		for (int j = 0; j < archive_size_initial; j++){ //initialize a chain of new models
			models[i][j] =  new_candidate_params_only(m,sim_params);
			//std::cout << models[i] << std::endl;
			f_val[i][j] = cost_main((*models[i][j]),protos);	
			//std::cout << "model inside of complicated SA" << std::endl;
			//std::cout << "models i" << (models[i]) << std::endl;
			//model_vals[i].push_back(cost_comp((*models[i][j]),protos));
			//std::cout << "starting cost" << f_val[i] << std::endl;
		}
	}
	
	//print_vec_vec_vec(model_vals);
	
	
	//implememt bottom sliding to initialize archive
	//only accept models that completely dominate 
	#pragma omp parallel for
	for (int i = 0; i < n_chains; i++){
		for(int j = 0; j < archive_size_initial; j++){
			for (int k = 0; k < archive_iter_initial; k++){
			Model* model_candidate =  new_candidate_params_only(m,sim_params);
			std::vector<double> costs = cost_comp(*model_candidate,protos);
			double total_cost = sum_vector(costs);
			int better = 0;
				for (int l = 0; l < model_vals[i][j].size(); l++){//iterate over protocols
					if(costs[l] <= model_vals[i][j][l]) better++;
				}
				
				if( better == model_vals[i][j].size()){ //solution is better for all protos
					delete models[i][j];
					models[i][j] = new Model(model_candidate);
					f_val[i][j] = total_cost;
					delete model_candidate;
				}
				else{
				
				delete model_candidate;
				
				}
			}
		}
	}	
	//print_vec_vec(f_val); 
		
	//cluster here	
	//bring archive down to HL	
	for (int i = 0; i < n_chains; i++){
		if (models[i].size() > SL) {
			single_linkage_cluster(models[i], f_val[i], HL, protos);
		} 	
	}
		
	std::vector<double> best_pareto_so_far = f_val[0];
	print_vec(best_pareto_so_far);
	double mini = best_pareto_so_far[vec_min(best_pareto_so_far)];
	//std::cout << "mini of f_val zero" << "\t" << mini << std::endl;
	std::vector<Model*> best_paretom_so_far (best_pareto_so_far.size()); 
	for (int i = 0; i < best_paretom_so_far.size(); i++){
		best_paretom_so_far[i] = new Model(models[0][i]);
	}
	
	
	std::vector<int> a_indices(n_chains);
	for (int i =0; i < n_chains; i++){
		a_indices[i] = Math::rng_int(0,models[i].size());
	}
	double T_j = T0;
	for (int i = 0; i <= k_max; i+=display){
		
			//std::cout << "archive size at iteration:" << "\t" << i << "\t" << szof(models) << std::endl;
			
			//std::cout << "best pareto front seen so far" << std::endl; 
			// update best pareto;
			int min_chain = min_cost_idx(f_val);
			//std::cout << "min_chain" << min_chain << std::endl;
			if(f_val[min_chain][vec_min(f_val[min_chain])] < mini){
				mini = f_val[min_chain][vec_min(f_val[min_chain])];
				//std::cout << "new mini" << mini << std::endl;
				best_pareto_so_far.clear();
				best_pareto_so_far = f_val[min_chain];
				for(int k = 0; k < best_paretom_so_far.size(); k++){
					delete(best_paretom_so_far[k]);
				}
				best_paretom_so_far.resize(models[min_chain].size());
				for(int k = 0; k < best_paretom_so_far.size(); k++){
					best_paretom_so_far[k] = new Model(models[min_chain][k]);
				}
				
			}
			/* if(mini < 0.04){
				std::cout << "Good pareto front found at iteration" << "\t" << i << "\t" << "with cost of" << "\t" << mini << std::endl;
				out << "Good pareto front found at iteration" << "\t" << i << "\t" << "with cost of" << "\t" << mini << std::endl;
				print_vec(best_pareto_so_far);
				break;
			} */
			//print_vec(best_pareto_so_far);
			//std::cout << "current temp together" << "\t" << T_j << std::endl;
		
		//print_vec_vec(models);
		//print_vec_vec(f_val);
		//print_vec_vec(models);
		/* for (int j =0; j < n_chains; j++){
			std::cout << a_indices[j] << std::endl;
		}  */
			
		#pragma omp parallel for
		for (int j = 0; j < n_chains; j++){
			T_j=T0;
			if (i % step == 0){
				for (int s = 0; s < i/step+1; s++){ 
					T_j*= gamma;
					//std::cout << "T-J gamma" << "\t" << T_j << std::endl;
					
				}
			}
			//std::cout << "current temp" << "\t" << T_j << std::endl;
			for (int z = 0 ; z< display; z++){
				//std::cout << "current temp " << "\t" << T_j << std::endl;
				Model* current_model = models[j][a_indices[j]]; 
				//std::cout << current_model << std::endl;
				Model* new_model =  new_candidate_params_only(current_model,sim_params);
				//std::cout << new_model << std::endl;
				std::vector<double> costs_current (protos.size());
				std::vector<double> costs_new (protos.size());
		
	
				int dom = domination(current_model,new_model,protos,costs_current,costs_new);
				//std::cout << "costs current" << std::endl;
				//print_vec(costs_current);
				//std::cout << "costs new" << std::endl;
				//print_vec(costs_new);
				double domq = delta_dom(costs_current,costs_new,(protos.size()));
			
	
				//std::cout << "dom" << dom << std::endl;
				//std::cout << "amount of domination" << domq << std::endl;
				if(dom == 0){ // current dominates new
				//std::cout << "dom" << dom << std::endl;
					//find k: no of points in archive which dominate new pt
					std::vector<int> k_vec (models[j].size());
					 for(int k = 0; k < models[j].size();k++){ //look over archive
						k_vec[k] = domination(new_model,models[j][k],protos);
					} 
					//std::cout << "k_vec" << std::endl;
					//print_vec(k_vec);
					double dom_archive_total = 0;
					auto vec = find_indxs(k_vec,2);
					for(int k = 0; k < vec.size(); k++){
				
						std::vector<double> costs_current_archive_temp (protos.size());
						std::vector<double> costs_new_archive_temp (protos.size());
				
						int dom0 = domination(models[j][vec[k]],new_model,protos,costs_current_archive_temp,costs_new_archive_temp);
						double domq0 = delta_dom(costs_current_archive_temp,costs_new_archive_temp,(protos.size()));
						dom_archive_total = dom_archive_total+domq0;
					} 
			
					double dom_avg = (dom_archive_total+ domq)/(vec.size()+1);
					double prob = 1/(1+exp(dom_avg*T_j)); 
					//double prob = exp(-(dom_avg)/T_j);
					//std::cout << "prob in dom 0" << "\t" << prob << std::endl; */
					if( prob < Math::rng_uniform(0,1)){
						//std::cout << "new model accepted that is worse" << std::endl;
						//std::cout << "models before delte" << std::endl;
						//print_vec_vec(models);
						delete (models[j][a_indices[j]]);
						//std::cout << "j" << j << "a_indices[j]" << a_indices[j] << std::endl;
						//std::cout << "models after delte" << std::endl;
						//print_vec_vec(models);
						models[j][a_indices[j]] = new Model(new_model);
						//std::cout << "models after replace new model" << std::endl;
						//print_vec_vec(models);
						f_val[j][a_indices[j]] = sum_vector(costs_new);
						delete (new_model);
					}
					else{
						delete (new_model);
					}
				}
		
		
		
				else if(dom ==1){ //current and new are ND 
					//std::cout << "dom" << dom << std::endl;
					//populate k_vec
					std::vector<int> k_vec (models[j].size());
					for(int k = 0; k < models[j].size();k++){ //look over archive
						k_vec[k] = domination(new_model,models[j][k],protos);	
					} 
			
					//std::cout << "k_vec" << std::endl;
					//print_vec(k_vec);
					auto vec = find_indxs(k_vec,2);
					auto vec1 = find_indxs(k_vec,1);
					if(vec.size() >= 1){  
						//case 2A //k_vec contains at least one 2
						double dom_archive_total = 0;
						for(int k = 0; k < vec.size(); k++){
				
							std::vector<double> costs_current_archive_temp (protos.size());
							std::vector<double> costs_new_archive_temp (protos.size());
				
							int dom1 = domination(models[j][vec[k]],new_model,protos,costs_current_archive_temp,costs_new_archive_temp);
							double domq1 = delta_dom(costs_current_archive_temp,costs_new_archive_temp,(protos.size()));
							dom_archive_total = dom_archive_total+domq1;
						}
						double dom_avg = (dom_archive_total+ domq)/(vec.size());
						double prob = 1/(1+exp(dom_avg/T_j)); 
						//	std::cout << "prob in dom 1" << "\t" << prob << std::endl;
						if( prob < Math::rng_uniform(0,1)){
							//std::cout << "new model accepted that is nondom" << std::endl;
							//std::cout << "models before delte" << std::endl;
							//print_vec_vec(models);
							delete (models[j][a_indices[j]]);
							//std::cout << "j" << j << "a_indices[j]" << a_indices[j] << std::endl;
							//std::cout << "models after delte" << std::endl;
							//print_vec_vec(models);
							models[j][a_indices[j]] = new Model(new_model);
							//std::cout << "models after replace new model" << std::endl;
							//print_vec_vec(models);
							f_val[j][a_indices[j]] = sum_vector(costs_new);
							delete (new_model);
						}
						else{
							delete (new_model);
						}
					
					}
					else if (vec1.size() == models[j].size()){
						// case 2B //new point is ND to all of archive
						//std::cout << "K-vec == 1" << std::endl;
						//std::cout << "model push back" << std::endl;
						
						//std::cout << "models after add" << std::endl;
						//print_vec_vec(models);
						if (models[j].size() > SL) {
							single_linkage_cluster(models[j], f_val[j], HL, protos);
						}  
						models[j].push_back(new Model(new_model));
						a_indices[j] = models[j].size()-1;
						f_val[j].push_back(sum_vector(costs_new));
						delete new_model;
						//std::cout << "archive_indx updated" << a_indices[j]  << std::endl;
						//std::cout << "model added" << std::endl;
						 
					}
				
					else{ // case 2C new-pt dominates k >= 1 points in archive 
				
						//std::cout << "k vec has at least one zero" << std::endl;
						//print_vec(k_vec);
						std::vector<Model*> models_temp_new;
						std::vector<double> fvals_temp;
						int count = 0;
						//std::cout << "models before mass delte" << std::endl;
						//print_vec_vec(models);
						for (int k = 0; k < models[j].size(); k++){
							if(k_vec[k] != 0){ 
								models_temp_new.push_back(models[j][k]);
								fvals_temp.push_back(f_val[j][k]);
							}
							else{
								delete models[j][k];
								count++;
							}
						}			
						models[j].swap(models_temp_new);
						f_val[j].swap(fvals_temp);
						//std::cout << "models after mass delte" << std::endl;
						//print_vec_vec(models);
						//std::cout << count << "models deleted" << std::endl;
						//set current point as new point and add to archive
						models[j].push_back(new Model(new_model));
						delete (new_model);
						a_indices[j] = models[j].size()-1;
						f_val[j].push_back(sum_vector(costs_new));
						
						//std::cout << "model added" << std::endl;
						//std::cout << "archive_indx updated" << a_indices[j]  << std::endl
						
					}   
				}   
		
				else{ // Case 3 new dominates current
					//std::cout << "dom" << dom << std::endl;
					//populate k_vec
					std::vector<int> k_vec (models[j].size());
					for(int k = 0; k < models[j].size();k++){ //look over archive
							k_vec[k] = domination(new_model,models[j][k],protos);	
						}
						auto vec = find_indxs(k_vec,2);
						auto vec1 = find_indxs(k_vec,1); 
					if(vec.size() >= 1){ 
						//case 3A //k_vec contains at least one 2
						std::vector<double> dom_archive (vec.size());
						for(int k = 0; k < vec.size(); k++){
				
							std::vector<double> costs_current_archive_temp (protos.size());
							std::vector<double> costs_new_archive_temp (protos.size());
				
							int dom2 = domination(models[j][vec[k]],new_model,protos,costs_current_archive_temp,costs_new_archive_temp);
							double domq2 = delta_dom(costs_current_archive_temp,costs_new_archive_temp,(protos.size()));
							dom_archive[k] = domq2;
						}
						int min_idx = vec_min(dom_archive);
						double prob = 1/(1+exp(-(dom_archive[min_idx]))); 
						//std::cout << "prob in dom 2" << "\t" << prob << std::endl; */
						if( prob < Math::rng_uniform(0,1)){
							//std::cout << "new dominates current, pointer updated to min" << std::endl;
							a_indices[j] = vec[min_idx];
							//std::cout << "archive_indx updated" << archive_indx << std::endl;
							delete (new_model);
						}
						else{
							delete (models[j][a_indices[j]]);
							models[j][a_indices[j]] = new Model(new_model);
							f_val[j][a_indices[j]] = sum_vector(costs_new);
							delete (new_model);
						}
					}
			
					else if(vec1.size() == models[j].size()){	
						if (models[j].size() > SL) {
							single_linkage_cluster(models[j], f_val[j], HL, protos);
						}  
						models[j].push_back(new Model(new_model));
						delete (new_model);
						a_indices[j] = models[j].size()-1;
						f_val[j].push_back(sum_vector(costs_new));
					}
					else{ // new_point dominates k other points in archive
						int count = 0;
						std::vector<Model*> models_temp_new;
						std::vector<double> fvals_temp;
						for (int k = 0; k < models[j].size(); k++){
							if(k_vec[k] != 0){ 
								models_temp_new.push_back(models[j][k]);
								fvals_temp.push_back(f_val[j][k]);
							}
							else{
								delete models[j][k];
								count++;
							}
						}
						models[j].swap(models_temp_new);
						f_val[j].swap(fvals_temp);
						//std::cout << count << "models deleted" << std::endl;
						models[j].push_back(new Model(new_model));
						delete (new_model);
						a_indices[j] = models[j].size()-1;
						f_val[j].push_back(sum_vector(costs_new));
						//std::cout << "archive_indx updated" << archive_indx << std::endl;
						//std::cout << "model added" << std::endl;
					} 
				} //dominate cases */
			} //z times per display
		}//n_chains
	} // end_i	 
	//std::cout << "come together now" << std::endl;
	//print_vec_vec(models);
	//print_vec_vec(f_val);
	
	int min_chain = min_cost_idx(f_val);
	//std::cout << "min_chain" << min_chain << std::endl;
	if(f_val[min_chain][vec_min(f_val[min_chain])] < mini){
		mini = f_val[min_chain][vec_min(f_val[min_chain])];
		//std::cout << "new mini" << mini << std::endl;
		best_pareto_so_far.clear();
		best_pareto_so_far = f_val[min_chain];
		for(int k = 0; k < best_paretom_so_far.size(); k++){
			delete(best_paretom_so_far[k]);
		}
		best_paretom_so_far.resize(models[min_chain].size());
		for(int k = 0; k < best_paretom_so_far.size(); k++){
			best_paretom_so_far[k] = new Model(models[min_chain][k]);
		}
				
	}
	print_vec(best_pareto_so_far);
	print_vec(best_paretom_so_far); 
	
	
	//cluster 
	if (best_paretom_so_far.size() > HL) single_linkage_cluster(best_paretom_so_far,best_pareto_so_far,HL, protos);
	
	std::vector<double> min_costs (best_paretom_so_far.size(), 0);
	for(int i = 0; i < best_paretom_so_far.size(); i++){
		
		min_costs[i] = cost_main(*best_paretom_so_far[i],protos);
		auto costs = cost_comp((*best_paretom_so_far[i]),protos);	
		std::stringstream ss;
		for(int k = 0; k < costs.size(); k++){
					
			ss << costs[k] << "\t";
				
		}
		std::string strVal= ss.str();
		out << i+1 << ":" << "\t" << min_costs[i] << "\t" << "(" << strVal << ")" << std::endl;
		std::ostringstream iss; 
		std::string snapshot_file;
		iss << path << i+1 << ".model";
		snapshot_file = iss.str();

		std::ostringstream mss; 
		std::string modelfit_file;
		mss << path << i+1 << ".txt";
		modelfit_file = mss.str();

		std::ofstream fss(snapshot_file.c_str());
		fss << (*best_paretom_so_far[i]) << std::endl;
		fss.close();
		std::ofstream os(modelfit_file.c_str());

		std::vector<double> modelval1;
		std::vector<double> data1;
		os << "cost" << min_costs[i] << std::endl;
		
		for(int k = 0; k < costs.size(); k++){
					
			os << costs[k] << "\t";
				
		}
			os << std::endl;
		for (int k = 0; k < protos.size(); k++) {
		
			os << protos[k].name << "\n";
		
			cost((*best_paretom_so_far[i]), protos[k], data1, modelval1);
			os << "Data" << "\t" << "Model" << "\n";
			for (int l = 0; l < protos[k].steps.size(); l++) {

				os << data1[l] << "\t" << modelval1[l] << "\n" ;


			}

			data1.clear();
			modelval1.clear();
			os << std::endl;
			} 
	  
		os.close(); 
		delete best_paretom_so_far[i];
	}  
	
	
	for(int i = 0; i < models.size(); i++){
		for (int j = 0; j < models[i].size(); j++){
			delete models[i][j];
		}
	}
	int min_idx_to_return = vec_min(min_costs);
	return min_costs[min_idx_to_return]; 
} 

/* double SAmosa(std::string path,Model::Model* m, std::vector<ProtocolParameter> protos, SimulationParameters sim_params, std::ofstream& out){


	const int k_max = sim_params.k_max;
	const int step = sim_params.step;
	const int n_chains = sim_params.n_chains;
	const double T0 = sim_params.t0;
	const double gamma = sim_params.gamma;
	const int display = sim_params.display;
	const int snapshot = sim_params.snapshot;
	const std::string snapshotdir = sim_params.snapshotdir;
	int N =  m->G.N;
	//omp_set_dynamic(0);

	//omp_set_num_threads(numThreads);

	int archive_size_initial = 25; //needs to be y*HL
	int archive_iter_initial = 5;
	int SL = 20;
	int HL = 10;
	
	
	
	std::vector<std::vector<Model*>> models (n_chains, std::vector<Model*> (archive_size_initial));
	std::vector<std::vector<double>> f_val (n_chains, std::vector<double> (archive_size_initial)); //f_val is not changing
	std::vector<std::vector<std::vector<double>>> model_vals(n_chains, std::vector<std::vector<double>> (archive_size_initial));
	
	//#pragma omp parallel for
	for (int i = 0; i < n_chains; i++){
		for (int j = 0; j < archive_size_initial; j++){ //initialize a chain of new models
			models[i][j] =  new_candidate_params_only(m,sim_params);
			//std::cout << models[i] << std::endl;
			//f_val[i][j] = cost_main((*models[i][j]),protos);
			f_val[i][j] = Math::rng_uniform(0,10);		
			//std::cout << "model inside of complicated SA" << std::endl;
			//std::cout << "models i" << (models[i]) << std::endl;
			//model_vals[i].push_back(cost_comp((*models[i][j]),protos));
			//std::cout << "starting cost" << f_val[i] << std::endl;
		}
	}
	
	//print_vec_vec_vec(model_vals);
	
	
	//implememt bottom sliding to initialize archive
	//only accept models that completely dominate 
	/* #pragma omp parallel for
	for (int i = 0; i < n_chains; i++){
		for(int j = 0; j < archive_size_initial; j++){
			for (int k = 0; k < archive_iter_initial; k++){
			Model* model_candidate =  new_candidate_params_only(m,sim_params);
			std::vector<double> costs = cost_comp(*model_candidate,protos);
			double total_cost = sum_vector(costs);
			int better = 0;
				for (int l = 0; l < model_vals[i][j].size(); l++){//iterate over protocols
					if(costs[l] <= model_vals[i][j][l]) better++;
				}
				
				if( better == model_vals[i][j].size()){ //solution is better for all protos
					delete models[i][j];
					models[i][j] = new Model(model_candidate);
					f_val[i][j] = total_cost;
					delete model_candidate;
				}
				else{
				
				delete model_candidate;
				
				}
			}
		}
	}	 
	print_vec_vec(f_val); 
		
	//cluster here	
	//bring archive down to HL	
	for (int i = 0; i < n_chains; i++){
		if (models[i].size() > SL) {
			single_linkage_cluster(models[i], f_val[i], HL, protos);
		} 	
	}
		
	std::vector<double> best_pareto_so_far = f_val[0];
	print_vec(best_pareto_so_far);
	double mini = best_pareto_so_far[vec_min(best_pareto_so_far)];
	std::cout << "mini of f_val zero" << "\t" << mini << std::endl;
	std::vector<Model*> best_paretom_so_far (best_pareto_so_far.size()); 
	for (int i = 0; i < best_paretom_so_far.size(); i++){
		best_paretom_so_far[i] = new Model(models[0][i]);
		
	}
	
	
	
	
	
	std::vector<int> a_indices(n_chains);
	for (int i =0; i < n_chains; i++){
		a_indices[i] = Math::rng_int(0,models[i].size());
		//a_indices[i] = models[i].size()-1;
		//std::cout << a_indices[i] << std::endl;
	}
	double T_j = T0;
	for (int i = 0; i <= k_max; i+=display){
		
			std::cout << "archive size at iteration:" << "\t" << i << "\t" << szof(models) << std::endl;
			
			std::cout << "best pareto front seen so far" << std::endl; 
			// update best pareto;
			int min_chain = min_cost_idx(f_val);
			std::cout << "min_chain" << min_chain << std::endl;
			if(f_val[min_chain][vec_min(f_val[min_chain])] < mini){
				mini = f_val[min_chain][vec_min(f_val[min_chain])];
				std::cout << "new mini" << mini << std::endl;
				best_pareto_so_far.clear();
				best_pareto_so_far = f_val[min_chain];
				for(int k = 0; k < best_paretom_so_far.size(); k++){
					delete(best_paretom_so_far[k]);
				}
				best_paretom_so_far.resize(models[min_chain].size());
				for(int k = 0; k < best_paretom_so_far.size(); k++){
					best_paretom_so_far[k] = new Model(models[min_chain][k]);
				}
				
			}
			//if(mini <= 0.05){
				//std::cout << "Good pareto front found at iteration" << "\t" << i << "\t" << "with cost of" << "\t" << mini << std::endl;
				//print_vec(best_pareto_so_far);
				//break;
			//}
			print_vec(best_pareto_so_far);
			//std::cout << "current temp together" << "\t" << T_j << std::endl;
		
		//print_vec_vec(models);
		//print_vec_vec(f_val);
		//print_vec_vec(models);
		for (int j =0; j < n_chains; j++){
			std::cout << a_indices[j] << std::endl;
		} 
			
		//#pragma omp parallel for
		for (int j = 0; j < n_chains; j++){
			T_j=T0;
			if (i % step == 0){
				for (int s = 0; s < i/step+1; s++){ 
					T_j*= gamma;
					//std::cout << "T-J gamma" << "\t" << T_j << std::endl;
					
				}
			}
			//std::cout << "current temp" << "\t" << T_j << std::endl;
			for (int z = 0 ; z< display; z++){
				//std::cout << "current temp " << "\t" << T_j << std::endl;
				Model* current_model = models[j][a_indices[j]]; 
				//std::cout << current_model << std::endl;
				Model* new_model =  new_candidate_params_only(current_model,sim_params);
				//std::cout << new_model << std::endl;
				std::vector<double> costs_current (protos.size());
				std::vector<double> costs_new (protos.size());
		
	
				//int dom = domination(current_model,new_model,protos,costs_current,costs_new);
				//std::cout << "costs current" << std::endl;
				//print_vec(costs_current);
				//std::cout << "costs new" << std::endl;
				//print_vec(costs_new);
				//double domq = delta_dom(costs_current,costs_new,(protos.size()));
			
	
				//std::cout << "dom" << dom << std::endl;
				//std::cout << "amount of domination" << domq << std::endl;
				int dom = Math::rng_int(0,3);
				if(dom == 0){ // current dominates new
				//std::cout << "dom" << dom << std::endl;
					//find k: no of points in archive which dominate new pt
					std::vector<int> k_vec (models[j].size());
					/* for(int k = 0; k < models[j].size();k++){ //look over archive
						k_vec[k] = domination(new_model,models[j][k],protos);
					} */
					//std::cout << "k_vec" << std::endl;
					//print_vec(k_vec);
					/* double dom_archive_total = 0;
					auto vec = find_indxs(k_vec,2);
					for(int k = 0; k < vec.size(); k++){
				
						std::vector<double> costs_current_archive_temp (protos.size());
						std::vector<double> costs_new_archive_temp (protos.size());
				
						int dom0 = domination(models[j][vec[k]],new_model,protos,costs_current_archive_temp,costs_new_archive_temp);
						double domq0 = delta_dom(costs_current_archive_temp,costs_new_archive_temp,(protos.size()));
						dom_archive_total = dom_archive_total+domq0;
					} 
			
					double dom_avg = (dom_archive_total+ domq)/(vec.size()+1);
					double prob = 1/(1+exp(dom_avg*T_j)); 
					//double prob = exp(-(dom_avg)/T_j);
					//std::cout << "prob in dom 0" << "\t" << prob << std::endl; 
					double prob = 0.1;
					if( prob < Math::rng_uniform(0,1)){
						//std::cout << "new model accepted that is worse" << std::endl;
						//std::cout << "models before delte" << std::endl;
						//print_vec_vec(models);
						delete (models[j][a_indices[j]]);
						//std::cout << "j" << j << "a_indices[j]" << a_indices[j] << std::endl;
						//std::cout << "models after delte" << std::endl;
						//print_vec_vec(models);
						models[j][a_indices[j]] = new Model(new_model);
						//std::cout << "models after replace new model" << std::endl;
						//print_vec_vec(models);
						//f_val[j][a_indices[j]] = sum_vector(costs_new);
						std::cout << "f_val before" << std::endl;
						print_vec_vec(f_val);
						std::cout << "a-J indx" << std::endl;
							std::cout << a_indices[j] << std::endl;
						f_val[j][a_indices[j]]= Math::rng_uniform(0,10);
						std::cout << "f_val after" << std::endl;
						print_vec_vec(f_val);
						delete (new_model);
					}
					else{
						delete (new_model);
					}
				}
		
		
		
				else if(dom ==1){ //current and new are ND 
					//std::cout << "dom" << dom << std::endl;
					//populate k_vec
					/* std::vector<int> k_vec (models[j].size());
					for(int k = 0; k < models[j].size();k++){ //look over archive
						k_vec[k] = domination(new_model,models[j][k],protos);	
					} */
			
					//std::cout << "k_vec" << std::endl;
					//print_vec(k_vec);
					/* auto vec = find_indxs(k_vec,2);
					auto vec1 = find_indxs(k_vec,1);
					if(vec.size() >= 1){  
					if (0.5 < Math::rng_uniform(0,1)){
								//case 2A //k_vec contains at least one 2
					//std::cout << "K-vec contains at least one 2" << std::endl;
						/*  double dom_archive_total = 0;
						for(int k = 0; k < vec.size(); k++){
				
							std::vector<double> costs_current_archive_temp (protos.size());
							std::vector<double> costs_new_archive_temp (protos.size());
				
							int dom1 = domination(models[j][vec[k]],new_model,protos,costs_current_archive_temp,costs_new_archive_temp);
							double domq1 = delta_dom(costs_current_archive_temp,costs_new_archive_temp,(protos.size()));
							dom_archive_total = dom_archive_total+domq1;
						}
						double dom_avg = (dom_archive_total+ domq)/(vec.size());
						double prob = 1/(1+exp(dom_avg/T_j));   
					//	std::cout << "prob in dom 1" << "\t" << prob << std::endl;
						double prob = 0.1;
						if( prob < Math::rng_uniform(0,1)){
							//std::cout << "new model accepted that is nondom" << std::endl;
							//std::cout << "models before delte" << std::endl;
							//print_vec_vec(models);
							delete (models[j][a_indices[j]]);
							//std::cout << "j" << j << "a_indices[j]" << a_indices[j] << std::endl;
							//std::cout << "models after delte" << std::endl;
							//print_vec_vec(models);
							models[j][a_indices[j]] = new Model(new_model);
							//std::cout << "models after replace new model" << std::endl;
							//print_vec_vec(models);
							//f_val[j][a_indices[j]] = sum_vector(costs_new);
							std::cout << "f_val before" << std::endl;
						print_vec_vec(f_val);
						std::cout << "a-J indx" << std::endl;
							std::cout << a_indices[j] << std::endl;
						f_val[j][a_indices[j]]= Math::rng_uniform(0,10);
						std::cout << "f_val after" << std::endl;
						print_vec_vec(f_val);
							delete (new_model);
						}
						else{
							delete (new_model);
						}
					
					}
					//else if (vec1.size() == models[j].size()){
					else if(Math::rng_uniform(0,1) > 0.5){// case 2B //new point is ND to all of archive
						//std::cout << "K-vec == 1" << std::endl;
						//std::cout << "model push back" << std::endl;
						
						//std::cout << "models after add" << std::endl;
						//print_vec_vec(models);
						if (models[j].size() > SL) {
							single_linkage_cluster(models[j], f_val[j], HL, protos);
						}  
						models[j].push_back(new Model(new_model));
						a_indices[j] = models[j].size()-1;
						//f_val[j].push_back(sum_vector(costs_new));
						f_val[j].push_back(Math::rng_uniform(0,10));
						delete new_model;
						//std::cout << "archive_indx updated" << a_indices[j]  << std::endl;
						//std::cout << "model added" << std::endl;
						 
					}
				
					else{ // case 2C new-pt dominates k >= 1 points in archive 
				
						//std::cout << "k vec has at least one zero" << std::endl;
						//print_vec(k_vec);
						std::vector<Model*> models_temp_new;
						std::vector<double> fvals_temp;
						int count = 0;
						//std::cout << "models before mass delte" << std::endl;
						//print_vec_vec(models);
						for (int k = 0; k < models[j].size(); k++){
							//if(k_vec[k] != 0){ 
							double prob = 0.5;
							if(prob < Math::rng_uniform(0,1)){ 
								models_temp_new.push_back(models[j][k]);
								fvals_temp.push_back(f_val[j][k]);
							}
							else{
								delete models[j][k];
								count++;
							}
						}			
						models[j].swap(models_temp_new);
						f_val[j].swap(fvals_temp);
						//std::cout << "models after mass delte" << std::endl;
						//print_vec_vec(models);
						//std::cout << count << "models deleted" << std::endl;
						//set current point as new point and add to archive
						models[j].push_back(new Model(new_model));
						delete (new_model);
						a_indices[j] = models[j].size()-1;
						//f_val[j].push_back(sum_vector(costs_new));
						f_val[j].push_back(Math::rng_uniform(0,10));
						
						//std::cout << "model added" << std::endl;
						//std::cout << "archive_indx updated" << a_indices[j]  << std::endl
						
					}   
				}   
		
				else{ // Case 3 new dominates current
					//std::cout << "dom" << dom << std::endl;
					//populate k_vec
					/* std::vector<int> k_vec (models[j].size());
					for(int k = 0; k < models[j].size();k++){ //look over archive
							k_vec[k] = domination(new_model,models[j][k],protos);	
						}
						auto vec = find_indxs(k_vec,2);
						auto vec1 = find_indxs(k_vec,1);  
					//if(vec.size() >= 1){ 
					if (Math::rng_uniform(0,1) < 0.5){//case 3A //k_vec contains at least one 2
						//std::vector<double> dom_archive (vec.size());
						//int k = vec.size(); //actually find how many 2s but indexes??
						/* for(int k = 0; k < vec.size(); k++){
				
							std::vector<double> costs_current_archive_temp (protos.size());
							std::vector<double> costs_new_archive_temp (protos.size());
				
							int dom2 = domination(models[j][vec[k]],new_model,protos,costs_current_archive_temp,costs_new_archive_temp);
							double domq2 = delta_dom(costs_current_archive_temp,costs_new_archive_temp,(protos.size()));
							dom_archive[k] = domq2;
						}
						int min_idx = vec_min(dom_archive);
						double prob = 1/(1+exp(-(dom_archive[min_idx]))); 
						//std::cout << "prob in dom 2" << "\t" << prob << std::endl; 
						double prob = 0.1;
						if( prob < Math::rng_uniform(0,1)){
							//std::cout << "new dominates current, pointer updated to min" << std::endl;
							//a_indices[j] = vec[min_idx];
							int min_idx = vec_min(f_val[j]);
							a_indices[j] = min_idx;
							//std::cout << "archive_indx updated" << archive_indx << std::endl;
							delete (new_model);
						}
						else{
							delete (models[j][a_indices[j]]);
							models[j][a_indices[j]] = new Model(new_model);
							//f_val[j][a_indices[j]] = sum_vector(costs_new);
							std::cout << "f_val before" << std::endl;
						print_vec_vec(f_val);
						std::cout << "a-J indx" << std::endl;
							std::cout << a_indices[j] << std::endl;
						f_val[j][a_indices[j]]= Math::rng_uniform(0,10);
						std::cout << "f_val after" << std::endl;
						print_vec_vec(f_val);
							delete (new_model);
						}
					}
			
					//else if(vec1.size() == models[j].size()){	
						else if(Math::rng_uniform(0,1) < 0.5){	
						if (models[j].size() > SL) {
							single_linkage_cluster(models[j], f_val[j], HL, protos);
						}  
						models[j].push_back(new Model(new_model));
						delete (new_model);
						a_indices[j] = models[j].size()-1;
						//f_val[j].push_back(sum_vector(costs_new));
						f_val[j].push_back(Math::rng_uniform(0,5));
						
					}
					else{ // new_point dominates k other points in archive
						int count = 0;
						std::vector<Model*> models_temp_new;
						std::vector<double> fvals_temp;
						for (int k = 0; k < models[j].size(); k++){
							//if(k_vec[k] != 0){ 
							double prob = 0.5;
							if(prob < Math::rng_uniform(0,1)){ 
								models_temp_new.push_back(models[j][k]);
								fvals_temp.push_back(f_val[j][k]);
							}
							else{
								delete models[j][k];
								count++;
							}
						}
						models[j].swap(models_temp_new);
						f_val[j].swap(fvals_temp);
						//std::cout << count << "models deleted" << std::endl;
						models[j].push_back(new Model(new_model));
						delete (new_model);
						a_indices[j] = models[j].size()-1;
						//f_val[j].push_back(sum_vector(costs_new));
						f_val[j].push_back(Math::rng_uniform(0,5));
						//std::cout << "archive_indx updated" << archive_indx << std::endl;
						//std::cout << "model added" << std::endl;
					} 
				} //dominate cases 
			} //z times per display
		}//n_chains
	} // end_i	 
	//std::cout << "come together now" << std::endl;
	//print_vec_vec(models);
	//print_vec_vec(f_val);
	
	int min_chain = min_cost_idx(f_val);
	std::cout << "min_chain" << min_chain << std::endl;
	if(f_val[min_chain][vec_min(f_val[min_chain])] < mini){
		mini = f_val[min_chain][vec_min(f_val[min_chain])];
		std::cout << "new mini" << mini << std::endl;
		best_pareto_so_far.clear();
		best_pareto_so_far = f_val[min_chain];
		for(int k = 0; k < best_paretom_so_far.size(); k++){
			delete(best_paretom_so_far[k]);
		}
		best_paretom_so_far.resize(models[min_chain].size());
		for(int k = 0; k < best_paretom_so_far.size(); k++){
			best_paretom_so_far[k] = new Model(models[min_chain][k]);
		}
				
	}
	print_vec(best_pareto_so_far);
	print_vec(best_paretom_so_far); 
	
	
	//cluster 
	if (best_paretom_so_far.size() > HL) single_linkage_cluster(best_paretom_so_far,best_pareto_so_far,HL, protos);
	
	double min_cost = 0;

for(int i = 0; i < best_paretom_so_far.size(); i++){
		
		/* 	double obj = cost_main(*best_paretom_so_far[i],protos);
			if(obj < min_cost){
				min_cost = obj;
			}
			std::cout << i+1 << ":" << "\t" << obj << std::endl;
			std::ostringstream iss; 
			std::string snapshot_file;
			iss << path << i+1 << ".model";
			snapshot_file = iss.str();

			std::ostringstream mss; 
			std::string modelfit_file;
			mss << path << i+1 << ".txt";
			modelfit_file = mss.str();

			std::ofstream fss(snapshot_file.c_str());
			fss << (*best_paretom_so_far[i]) << std::endl;
			fss.close();
			std::ofstream os(modelfit_file.c_str());

			std::vector<double> modelval1;
			std::vector<double> data1;
			double fn2 = cost_main(*best_paretom_so_far[i],protos); 
			os << "cost" << fn2 << std::endl;
			auto costs = cost_comp((*best_paretom_so_far[i]),protos);
			for(int k = 0; k < costs.size(); k++){
					
					os << costs[k] << "\t";
					
			}
			os << std::endl;
			for (int k = 0; k < protos.size(); k++) {
		
				os << protos[k].name << "\n";
		
				cost((*best_paretom_so_far[i]), protos[k], data1, modelval1);
				os << "Data" << "\t" << "Model" << "\n";
				for (int l = 0; l < protos[k].steps.size(); l++) {

					os << data1[l] << "\t" << modelval1[l] << "\n" ;


				}

				data1.clear();
				modelval1.clear();
				os << std::endl;
			} 
	  
		os.close();  
		delete best_paretom_so_far[i];
	}  
	
	
	for(int i = 0; i < models.size(); i++){
		for (int j = 0; j < models[i].size(); j++){
			delete models[i][j];
		}
	}
	
	return min_cost; 

}  */

/* double SAmosa(std::string path,Model::Model* m, std::vector<ProtocolParameter> protos, SimulationParameters sim_params, std::ofstream& out){


	const int k_max = sim_params.k_max;
	const int step = sim_params.step;
	const int n_chains = sim_params.n_chains;
	const double T0 = sim_params.t0;
	const double gamma = sim_params.gamma;
	const int display = sim_params.display;
	const int snapshot = sim_params.snapshot;
	const std::string snapshotdir = sim_params.snapshotdir;
	int N =  m->G.N;
	//omp_set_dynamic(0);

	//omp_set_num_threads(numThreads);

	int archive_size_initial = 10; //needs to be y*HL
	int archive_iter_initial = 25;
	int SL = 25;
	int HL = 15;
	
	
	
	std::vector<Model*> models (archive_size_initial);
	std::vector<double> f_val (archive_size_initial);
	std::vector<std::vector<double>> model_vals(archive_size_initial);

	//#pragma omp parallel for
	for (int i = 0; i < archive_size_initial; i++){ //initialize a chain of new models

		models[i] =  new_candidate_params_only(m,sim_params);
		//std::cout << models[i] << std::endl;
		f_val[i] = cost_main((*models[i]),protos);
		//std::cout << "model inside of complicated SA" << std::endl;
		//std::cout << "models i" << (models[i]) << std::endl;
		model_vals.push_back(cost_comp((*models[i]),protos));
		//std::cout << "starting cost" << f_val[i] << std::endl;
	}
	
	
	//print_vec_vec(model_vals);
	//implememt bottom sliding to initialize archive
	//only accept models that completely dominate 
	
	
	for(int i = 0; i < archive_size_initial; i++){
			
		for (int j = 0; j < archive_iter_initial; j++){
			Model* model_candidate =  new_candidate_params_only(m,sim_params);
			std::vector<double> costs = cost_comp(*model_candidate,protos);
			double total_cost = sum_vector(costs);
			int better = 0;
			for (int k = 0; k < model_vals[i].size(); k++){
				if(costs[k] <= model_vals[i][k]) better++;
			}
				
			if( better == model_vals[i].size()){ //solution is better for all protos
				delete models[i];
				models[i] = new Model(model_candidate);
				f_val[i] = total_cost;
				delete model_candidate;
			}
			else{
				
				delete model_candidate;
				
			}
		}
	 
		
  
	}		
	
		
		
		//cluster here		
		//bring archive down to HL
		
		
		std::cout << "archive total costs" << std::endl;
		print_vec(f_val);
	
	int archive_indx = Math::rng_int(0,models.size());
	std::cout << "staring archive_indx" << archive_indx << std::endl;
	Model* current_model = models[archive_indx]; 
	
	for (int i = 0; i < k_max; i++){
		
		double T = T0;
		if(i % 100 == 0){
			
			std::cout << "archive size at iteration:" << "\t" << i << "\t" << models.size() << std::endl;
			double min_cost = 10000;
			for(int j = 0; j < models.size(); j++){
		
				double obj = cost_main(*models[j],protos);
				if(obj < min_cost){
					min_cost = obj;
				}
			}
			std::cout << "best cost seen so far" << "\t" << min_cost << std::endl;
		}
		if (i % step){
			for (int s = 0; s < i/step+1; s++){ 
				T*= gamma;
				//std::cout << T_j << std::endl;
			}
		}
	
		//std::cout << "new iteration!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		
		//std::cout << current_model << std::endl;
		Model* new_model =  new_candidate_params_only(current_model,sim_params);
		//std::cout << new_model << std::endl;
		std::vector<double> costs_current (protos.size());
		std::vector<double> costs_new (protos.size());
		
	
		int dom = domination(current_model,new_model,protos,costs_current,costs_new);
		//std::cout << "costs current" << std::endl;
		//print_vec(costs_current);
		//std::cout << "costs new" << std::endl;
		//print_vec(costs_new);
		double domq = delta_dom(costs_current,costs_new,(protos.size()));
			
	
		//std::cout << "dom" << dom << std::endl;
		//std::cout << "amount of domination" << domq << std::endl;
		
		if(dom == 0){ // current dominates new
			//find k: no of points in archive which dominate new pt
			std::vector<int> k_vec;
			for(int j = 0; j < models.size();j++){ //look over archive
				int dom0 = domination(new_model,models[j],protos);
				k_vec.push_back(dom0);
			}
			//print_vec(k_vec);
			double dom_archive_total = 0;
			auto vec = find_indxs(k_vec,2);
			for(int j = 0; j < vec.size(); j++){
				
				std::vector<double> costs_current (protos.size());
				std::vector<double> costs_new (protos.size());
				
				int dom0 = domination(models[vec[j]],new_model,protos,costs_current,costs_new);
				double domq0 = delta_dom(costs_current,costs_new,(protos.size()));
				dom_archive_total = dom_archive_total+domq0;
			}
			
			double dom_avg = (dom_archive_total+ domq)/(vec.size()+1);
			double prob = 1/(1+exp(dom_avg*T));
			//std::cout << "prob" << prob << std::endl;
			if( prob < Math::rng_uniform(0,1)){
				//std::cout << "new model accepted" << std::endl;
				delete (models[archive_indx]);
				models[archive_indx] = new Model(new_model);
				delete (new_model);
			}
			else{
				delete (new_model);
			}
		}
		
		//delete (new_model);
		
		else if(dom ==1){ //current and new are ND 
		
			//populate k_vec
			std::vector<int> k_vec;
			for(int j = 0; j < models.size();j++){ //look over archive
				int dom1 = domination(new_model,models[j],protos);
				k_vec.push_back(dom1);	
			}
			
			//std::cout << "k_vec" << std::endl;
			//print_vec(k_vec);
			auto vec = find_indxs(k_vec,2);
			auto vec1 = find_indxs(k_vec,1);
			if(vec.size() >= 1){ //case 2A //k_vec contains at least one 2
				//std::cout << "K-vec contains at least one 2" << std::endl;
				double dom_archive_total = 0;
				int k = vec.size(); //actually find how many 2s but indexes??
				for(int j = 0; j < vec.size(); j++){
				
					std::vector<double> costs_current (protos.size());
					std::vector<double> costs_new (protos.size());
				
					int dom1 = domination(models[vec[j]],new_model,protos,costs_current,costs_new);
					double domq1 = delta_dom(costs_current,costs_new,(protos.size()));
					dom_archive_total = dom_archive_total+domq1;
				}
				double dom_avg = (dom_archive_total+ domq)/(vec.size());
				double prob = 1/(1+exp(dom_avg*T));
				//std::cout << prob << std::endl;
				if( prob < Math::rng_uniform(0,1)){
					//std::cout << "model accepted" << std::endl;
					delete (models[archive_indx]);
					models[archive_indx] = new Model(new_model);
					delete (new_model);
				}
				else{
					delete (new_model);
					
				}
			} 
			
			else if(vec1.size() == models.size()){ // case 2B //new point is ND to all of archive
				//std::cout << "K-vec == 1" << std::endl;
				models.push_back(new Model(new_model));
				delete (new_model);
				archive_indx = models.size()-1;
				//std::cout << "archive_indx updated" << archive_indx << std::endl;
				//std::cout << "model added" << std::endl;
				//if (models.size() > SL) //cluster(models); down to HL
			}
			else{ // case 2C new-pt dominates k >= 1 points in archive 
				
				//std::cout << "k vec has at least one zero" << std::endl;
				//print_vec(k_vec);
				std::vector<Model*> models_temp_new;
				int count = 0;
				for (int j = 0; j < models.size(); j++){
					if(k_vec[j] != 0){ 
						models_temp_new.push_back(models[j]);
					}
					else{
						delete models[j];
						count++;
					}
					
				}
				models.swap(models_temp_new);
				//std::cout << count << "models deleted" << std::endl;
				//set current point as new point and add to archive
				models.push_back(new Model(new_model));
				delete (new_model);
				archive_indx = models.size()-1;
				//std::cout << "model added" << std::endl;
				//std::cout << "archive_indx updated" << archive_indx << std::endl;
			}
		}  
		
		else{ // Case 3 new dominates current
			
			//populate k_vec
			std::vector<int> k_vec;
			for(int j = 0; j < models.size();j++){ //look over archive
				int dom2 = domination(new_model,models[j],protos);
				k_vec.push_back(dom2);	
			}
			auto vec = find_indxs(k_vec,2);
			auto vec1 = find_indxs(k_vec,1);
			if(vec.size() >= 1){ //case 3A //k_vec contains at least one 2
				std::vector<double> dom_archive (vec.size());
				int k = vec.size(); //actually find how many 2s but indexes??
				for(int j = 0; j < vec.size(); j++){
				
					std::vector<double> costs_current (protos.size());
					std::vector<double> costs_new (protos.size());
				
					int dom2 = domination(models[vec[j]],new_model,protos,costs_current,costs_new);
					double domq2 = delta_dom(costs_current,costs_new,(protos.size()));
					dom_archive[j] = domq2;
				}
				int min_idx = vec_min(dom_archive);
				double prob = 1/(1+exp(-(dom_archive[min_idx])));
				if( prob < Math::rng_uniform(0,1)){
					archive_indx = vec[min_idx];
					//std::cout << "archive_indx updated" << archive_indx << std::endl;
					delete (new_model);
				}
				else{
					delete (models[archive_indx]);
					models[archive_indx] = new Model(new_model);
					delete (new_model);
				}
			}
			
			else if(vec1.size() == models.size()){	
				delete (models[archive_indx]);
				models[archive_indx] = new Model(new_model);
				delete (new_model);
			}
			else{ // new_point dominates k other points in archive
				int count = 0;
				std::vector<Model*> models_temp_new;
				for (int j = 0; j < models.size(); j++){
					if(k_vec[j] != 0){ 
						models_temp_new.push_back(models[j]);
					}
					else{
						delete models[j];
						count++;
					}
				}
				models.swap(models_temp_new);
				//std::cout << count << "models deleted" << std::endl;
				models.push_back(new Model(new_model));
				delete (new_model);
				archive_indx = models.size()-1;
				//std::cout << "archive_indx updated" << archive_indx << std::endl;
				//std::cout << "model added" << std::endl;
			}
		}
		
		
		current_model = models[archive_indx];
		//std::cout << "current archive indx" << archive_indx << std::endl;
		//std::cout << "models.size" << models.size() << std::endl;
	}
		

	
	
	if(models.size() > SL){
		
		
		//cluster down to HL
		
	}
	double min_cost = 10000;
	for(int i = 0; i < models.size(); i++){
		
		double obj = cost_main(*models[i],protos);
		if(obj < min_cost){
			min_cost = obj;
		}
		std::cout << obj << std::endl;
		std::ostringstream iss; 
		std::string snapshot_file;
		iss << path << i << ".model";
		snapshot_file = iss.str();

		std::ostringstream mss; 
		std::string modelfit_file;
		mss << path << i << ".txt";
		modelfit_file = mss.str();

		std::ofstream fss(snapshot_file.c_str());
		fss << models[i] << std::endl;
		fss.close();

	
	
		std::ofstream os(modelfit_file.c_str());

		std::vector<double> modelval1;
		std::vector<double> data1;
		double fn2 = cost_main(*models[i],protos); 
		os << "cost" << fn2 << std::endl;
		auto costs = cost_comp((*models[i]),protos);
		for(int j = 0; j < costs.size(); j++){
					
					os << costs[j] << "\t";
					
		}
		os << std::endl;
			for (int j = 0; j < protos.size(); j++) {
		
				os << protos[j].name << "\n";
		
				cost((*models[i]), protos[j], data1, modelval1);
				os << "Data" << "\t" << "Model" << "\n";
				for (int k = 0; k < protos[j].steps.size(); k++) {

					os << data1[k] << "\t" << modelval1[k] << "\n" ;


				}

				data1.clear();
				modelval1.clear();
				os << std::endl;
			} 
	  
		os.close(); 
		delete models[i];
		
	}

	return min_cost;

} */