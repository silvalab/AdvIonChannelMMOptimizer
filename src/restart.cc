#include "restart.hpp"

Restart_params::Restart_params(std::ifstream& restarttxtfile){
	
	std::string line;
	
	if (restarttxtfile.is_open()){
		
		while(getline(restarttxtfile,line)){
			if(line.compare(0,5,"Time:") == 0){
				times = std::stoi(line.substr(6));
				std::cout << "times" << "\t" << times << std::endl;
			}
			if(line.compare(0,1,"I") == 0){
				size_t colon = line.find(":");
				std::string subline = line.substr(colon+1);
				std::cout << subline << std::endl;
				std::stringstream ss;
				ss << subline;
				std::vector<int> nums;
				std::string temp;
				int found;
				while(!ss.eof()){
					ss >> temp;
					if(std::stringstream(temp)>> found)
						nums.push_back(found);
				}
				iterations = nums[0];
				std::cout << "iterations" << "\t" << iterations << std::endl;
				k_max = nums[1];
				std::cout << "k_max" << "\t" << k_max << std::endl;
			}
			
			if(line.compare(0,2,"So") == 0){
				Sobol_indx = std::stoi(line.substr(12));
				std::cout << "Sobol_indx" << "\t" << Sobol_indx << std::endl;
			}
			if(line.compare(0,6,"params") == 0){
				//std::cout << line.substr(17) << std::endl;
				random_counter_params = std::stoi(line.substr(30));
				std::cout << "random_counter_params" << "\t" << random_counter_params << std::endl;
			} 
			if(line.compare(0,5,"args_") == 0){
				//std::cout << line.substr(17) << std::endl;
				random_counter_params_args = std::stoi(line.substr(33));
				std::cout << "random_counter_params_args" << "\t" << random_counter_params_args << std::endl;
			} 
			if(line.compare(0,3,"uni") == 0){
				//std::cout << line.substr(17) << std::endl;
				random_counter_uni_real = std::stoi(line.substr(28));
				std::cout << "random_counter_params_uni_real" << "\t" << random_counter_uni_real << std::endl;
			} 
			if(line.compare(0,1,"V") == 0){
				//std::cout << "found validate history vector" << std::endl;
				getline(restarttxtfile,line);
				//std::cout << "should be vector line" << std::endl;
				//std::cout << line << std::endl;
				std::istringstream iss(line);
				double x;
				while(iss >> x){ cost_val_history.push_back(x);}
				//disp(cost_val_history);
				
			}
			if(line.compare(0,2,"Co") == 0){
				//std::cout << "found cost history vector" << std::endl;
				getline(restarttxtfile,line);
				//std::cout << "should be vector line" << std::endl;
				//std::cout << line << std::endl;
				std::istringstream iss(line);
				double x;
				while(iss >> x){ cost_tr_history.push_back(x);}
				//disp(cost_tr_history);
				
			}
			if(line.compare(0,2,"PQ") == 0){
				std::cout << "found PQ history vector" << std::endl;
				getline(restarttxtfile,line);
				//std::cout << "should be vector line" << std::endl;
				//std::cout << line << std::endl;
				std::istringstream iss(line);
				double x;
				while(iss >> x){ PQ_history.push_back(x);}
				disp(PQ_history);
				
			}
			
			if(line.compare(0,11,"Times_since") == 0){
				times_since_exchange = std::stoi(line.substr(22));
				std::cout << "times_since_exchange" << "\t" << times_since_exchange << std::endl;
			}
			
			
			if(line.compare(0,12,"Acceptance T") == 0){
				std::cout << "found acceptance times" << std::endl;
				getline(restarttxtfile,line);
				//std::cout << "should be vector line" << std::endl;
				//std::cout << line << std::endl;
				std::istringstream iss(line);
				double x;
				while(iss >> x){ acceptance_times.push_back(x);}
				disp(acceptance_times);
				
			}
			if(line.compare(0,12,"Acceptance R") == 0){
				std::cout << "found acceptance rates" << std::endl;
				getline(restarttxtfile,line);
				//std::cout << "should be vector line" << std::endl;
				//std::cout << line << std::endl;
				std::istringstream iss(line);
				double x;
				while(iss >> x){ acceptance_rates.push_back(x);}
				disp(acceptance_rates);
				
			}
			if(line.compare(0,5,"r vec") == 0){
				std::cout << "found r vec" << std::endl;
				getline(restarttxtfile,line);
				//std::cout << "should be vector line" << std::endl;
				//std::cout << line << std::endl;
				std::istringstream iss(line);
				double x;
				while(iss >> x){ r.push_back(x);}
				disp(r);
				
			}
			
			if(line.compare(0,1,"f") == 0){
			getline(restarttxtfile,line);
			getline(restarttxtfile,line);
			getline(restarttxtfile,line);
			//in fmin_model
			std::vector<int> indicators (3);
			root = std::stoi(line.substr(6));
				std::cout << "root in fmin_model" << "\t" << root << std::endl;
				while(getline(restarttxtfile,line)){
					if(line.compare(0,2,"rs") == 0){
						//std::cout << line << std::endl;
						getline(restarttxtfile,line);
						std::vector<double> rs_temp;
						//std::cout << line << std::endl;
						while(line.compare(1,2,"]") != 0){
							//std::cout << line << std::endl;
							line.erase(0,2); //remove left bracket
							size_t found = line.find("]");
							line.erase(found,1); //remove right bracket 
							//std::cout << line << std::endl;
							std::istringstream iss(line);
							double x;
							while(iss >> x) fmin_rs.push_back(x);
							getline(restarttxtfile,line);
						}
						//model_rs.push_back(rs_temp);
						indicators[0]=1;
						//disp(fmin_rs);
					}
					if(line.compare(0,2,"rk") == 0){
						//std::cout << line << std::endl;
						getline(restarttxtfile,line);
						std::vector<double> rk_temp;
						while(line.compare(1,2,"]") != 0){
							//std::cout << line << std::endl;
							line.erase(0,2); //remove left bracket
							size_t found = line.find("]");
							line.erase(found,1); //remove right bracket 
							//std::cout << line << std::endl;
							std::istringstream iss(line);
							double x;
							while(iss >> x){ fmin_rk.push_back(x);}
							getline(restarttxtfile,line);
						}
						//model_rk.push_back(rk_temp);
						indicators[1] =1;
						//disp(fmin_rk);
					}
					if(line.compare(0,2,"ar") == 0){
						//std::cout << line << std::endl;
						std::vector<double> args_temp;
						getline(restarttxtfile,line); //should be the args number line
						//std::cout << line << std::endl;
						size_t found = line.find("]");
						line.erase(found,1); //remove right bracket 
						//std::cout << line << std::endl;
						std::istringstream iss(line);
						double x;
						while(iss >> x) fmin_args.push_back(x);
						//model_args.push_back(args_temp);
						indicators[2] =1;
						//disp(fmin_args);
					}
					if(indicators[0] == 1 && indicators[1] == 1 && indicators[2] == 1) break;
				}
			
				
				
			}
			if(line.compare(0,1,"C") == 0){
				chains = std::stoi(line.substr(7));
				std::cout << "chains" << "\t" << chains << std::endl;
			}
			
			if(line.compare(0,4,"Root") == 0){
				root = std::stoi(line.substr(6));
				//std::cout << "root" << "\t" << root << std::endl;
				while(getline(restarttxtfile,line)){
					if(line.compare(0,2,"rs") == 0){
						//std::cout << line << std::endl;
						getline(restarttxtfile,line);
						std::vector<double> rs_temp;
						while(line.compare(1,2,"]") != 0){
							//std::cout << line << std::endl;
							line.erase(0,2); //remove left bracket
							size_t found = line.find("]");
							line.erase(found,1); //remove right bracket 
							//std::cout << line << std::endl;
							std::istringstream iss(line);
							double x;
							while(iss >> x) rs_temp.push_back(x);
							getline(restarttxtfile,line);
						}
						model_rs.push_back(rs_temp);
						//disp(model_rs);
					}
					if(line.compare(0,2,"rk") == 0){
						//std::cout << line << std::endl;
						getline(restarttxtfile,line);
						std::vector<double> rk_temp;
						while(line.compare(1,2,"]") != 0){
							//std::cout << line << std::endl;
							line.erase(0,2); //remove left bracket
							size_t found = line.find("]");
							line.erase(found,1); //remove right bracket 
							//std::cout << line << std::endl;
							std::istringstream iss(line);
							double x;
							while(iss >> x) rk_temp.push_back(x);
							getline(restarttxtfile,line);
						}
						model_rk.push_back(rk_temp);
						//disp(model_rk);
					}
					if(line.compare(0,2,"ar") == 0){
						//std::cout << line << std::endl;
						std::vector<double> args_temp;
						getline(restarttxtfile,line); //should be the args number line
						//std::cout << line << std::endl;
						size_t found = line.find("]");
						line.erase(found,1); //remove right bracket 
						//std::cout << line << std::endl;
						std::istringstream iss(line);
						double x;
						while(iss >> x) args_temp.push_back(x);
						model_args.push_back(args_temp);
						//disp(model_args);
					}
				}
			}
			
		}
		
	std::cout << model_rs.size() << std::endl;
	std::cout << model_rk.size() << std::endl;
	std::cout << model_args.size() << std::endl;
		
	}
	
	else {
		
		std::cout << "Warning: Restart file is not open" << std::endl;
	}
	
	
	
	
}