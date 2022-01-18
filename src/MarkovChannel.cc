#include "MarkovChannel.hpp"
#include <sstream>
#include <stdlib.h>
#include <limits.h>








std::string proto_parse(std::string line){
	int pos = line.find(":");
	int i = 0;
	std::string parsed = line.substr(0, pos);
	while(std::isspace(line[i])) i++;
	parsed = line.substr(i,pos-i);
	return parsed;
	
}

std::string var_parse(std::string line){
	int pos = line.find(":");
	pos++;
	while(std::isspace(line[pos])) pos++;
	auto parsed = line.substr(pos);
	return parsed;
}


std::string populate_var(std::ifstream& txtfile,std::string line){
	std::getline(txtfile, line);
	std::string parsed = var_parse(line);
	return parsed;
}



StepType string2StepType(std::string line){
												//NONE,PEAK,TAU,TRACE,MIN,TAUD,TAUDD,LATEINA
	std::vector<std::pair<std::string, StepType>> pairs={{ "NONE", NONE}, { "PEAK", PEAK }, 
	{ "TAU", TAU } , {"TRACE", TRACE}, {"MIN", MIN}, {"TAUD", TAUD}, {"TAUDD", TAUDD}, {"LATEINA", LATEINA},{"MOT", MOT},{"MAXPO",MAXPO}}; //look up mapping
	
	for (int i = 0; i < pairs.size(); i++){

		if (pairs[i].first == line){
			
			return pairs[i].second;
		}

	}	
}

SimulationParameters::SimulationParameters(std::ifstream& txtfile){
	
	std::string line;
	if (txtfile.is_open()){
		
		this->init_mu = atof(populate_var(txtfile,line).c_str());
		
		this->init_std = atof(populate_var(txtfile,line).c_str());
		
		this->mut_prob = atof(populate_var(txtfile,line).c_str());
			
		this->update_mu_rates = atof(populate_var(txtfile,line).c_str());
		
		this->update_std_rates = atof(populate_var(txtfile,line).c_str());
		
		this->update_mu_args = atof(populate_var(txtfile,line).c_str());
		
		this->update_std_args = atof(populate_var(txtfile,line).c_str());
		
		this->gamma = atof(populate_var(txtfile,line).c_str());
			
		this->t0 = atof(populate_var(txtfile,line).c_str());
		
		this->rate_min = atof(populate_var(txtfile,line).c_str());
		
		this->rate_max = atof(populate_var(txtfile,line).c_str());
		
		this->arg_min = atof(populate_var(txtfile,line).c_str());
		
		this->arg_max = atof(populate_var(txtfile,line).c_str());
		
		this->k_max = atoi(populate_var(txtfile,line).c_str());
			
		this->step = atoi(populate_var(txtfile,line).c_str());
			
		this->display = atoi(populate_var(txtfile,line).c_str());
			
		this->n_chains = atoi(populate_var(txtfile,line).c_str());
			
		this->restart = atoi(populate_var(txtfile,line).c_str()); 
		
		this->AWS_S3_path = populate_var(txtfile,line);
		//std::cout << "AWS_S3_path" << AWS_S3_path << std::endl;
		this->snapshot = atoi(populate_var(txtfile,line).c_str()); 
				
		
	}
}

ProtocolStep ProtocolParameter::get_step(std::ifstream& prototxt){
	has_dt = 0;
	has_vm = 0;
	has_dt_vm = 0;
	has_extra_args = 0;
	has_stepsize = 0;
	//first place is dt, second vm,third both dt and vm,fourth extra_args,5th stepsize  
	std::vector<double> args_storage;
	double dt, vm,stepsize;
	stepsize = 0.05;
	StepType stype;
	std::string line;
	std::getline(prototxt, line);
	while(line.compare(0,1,"}") != 0){
		std::string parsed = proto_parse(line); //need to get rid of white space for flexibility here
		//std::cout << parsed << std::endl;
		if(parsed == "dt"){ //optional but either or dt or vm needs to be provided
		
			dt = atof((var_parse(line)).c_str());	
			//indicator[0]++;
			has_dt = 1;
		}
		else if(parsed == "vm"){ //optional but either dt or vm needs to be provided

			vm = atof((var_parse(line)).c_str());	
			//std::cout << "vm" << vm << std::endl;
			//indicator[1]++;
			has_vm = 1;
		}
		else if(parsed == "stype"){ //required
		
			std::string parsed_step = var_parse(line);
			stype = string2StepType(parsed_step);
		}
		else if (parsed == "stepsize") { //optinal but default is defined in class
		
			stepsize = atof((var_parse(line)).c_str());
			has_stepsize = 1;
		}
		std::getline(prototxt, line);
	}
	
	if(has_stepsize) return ProtocolStep(dt,stepsize,stype,vm);

	else return ProtocolStep(dt,stype,vm);
	
			
	
}

ProtocolParameter::ProtocolParameter(std::string prototxt_path){
	
	std::ifstream prototxt;
	prototxt.open(prototxt_path);
	std::string line;
	
	int i = prototxt_path.length()-1;
	while (prototxt_path[i] != '/') i--;
	std::string pth = prototxt_path.substr(0, i+1);
	
	//extract base path
	if (prototxt.is_open()){
		this->name = populate_var(prototxt, line);
		this->source = populate_var(prototxt, line);
		this->v0 = atof(populate_var(prototxt, line).c_str());
		this->normalize = std::stoi(populate_var(prototxt, line).c_str());
		this->FLAG = std::stoi(populate_var(prototxt, line).c_str());
		this->weight = std::stoi(populate_var(prototxt, line).c_str());
		this->has_validation_points = std::stoi(populate_var(prototxt, line).c_str());
		
		std::ifstream data;
		std::string data_path = pth+source;
		data.open(data_path);
		//std::cout << pth+source << std::endl;
		
		while(std::getline(data, line)){
			//std::cout << line << std::endl;
			std::istringstream iss(line);
			std::vector<double> v_read;
			double x;
			while(iss >> x) v_read.push_back(x);
			if(v_read.size() == 3){ //SEM included
				this->vars.push_back(v_read[0]);
				this->data.push_back(v_read[1]); 
				this->SE.push_back(v_read[2]);
			}
			else if(v_read.size() == 2){
				this->vars.push_back(v_read[0]);
				this->data.push_back(v_read[1]); 
				this->SE.push_back(0);
			}
			else{
				
				std::cout << "Error: data not read in" << std::endl;
			}
		}
		
		if(has_validation_points){
			
			
			int i = source.length()-1;
			while (source[i] != '.') i--;
			name = source.substr(0, i);
			//std::cout << name << std::endl;
			std::string data_path = pth+name+"_val.dat";
			//std::cout << data_path << std::endl;
			std::ifstream data2;
			data2.open(data_path);
			std::string line;
			this->full_vars = this->vars;
			this->full_data = this->data;
			this->full_SE = this->SE;
			//populate full_data,full_vars,full_SE
			
			while(std::getline(data2, line)){
				
				//std::cout << line << std::endl;
				std::istringstream iss(line);
				std::vector<double> v_read;
				double x;
				while(iss >> x) v_read.push_back(x);
				if(v_read.size() == 3){ //SEM included
					this->full_vars.push_back(v_read[0]);
					this->full_data.push_back(v_read[1]); 
					this->full_SE.push_back(v_read[2]);
					
				}
				else if(v_read.size() == 2){
					this->full_vars.push_back(v_read[0]);
					this->full_data.push_back(v_read[1]); 
					this->full_SE.push_back(0);
					
				}
				else{
					
					std::cout << "Error: full data not read in" << std::endl;
				}
			}
			data2.close();
			
			
			
			
			
			
		}
		
		this->n_traces = this->data.size();
		this->steps.resize(n_traces);
		if(has_validation_points){
			this->n_traces_full = full_data.size();
			this->full_steps.resize(n_traces_full);
		}
		
		std::vector<ProtocolStep> repeat_steps;
		
		while(std::getline(prototxt, line)){
			if (line.compare(0,1,"R") == 0){ //in a repeat step  
					//std::cout << "R" << std::endl;
					repeat_steps.push_back(get_step(prototxt));
			}
		}
		
		
		if(repeat_steps.size() != 0){
			if(has_validation_points){
				for(int i = 0; i < n_traces_full; i++){
					for (int pulse = 0; pulse < 99; pulse++){
						for(int j = 0; j < repeat_steps.size(); j++){
							full_steps[i].push_back(repeat_steps[j]);
						}
					}
				}
			}
			for(int i = 0; i < n_traces; i++){
				for (int pulse = 0; pulse < 99; pulse++){
					for(int j = 0; j < repeat_steps.size(); j++){
							
							steps[i].push_back(repeat_steps[j]);
						}
				}
			}
		}
		prototxt.clear();
		prototxt.seekg (0, std::ios::beg);
		while(std::getline(prototxt, line)){
			if (line.compare(0,2,"st") == 0){ // we are in a step, either dt or vm is not given
				//std::cout << "in regular step" << std::endl;
				has_dt = 0;
				has_vm = 0;
				has_dt_vm = 0;
				has_extra_args = 0;
				//first place is dt, second vm,third both dt and vm,fourth extra_args,5th stepsize  
				std::vector<double> args_storage;
				double dt = 0;
				double vm = -500;
				double stepsize = 0.05;
				StepType stype;
				std::getline(prototxt, line);
				while(line.compare(0,1,"}") != 0){
					std::string parsed = proto_parse(line); //need to get rid of white space for flexibility here
					//std::cout << parsed << std::endl;
					if(parsed == "dt"){ //optional
					
						dt = atof((var_parse(line)).c_str());	
						//indicator[0]++;
						has_dt = 1;
					}
					else if(parsed == "vm"){ //optional

						vm = atof((var_parse(line)).c_str());	
						//std::cout << "vm" << vm << std::endl;
						//indicator[1]++;
						has_vm = 1;
					}
					else if(parsed == "stype"){ //required
					
						std::string parsed_step = var_parse(line);
						stype = string2StepType(parsed_step);
					}
					else if (parsed == "stepsize") { //required
					
						stepsize = atof((var_parse(line)).c_str());
					}
					else{ //extra_args optional
					
						args_storage.push_back(atof((var_parse(line)).c_str()));
						//indicator[3]++;
						has_extra_args = 1;
					}
					std::getline(prototxt, line);
				}
				
		
		
		if( stype == TRACE ){
					this->sweeps = 1;
					this->steps.resize(sweeps);
				if(has_validation_points){
					std::vector<double> trace_val_times;
					this->sweeps_full = 1;
					this->full_steps.resize(sweeps_full);
					
					for (int i = n_traces; i < n_traces_full; i++){
						trace_val_times.push_back(this->full_vars[i]);
					}
					//std::cout << "trace_val_times" << std::endl;
					//disp(trace_val_times);
					
					for(int i = 0; i < trace_val_times.size(); i++){
						
						this->trace_val_indxs.push_back(find(this->vars, trace_val_times[i]));
							
					} 
					
					//disp(trace_val_indxs);
				}
				
				
		}
		else{
			
			this->sweeps = this->data.size();
			this->sweeps_full = full_data.size();
		}
		
				
				if (has_dt && has_vm) has_dt_vm = 1;
				 if(has_validation_points){
					
						for (int i= 0; i < sweeps_full; i++){
							
							if(has_extra_args){ //extra arg present 
							std::cout << "extra args" << std::endl;
								if(has_dt_vm) //both dt and vm populated
									full_steps[i].push_back(ProtocolStep(dt,stepsize,stype,vm,args_storage));
								else{ //either dt or vm is missing
									if(has_dt){ //dt is present, need to import vm from vars
										full_steps[i].push_back(ProtocolStep(dt,stepsize,stype,this->full_vars[i],args_storage));
									}
									else{ //dt is missing so input from vars
										full_steps[i].push_back(ProtocolStep(this->full_vars[i],stepsize,stype,vm,args_storage));
										
									}
								}
							}
							else{ // no extra args present
								if(has_dt_vm) {//both dt and vm populated
									
									full_steps[i].push_back(ProtocolStep(dt,stepsize,stype,vm));
								}
								else{ //either dt or vm is missing
									if(has_dt){ //dt is present, need to import vm from vars
										full_steps[i].push_back(ProtocolStep(dt, stepsize,stype,this->full_vars[i]));
									}
									else{ //dt is missing so input from vars
										full_steps[i].push_back(ProtocolStep(this->full_vars[i], stepsize,stype,vm));
									}
								}
							}
					}
					
				} 
				
				for (int i= 0; i < sweeps; i++){
					if(has_extra_args){ //extra arg present 
						if(has_dt_vm) //both dt and vm populated
							steps[i].push_back(ProtocolStep(dt,stepsize,stype,vm,args_storage));
						else{ //either dt or vm is missing
							if(has_dt){ //dt is present, need to import vm from vars
								steps[i].push_back(ProtocolStep(dt,stepsize,stype,this->vars[i],args_storage));
							}
							else{ //dt is missing so input from vars
								steps[i].push_back(ProtocolStep(this->vars[i],stepsize,stype,vm,args_storage));
								
							}
						}
					}
					else{ // no extra args present
						if(has_dt_vm) {//both dt and vm populated
						
							steps[i].push_back(ProtocolStep(dt,stepsize,stype,vm));
						}
						else{ //either dt or vm is missing
							if(has_dt){ //dt is present, need to import vm from vars
								steps[i].push_back(ProtocolStep(dt, stepsize,stype,this->vars[i]));
							}
							else{ //dt is missing so input from vars
								steps[i].push_back(ProtocolStep(this->vars[i], stepsize,stype,vm));
							}
						}
					}
				} 
			} 
		}
	}
}


void SimulationParameters::print_sim_params(std::ofstream& out){
	
	out << "k_max" << "\t" << k_max << std::endl;
	out << "step" << "\t" << step << std::endl;
	out << "n_chains" << "\t" << n_chains << std::endl;
	out << "T0" << "\t" << t0 << std::endl;
	out << "gamma" << "\t" << gamma << std::endl;
	out << "mut_prob" << "\t" << mut_prob << std::endl;
	out << "update_mu_rates" << "\t" << update_mu_rates << std::endl;
	out << "update_std_rates" << "\t" << update_std_rates << std::endl;
	out << "update_mu_args" << "\t" << update_mu_args << std::endl;
	out << "init_std_args" << "\t" << update_std_args << std::endl;
	out << "display" << "\t" << display << std::endl;
	out << "snapshot" << "\t" << snapshot << std::endl;
	out << "rate_min"  << "\t" << rate_min << std::endl;
	out << "rate_max"  << "\t" << rate_max << std::endl;
	out << "arg_min"  << "\t" << arg_min << std::endl;
	out << "arg_max"  << "\t" << arg_max << std::endl;
}


std::ostream& operator<< (std::ostream &os,ProtocolStep protostep){
  os << "dt:\t" << protostep.dt << std::endl;
  os << "vm:\t" << protostep.vm << std::endl;
  os << "stepsize:\t" << protostep.stepsize << std::endl;
  
  if (protostep.extra_args.size() > 0){
	  
	  os << "extra_args:\t" << protostep.extra_args[0] << std::endl;
	  os << "extra_args:\t" << protostep.extra_args[1] << std::endl;
  }
  //{NONE,PEAK,TAU,TRACE,MIN,TAUD,TAUDD,LATEINA};

  //os << "StepType:\t" << protostep.stype << std::endl;
  switch ( protostep.stype ) {
    case 0: os << "NONE " << std::endl; break;
    case 1: os << "PEAK " << std::endl; break;
    case 2: os << "TAU  " << std::endl; break;
    case 3: os << "TRACE" << std::endl; break;
	case 4: os << "MIN" << std::endl; break;
	case 5: os << "TAUD" << std::endl; break;
	case 6: os << "TAUDD" << std::endl; break;
	case 7: os << "LATEINA" << std::endl; break;
	case 8: os << "MOT" << std::endl; break;
	case 9: os << "MAXPO" << std::endl; break;
  }

  return os;
}

std::ostream& operator<<(std::ostream &os,ProtocolParameter protoparam){
	os << "name:\t" << protoparam.name << std::endl;
	os << "source:\t" << protoparam.source << std::endl;
	os << "v0:\t" << protoparam.v0 << std::endl;
	os << "temperature:\t" << protoparam.temperature << std::endl;
	os << "normalize:\t" << protoparam.normalize << std::endl;
	os << "weight:\t" << protoparam.weight << std::endl;
	os << "has val:\t" << protoparam.has_validation_points << std::endl;
	os << "FLAG:\t" << protoparam.FLAG << std::endl;
	os << "C:\t" << protoparam.CONDUCTANCE << std::endl;
	
	for (int i = 0; i < protoparam.steps.size(); i++){

		for (int j = 0; j < protoparam.steps[i].size(); j++){
			
			os << protoparam.steps[i][j] << std::endl;
		}

	}	
	
	
	if(protoparam.has_validation_points){
		for (int i = 0; i < protoparam.full_steps.size(); i++){

			for (int j = 0; j < protoparam.full_steps[i].size(); j++){
				
				os << protoparam.full_steps[i][j] << std::endl;
			}

		}	
		
		
		
		for(int i = 0; i < protoparam.full_data.size(); i++){
			
			std::cout << protoparam.full_vars[i] << "\t" << protoparam.full_data[i] << "\t" << protoparam.full_SE[i] << std::endl;
		}
		
	}
	
	
	return os;
}


