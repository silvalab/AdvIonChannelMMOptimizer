#include <fcntl.h>
#include "MarkovChannel.hpp"
#include "graph.hpp"
#include "model.hpp"
#include <stdlib.h>
#include <string>
#include <sstream>
#include "time.h"
#include <algorithm>
#include "math.hpp"
#include "SA.hpp"
#include <sys/stat.h>
#include <sys/types.h>
#include <omp.h>
#include "restart.hpp"
#include <limits.h>
#include <unistd.h>
#include "setup.hpp"


//******COMMAND LINE OPTIONS KEY*********
// argv[1] solver.txt contains the optimization specifics such as n_chain, iteration #s, parameter ranges, etc
// argv[2] path to protocol data list
// argv[3] path to validation data list
// argv[4] # of states of model to optimize
// argv[5] path to parsed model topology list created with Python
// argv[6] model index with # of states to optimize	
// argv[7] path to Sobol direction numbers (please see includeH/modified_sobol.hpp)	
// argv[8] number of Sobol starts requested
// argv[9] optimization name as string	
// OPTIONAL argv[10] path to iter_state.txt produced from the previously stopped optimization needed to restart

//***********Please see included file.sh for examples of the command line arguments when calling the program*******************

	
int main(int argc, char *argv[]) {
	
	 if (argc == 11){ // restart the simulation if not finished

		std::ifstream simtxt(argv[1]); //simulation parameters read in
		SimulationParameters sim_params(simtxt); //create simulation parameters object based on input solver.txt
		
		
		Setup init(argv[2],argv[3],argv[4],argv[5],argv[6],argv[8],argv[9]); //read in simulation options from command line and store in Setup object
		
		std::ifstream statetxt(argv[10]); //load in restart params from iter_state.txt produced the last optimization display point
		Restart_params restart_params(statetxt);
		
		Math math_params(restart_params.get_random_counter_uni_real(),restart_params.get_Sobol_indx(),(((init.get_N()-1)*2)+(init.edges[init.get_model()-1].size())+2),argv[7]);
		std::cout << "random counter_uni_real" << math_params.get_random_counter_uni_real() << std::endl;
		std::cout << "Sobol indx" << math_params.get_sobol_indx() << std::endl;
		
		std::string path = init.get_model_pwd_and_create_storage_directories(); //get storage path and make sure directories are there
		
		
		
		std::ostringstream ss; //read in the previous Model cost progress to get ready to write again
		ss << "Model" << init.get_model() << ".txt";
		std::string path_model = path+ss.str();
		std::ofstream out;
		out.open(path_model,std::ios_base::app);
		std::cout << "Examining Model #" << "\t" << init.get_model() << std::endl;
		ss.str("");
		
		
		ss << "Time" << restart_params.get_times() << "/";
		std::string path_restart = path+ss.str();
		ss.str("");
		
		if(restart_params.get_iterations() != (restart_params.get_k_max())){ //if where we left off is not fully finished with a start, finish that start
			SA sa(path_restart,restart_params.get_times(),restart_params); //create structure to finish the start 
			sa.anneal_restart(math_params,init,sim_params,restart_params,out);
		}
		
		
		Model working_model(init.get_model(),init.get_N(),init.edges[init.get_model()-1],init.roots[init.get_model()-1],sim_params);
		for( int k = restart_params.get_times(); k < init.get_times(); k++){ //finish the rest of times(starts) like before
			out << "Time:" << "\t" << k+1 << std::endl;
			ss << "Time" << k+1 << "/";
			std::string path_time = path+ss.str();
			int result = mkdir(path_time.c_str(), 0777); //make time directories under the path of the simulation
			ss.str("");
			if(result != 0) std::cout << "possible problem making time directory or exists already" << std::endl;
			working_model.sobol(math_params,sim_params); //print sobol transformed rate rate parameters for the current topology
			std::cout << working_model << std::endl;
			SA sa(path_time,k+1); //Setup simulated annealing structure
			sa.anneal(math_params,working_model,init,sim_params,out); //anneal and find fminmodel while printing and saving optimization history, saving to Amazon S3 happens here
	
		}
		
	} 
	else{ // first time optimization ran		
		
		
        std::ifstream simtxt(argv[1]); //simulation parameters read in
		SimulationParameters sim_params(simtxt); //create simulation parameters object based on input solver.txt
		
		
		Setup init(argv[2],argv[3],argv[4],argv[5],argv[6],argv[8],argv[9]); //read in simulation options from command line and store in Setup object
		
		
		Math math_params((((init.get_N()-1)*2)+(init.edges[init.get_model()-1].size())+2),argv[7]); //set up math object that handles random number generation
		
		std::stringstream ss;
		std::string path = init.get_model_pwd_and_create_storage_directories(); //create version and model optimization storage directories;
		ss << "Model" << init.get_model() << ".txt";
		std::string path_model = path+ss.str(); //create std::string for output model file
		ss.str("");
		std::ofstream out(path_model);
		out << "Examining Model #" << "\t" << init.get_model() << std::endl;
		std::cout << "Examining Model #" << "\t" << init.get_model() << std::endl;
		sim_params.print_sim_params(out);
		Model working_model(init.get_model(),init.get_N(),init.edges[init.get_model()-1],init.roots[init.get_model()-1],sim_params); //load in model topology to optimize
		
		
			  for( int k = 0; k < init.get_times(); k++){
				out << "Time:" << "\t" << k+1 << std::endl;
				ss << "Time" << k+1 << "/";
				
				std::string path_time = path+ss.str();
				int result = mkdir(path_time.c_str(), 0777); //make time directories under the path of the simulation
				ss.str("");
				if(result != 0) std::cout << "possible problem making time directory or exists already" << std::endl;
				working_model.sobol(math_params,sim_params); //print sobol transformed rate rate parameters for the current topology
				std::cout << working_model << std::endl;
				//model_penalty(working_model);
				
				SA sa(path_time,k+1); //Setup simulated annealing structure
				sa.anneal(math_params,working_model,init,sim_params,out); //anneal and find fminmodel while printing and saving optimization history,
				
			 } 
	}
	
	return 0;
}
