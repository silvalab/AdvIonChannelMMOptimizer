#include "setup.hpp"


Setup::Setup(char* proto_path,char* valid_path,char* N, char* model_list_path,char* model,char* times,char* version){
	
	this->N = atoi(N); 
	this->times = atoi(times);
	this-> model = atoi(model);
	std::ifstream model_list(model_list_path);
	this->version = version;
	extract_model_params(model_list);
	load_protocols(proto_path,protos);
	//std::cout << proto_path << std::endl;
	//std::cout << valid_path << std::endl;
	load_protocols(valid_path,valids);
	std::cout << protos.size() << std::endl;
	std::cout << valids.size() << std::endl;
  }
  
 
Setup::~Setup(){

	
  }
  
  
void Setup::extract_model_params(std::ifstream& model_list){

	std::string line;
	
	if (model_list.is_open()){
		
		int counter = 0;
		
		while(getline(model_list,line)){
			if (line.compare(0,1,"R") == 0){ //start reading in the model
				//std::cout << line << std::endl;
				int root = std::stoi(line.substr(6));
				roots.push_back(root);
				getline(model_list,line);
				std::vector<int> edge;
				while(line.length() != 1){
					//std::cout << line << std::endl;//read in the edges until the white space 
					//std::cout << "line length" << line.length() << std::endl;
					std::stringstream ss(line);
					//in the numbers now
					for(int i = 0; ss >> i; ) {
						edge.push_back(i);
					}
					getline(model_list,line);
				}
				edges.push_back(edge); 
				counter++;
			}	
		}
				
		/* std::cout << "edges" << std::endl;	
		std::cout << ic.size() << std::endl;
		disp(ic);
		std::cout << "roots" << std::endl;
		disp(roots); */
	}
	
	else {
		
		std::cout << "Warning: Model list is not open" << std::endl;
	}
			
}  

void Setup::load_protocols(char *protolst, std::vector<ProtocolParameter>& ps) {
  char datapath[PATH_MAX+1]; //char array for the path of list file, need to include limits for PATH_MAX
  realpath(protolst, datapath); //retrieve directory path from protolst store in datapath
  std::ifstream protofile; //start an import file stream on list file
  protofile.open(protolst);

  std::string p1(datapath); 
  int i = p1.size() - 1;
  while (p1[i] != '/') i--; //remove the .lst name to get ready to append protocol names
  
  std::string pth = p1.substr(0, i+1); //retrieve the home directory sequence common to al protocols
  std::string line;
  std::string proto_file;
  while ( protofile >> line ) { //now cycle through protcols in the list and import
    std::cout << pth+line << std::endl;
    proto_file = pth + line;
	//ProtocolParameter p(proto_file);
   ps.push_back(ProtocolParameter(proto_file));
   /* ProtocolParameter(proto_file);
   //ps.push_back(p); */
	}
	
	
}


std::string Setup::getexepath(){ //https://stackoverflow.com/questions/143174/how-do-i-get-the-directory-that-a-program-is-running-from
		char result[ PATH_MAX ];
		ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
		std::string s( result, (count > 0) ? count : 0 );
		
		int i = s.size() - 2;
		while (s[i] != '/') i--; //remove the .lst name to get ready to append protocol names
  
		std::string path = s.substr(0, i+1); 
		return path;
}



std::string Setup::get_model_pwd_and_create_storage_directories(){

		std::string pwd = getexepath();
		
		std::cout << pwd << std::endl;
		
		
		std::stringstream ss;
		ss << pwd << "State" << N << "/";
		ss << version << "/";
		
		
		std::string path_version = ss.str();
		std::cout << path_version << std::endl;
		struct stat statbuf;
		int isDir = 0;
		if (stat(path_version.c_str(), &statbuf) != -1) {
			if (S_ISDIR(statbuf.st_mode)) {
				isDir = 1;
			}
		}
		if(!isDir){
			int result = mkdir(path_version.c_str(), 0777);
			if(result != 0) std::cout << "problem making version directory" << std::endl;
		}
		
		
		ss << "Model" << model << "/";
		std::string path = ss.str();
		
		int result = mkdir(path.c_str(), 0777); //make model directories
		if(result != 0) std::cout << "problem making model storage directory" << std::endl;
		return path;
	
		
		
}

