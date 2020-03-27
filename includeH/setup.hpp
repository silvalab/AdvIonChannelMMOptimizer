#ifndef SETUP_HPP_
#define SETUP_HPP_

#include <iostream>
#include "MarkovChannel.hpp"
#include <string>
#include <sstream>
#include "math.hpp"
#include "restart.hpp"
#include <vector>
#include <limits.h>
#include <sys/stat.h>
#include <sys/types.h>

class Setup{
	private:
	int N;
	int times;
	int model;
	std::string version;
	public:
	std::vector<ProtocolParameter> protos;
	std::vector<ProtocolParameter> valids;
	std::vector<std::vector<int>> edges;
	std::vector<int> roots;	
	Setup(char*,char*,char*,char*,char*,char*,char*);
	~Setup();
	void load_protocols(char *protolst, std::vector<ProtocolParameter>& p);
	void extract_model_params(std::ifstream& model_list);
	std::string getexepath();
	std::string get_model_pwd_and_create_storage_directories();
	int get_N() {return N;}
	int get_times() {return times;}
	int get_model() {return model;}
	std::string get_version() {return version;}
};

#endif