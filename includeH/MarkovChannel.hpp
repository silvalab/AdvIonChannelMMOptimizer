#ifndef MarkovChannel_HPP_
#define MarkovChannel_HPP_

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "helper.hpp"

class SAParameter;
class ProtocolStep;
class ProtocolParameter;


class MarkovChannel{
	public:
	MarkovChannel(std::ifstream& txtfile);
	std::vector<ProtocolParameter> protos;

};

enum StepType {NONE,PEAK,TAU,TRACE,MIN,TAUD,TAUDD,LATEINA,MOT,MAXPO};
class SimulationParameters {
	
	public:
	SimulationParameters(std::ifstream& txtfile);
	SimulationParameters() {};	
	void print_sim_params(std::ofstream& out);
	double init_mu;
	double init_std;
	double mut_prob;
	double update_std_rates;
	double update_mu_rates;
	double update_std_args;
	double update_mu_args;
	double gamma=0.9;
	double t0;
	double rate_min;
	double rate_max;
	double arg_min;
	double arg_max;
	int k_max;
	int step;
	int display;
	int n_chains;
	int restart;
	std::string AWS_S3_path;
	int snapshot;
};



class ProtocolParameter {
  private:
	int  CONDUCTANCE = 1;
	int  FLUORIMETRY = 0;
	std::string name;
	std::string source;
	double temperature =293;
  public:
  ProtocolParameter(std::string prototxt);
  int has_dt = 0;
  int has_vm = 0;
  int has_dt_vm = 0;
  int has_stepsize = 0;
  int has_extra_args = 0;
  int n_traces = 0;
  int n_traces_full = 0;
  int sweeps = 0;
  int sweeps_full = 0;
  std::vector<double> vars;
  std::vector<double> data;
  std::vector<double> SE;
  int has_validation_points = 0;
  std::vector<double> full_data;
  std::vector<double> full_vars;
  std::vector<double> full_SE;
  std::vector<int> trace_val_indxs;
  double v0 = -500;
 
  int normalize = 0;
  double weight = 0;
  int FLAG = 0; 
  friend std::ostream& operator<<(std::ostream &os,ProtocolParameter protoparam); 
  std::vector<std::vector<ProtocolStep>> steps;
  std::vector<std::vector<ProtocolStep>> full_steps;
  std::string get_name() {return name;}
  ProtocolStep get_step(std::ifstream& prototxt);
  
};


class ProtocolStep {
 private:
	
	double dt;
	double vm;
	StepType stype;
	double stepsize = 0.05;
	std::vector<double> extra_args;
 public:
  
  double get_dt() {return dt;}
  double get_vm() {return vm;}
  StepType get_stype() {return stype;}
  double get_stepsize() {return stepsize;}
  std::vector<double> get_extra_args() {return extra_args;}
  friend std::ostream& operator<< (std::ostream &os,ProtocolStep protostep);
  ProtocolStep() {}
  ProtocolStep(double dt, StepType stype, double vm): dt(dt), stype(stype), vm(vm){}
  ProtocolStep(double dt, StepType stype, double vm, std::vector<double> extra_args): dt(dt),stype(stype), vm(vm), extra_args(extra_args){}
  ProtocolStep(double dt, double stepsize, StepType stype, double vm): dt(dt), stepsize(stepsize), stype(stype), vm(vm){}
  ProtocolStep(double dt, double stepsize, StepType stype, double vm, std::vector<double> extra_args): dt(dt), stepsize(stepsize),stype(stype), vm(vm), extra_args(extra_args){}
  //Solver<T>(double tinit, T yinit, double tstep): told(tinit), yold(yinit), h(tstep) {}
};

#endif