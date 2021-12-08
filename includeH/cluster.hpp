#ifndef CLUSTER_HPP_
#define CLUSTER_HPP_


#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <list>
#include <limits>
#include "model.hpp"
#include "MarkovChannel.hpp"

static double inf = std::numeric_limits<double>::infinity();
void single_linkage_cluster(std::vector<Model*>& modelsj, std::vector<double>& f_val, int HL, std::vector<ProtocolParameter> protos, int idx);

void single_linkage_cluster(std::vector<Model*>& modelsj, std::vector<double>& f_val, int HL, std::vector<ProtocolParameter> protos);

#endif