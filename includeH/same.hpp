#ifndef SAME_HPP_
#define SAME_HPP_

#include <iostream>
#include <vector>
#include "model.hpp"
#include "math.hpp"
#include <set>
#include <iterator>
#include <algorithm>


class Graph_Structure{
	
	public:
	double* adj;
	Graph_Structure(int N);
	~Graph_Structure();
	std::vector< std::vector<int> > children;
	//std::vector<int> child_connections;
	std::vector<int> arr;
	//std::vector<int> max_child_seen;
	//std::vector< std::vector<int> > child_count;
	void build_array1(Graph G);
	double*  build_adj(Graph G,double * adj);
	void build_array2(Graph G);
	void build_array3();
	void build_array4();
	int loop;
	int gen;
	int counter;
	void print_array();
	void reset(int N);
	int E;
	int N_states;
	int root;
	std::vector<std::vector<int>> neighbors;
	std::vector<std::vector<int>> back_edges;
	std::vector<std::vector<int>> forward_edges;
	std::vector<std::vector<int>> current_edges;
};

bool include_f(std::vector<int> v, int x);

void print_adj(double * adj, int N);

std::vector<std::vector<int> > printCombination(std::vector<int> v, int r, bool& dup);

void combinationUtil(std::vector<int> v, int data[], int start, int end,
                     int index, int r, std::vector<std::vector<int> >& edgecombo);
					 
int factorial(int N);	

bool compare(Model * m, Model * m1, Graph_Structure * g, Graph_Structure * g1);

void print_vec(std::vector<int> v);

void print_vec_vec(std::vector<std::vector<int>> v);

int find_gen(std::vector<std::vector<int>> children, int k);

bool count_perm(std::vector<int> g, std::vector<int> g1, int N);



#endif