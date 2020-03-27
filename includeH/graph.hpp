#ifndef GRAPH_FUNCTIONS_HPP_
#define GRAPH_FUNCTIONS_HPP_

#include <vector>
#include <string>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <utility>


class Edge {
	public:
    Edge(int v1, int v2) : V1(v1), V2(v2) {};
    Edge(const Edge& e) : V1(e.V1), V2(e.V2) {};
    int V1, V2;
  };

class Graph {
	public:
    Graph() : E(0), N(0) {};
    std::vector<Edge> edges;
	std::vector<std::pair<int,int>> edge_pairs;
    int E, N, root;
	friend std::ostream& operator<< (std::ostream& os, Graph& G);
  
	void make_pair(std::vector<int>);

    void set_topology(int N, int root,std::vector<int>);
  
  
  };







#endif