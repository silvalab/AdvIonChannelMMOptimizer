#include "graph.hpp"



	

	
void Graph::set_topology(int N, int root,std::vector<int> edgesinputvector){
 
	this->N = N;
	this->root = root;
	this->make_pair(edgesinputvector); //populate edge_pairs vector
	for (int i = 0; i < edge_pairs.size(); i++) {

		edges.push_back(Edge(edge_pairs[i].first,edge_pairs[i].second));
		E++;  

	}


 }



void Graph::make_pair(std::vector<int> edgesinputvector){

	if (edgesinputvector.size() % 2 != 0) {
		std::cout << "error: missing edge" << std::endl;
	}

	int twins = edgesinputvector.size()/2;
	for (int i = 1; i <= twins; i++){
		edge_pairs.push_back(std::make_pair(edgesinputvector[(2*i)-2],edgesinputvector[(2*i)-1]));

	}
}



 std::ostream& operator<< (std::ostream& os, Graph& G) {
    for (int i=0; i<G.E; i++) {
      Edge& e = G.edges[i];
      os << e.V1 << "\t" << e.V2 << "\n";
    }
    return os;
  }

