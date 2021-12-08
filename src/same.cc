#include "same.hpp"

using namespace std;

void print_vec(std::vector<int> v){
	
	
	for (int i = 0; i < v.size(); i++){
		
		
		std::cout << v[i] << "\t";
		
	}
	
	std::cout << std::endl;
}

void print_vec_vec(std::vector<std::vector<int>> v){
	
	
	for (int i = 0; i < v.size(); i++){
		std::cout << i << "|";
		for (int j = 0; j < v[i].size(); j++){
		
		std::cout << v[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	
	
}



bool count_perm(vector<int> g, vector<int> g1, int N){ // g1 and g are aleady the same size
	
	vector<int> counts (N,0); 
	vector<int> counts1 (N,0);
	
	for (int i = 0; i < g.size(); i++){
		counts[g[i]]++;
		counts1[g1[i]]++;		
	}
	
	if (counts == counts1){
		return true;
	}
	else {
		return false;
	}
}

bool include_f(vector<int> v, int x){
	
	for (int i = 0; i < v.size(); i++){
		
		if(v[i] == x){
			return true;
		}
	}
	return false;
}

bool include_f_set(vector<int> v){
	
	for (int i = 0; i < v.size(); i++){
		
		if (v[i] == 1){
			
			return true;
			
		}
	}
	
	return false;
	
}

int factorial(int N){
	
	if (N ==1 || N==0) {

	return 1;

	}
	
	else {
	return factorial(N-1)*N;
	}
	
}

void combinationUtil(vector<int> v, int data[], int start, int end,
                     int index, int r, vector<vector<int> >& edgecombo)
{
    // Current combination is ready to be printed, print it
    if (index == r)
    {	vector<int> temp;
        for (int j=0; j<r; j++){
			//printf("%d ", data[j]);
			
			temp.push_back(data[j]);	
			
		}
		//printf("\n");
        edgecombo.push_back(temp); 
        return;
    }
 
	
	///printf("%d ", data[j]);
        //printf("\n");
    // replace index with all possible elements. The condition
    // "end-i+1 >= r-index" makes sure that including one element
    // at index will make a combination with remaining elements
    // at remaining positions
	//&& end-i+1 >= r-index
    for (int i=start; i < end; i++)
    {
        data[index] = v[i];
        combinationUtil(v, data, i+1, end, index+1, r, edgecombo);
		
		
		if (i < (end-1)){
			// Remove duplicates
			if (v[i] == v[i+1])
				i++;
		}
    }
}

vector<vector<int> > printCombination(vector<int> v, int r, bool& dup)
{	

	//http://www.geeksforgeeks.org/print-all-possible-combinations-of-r-elements-in-a-given-array-of-size-n/
	/* arr[]  ---> Input Array
   data[] ---> Temporary array to store current combination
   start & end ---> Staring and Ending indexes in arr[]
   index  ---> Current index in data[]
   r ---> Size of a combination to be printed */
   
	int n = v.size();
	sort(v.begin(),v.end());
	
	for (int i = 0; i < v.size()-1; i++){
		
		if (v[i] == v[i+1]){
		
			dup = 1;
		}
		
		
	}
	
	
	//int num_combos = factorial(n)/((2)*factorial(n-2));
	
	//std::cout << "num_combos" << num_combos << std::endl;
	vector<vector<int> > edgecombos;
    // A temporary array to store all combination one by one
    int data[r];
 
    // Print all combination using temprary array 'data[]'
    combinationUtil(v, data, 0, n, 0, r, edgecombos);
	
	
	return edgecombos;
}

void print_adj(double * adj, int N){
	
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
		
		std::cout << adj[i*N + j] << "\t";
		
		
		}
		std::cout << std::endl;
	}
}



int find_gen(std::vector<std::vector<int>> children, int k){
	
	for (int i = 1; i < children.size(); i++){
		
		bool inc = include_f(children[i],k);
		if (inc){
			return i;
		}
		
		
	}
	
	return 0;
}

int contain_vector(std::vector<int> v, std::vector<std::vector<int>> col){
		
	for (int i = 0; i < col.size(); i++){
		
		if (v == col[i]){
			//std::cout << "vector should be equal to collection" << std::endl; 
			//print_vec(v);
			//print_vec(col[i]);
			return i;
		}		
	}
		return -1;
}


bool vec_perm(std::vector<std::vector<int>> Edges,std::vector<std::vector<int>> Edges1){
	
	for (int i = 0; i < Edges.size(); i++){
		
		int indx_check = contain_vector(Edges[i], Edges1);
		//print_vec(Edges[i]);
		//print_vec_vec(Edges1);
		//cout << "indx_check" << indx_check << std::endl;
			if (indx_check < 0){
			
				return false; //can't find find an exact vector match so it is not included
			
			}
		
	}
	
	
		return true; //find indices for all the vectors the so the collection is a permutation
	
}


bool col_perm(std::vector<std::vector<int>> Edges,std::vector<std::vector<int>> Edges1, int gen, std::vector<std::vector<int>> children, std::vector<std::vector<int>> children1){
	
	// tests whether the Edges vector collections are permutations of each other
	std::vector<std::vector<int>> Edgestorage;
	std::vector<std::vector<int>> Edge1storage;
	for (int k = 1; k <= gen; k++){ 
		Edgestorage.resize(children[k].size(), vector<int> (3));
		Edge1storage.resize(children1[k].size(), vector<int> (3));
			for (int i = 0; i < children[k].size(); i++){
			
				for (int j =0; j < 3; j++){
				
				
					//std::cout << Edges[k-1][i*3+j] << "\t";
					//std::cout << std::endl;
					//std::cout << Edges1[k-1][i*3+j] << "\t";
					Edgestorage[i][j] = Edges[k-1][i*3+j];
					Edge1storage[i][j] = Edges1[k-1][i*3+j];
					//std::cout << "testing storage" << std::endl;
					//std::cout << Edgestorage[i][j] << std::endl;
					//std::cout << Edge1storage[i][j] << std::endl;
				}
			
				//std::cout << std::endl;
			
			}
			
			
			/* std::cout << "testing storage" << std::endl;
			for (int i = 0; i < children[k].size(); i++){
			
				for (int j =0; j < 3; j++){
					
					std::cout << Edgestorage[i][j] << "\t";
					//std::cout << Edge1storage[i][j] << std::endl;
				}
			
				std::cout << std::endl;
			
			}
			
			
			std::cout << "testing storage1" << std::endl;
			for (int i = 0; i < children[k].size(); i++){
			
				for (int j =0; j < 3; j++){
					
					//std::cout << Edgestorage[i][j] << std::endl;
					std::cout << Edge1storage[i][j] << "\t";
				}
			
				std::cout << std::endl;
			
			} */
			
			
			//call function here to determine whether collection of vectors are same
				if(!vec_perm(Edgestorage,Edge1storage)){
					return false;
				}
			//std::cout << "Edgestorage.size()" << Edgestorage.size()  << std::endl;
		}
	
		
	
		
	
	return true;
	
}





Graph_Structure::Graph_Structure (int N){
	
	adj = new double[N*N];
	children.resize(N);
	back_edges.resize(N);
	current_edges.resize(N);
	forward_edges.resize(N);
	neighbors.resize(N);
	loop = 0;
	gen = 0;
	E = 0;
	N_states = N;
	counter = 0;
}

Graph_Structure:: ~Graph_Structure(){
	
	
	delete [] adj;
	
	
	
}

void Graph_Structure:: reset(int N) {
	arr.clear();
	children.resize(N);
	for (int i = 0; i < N; i++){
		
		children[i].clear();
		
	}
	//back_edges.resize(N);
	for (int i = 0; i < N; i++){
		
		back_edges[i].clear();
		
	}
	//current_edges.resize(N);
	for (int i = 0; i < N; i++){
		
		current_edges[i].clear();
		
	}
	//forward_edges.resize(N);
	for (int i = 0; i < N; i++){
		
		forward_edges[i].clear();
		
	}
	for (int i = 0; i < N; i++){
		
		neighbors[i].clear();
		
	}
	loop = 0;
	gen = 0;
	E = 0;
	N_states = N;
	counter = 0;
}

void Graph_Structure:: print_array(){
		

	 std::cout << arr[0] << "\t";
	 std::cout << "|" << std::endl;
	
	vector<int> stops;
	int idx = 0;
	for (int i = 1; i <= gen; i++){
		
		idx = idx+children[i].size();
		stops.push_back(idx);
		
		/* if (child_connections[i] != 0){
		
			idx = idx+(child_connections[i]*2);
			stops.push_back(idx);
		
		} */
	}
	
	/* for (int i = 0; i < stops.size(); i++){
		
		
		
		cout << stops[i] << endl;
		
		
	} */
	
	int iter_step = 0;
	for (int i = 1; i < arr.size(); i++){
		
		
		
		std::cout << arr[i] << "\t";
	
		if (i == stops[iter_step]){
			std::cout << "|" << std::endl;
			iter_step++;
		}
		
		
		
	}
	
	
	cout << endl;
	cout << "Edges" << E << endl;
	cout << "gen" << gen << endl;
	
//cout << "N_states" << N_states << endl;
for (int i = 0; i <N_states; i++){
	 
	cout << back_edges[i].size() << endl;
	cout << forward_edges[i].size() << endl;
	cout << current_edges[i].size() << endl;
}  
}


std::vector<std::pair <int,int>> edges_to_loop(int N){
	
std::vector<std::pair <int,int>> edges;

	
	for (int i = 0; i < N-1; i++){

		std::pair <int,int> foo;
		edges.push_back(foo = std::make_pair(i,i+1));


	}
	
	
	return edges;	
	
}
	


//need to make separate functions


double*  Graph_Structure::build_adj(Graph G,double * adj){
int N = G.N;	
int tot_num_edges = ((N*N)-N)-2;
	
	
for (int i = 0; i < N*N; i++){
	
	adj[i] = 0;
	
}

for (int i = 0; i < G.E; i++){
	
	
	
	int e1 = G.edges[i].V1;
	int e2 = G.edges[i].V2;
	
	adj[e1*N + e2] = 1;
	adj[e2*N + e1] = 1;
	
	
}

//print_adj(adj,N);



return adj;	

}

void Graph_Structure::build_array1(Graph G){
	
	//std::cout << "build array called" << std::endl;
	//std::cout << (*m) << std::endl;
	gen = 0;
	N_states= G.N;
	E = G.E;
	arr.push_back(G.root);	
	children[gen].push_back(G.root);
	adj = build_adj(G,adj);
	
}

void Graph_Structure::build_array2(Graph G){
	//print_adj(adj,N);
	gen++;
	counter = 1;	
	//std::cout << "children[1].size()" << children[1].size() << std::endl;
	//children[1].clear();
	//std::cout << "children[1].size()" << children[1].size() << std::endl;
	for (int i = (N_states*G.root); i < ((G.root*N_states)+N_states); i++){ //populate the 1st generation vector 
		if (adj[i] == 1){
			children[gen].push_back(i%N_states);
			//std::cout << "i" << i << std::endl;
			arr.push_back(i);
			counter++; //nodes classified
		}

	}
	
}	


void Graph_Structure::build_array3(){	
vector<int> include_f_vector;
while( counter < N_states){ //populate the rest of the generations of children
	gen++;
	//std::cout << "gen" << gen << "counter" << counter << std::endl;
	//include_f_vector.resize(gen);
	for (int i = 0; i < children[gen-1].size(); i++){
			int child_counter = 0;
			for (int j = 0; j < N_states; j++){
					include_f_vector.resize(gen+1);
					for (int k = 0; k <= gen; k++){
						include_f_vector[k] = include_f(children[k],j);
						}
					if (!include_f_set(include_f_vector)){
						int e1 = children[gen-1][i];
						if ((adj[(e1*N_states)+ j] == 1)){
							children[gen].push_back(j);
							arr.push_back(j);
							counter++; //nodes classified
							child_counter++;
						}
					}
			}
	}
}
}



void Graph_Structure::build_array4(){
//print_adj(adj,N_states);
int node;
for (int i = 0; i <= gen; i++){
	for (int j = 0; j < children[i].size(); j++){
		//std::cout << children[i][j] << std::endl;
		node = children[i][j];
		for (int k = 0; k < N_states; k++){
			
			if ((adj[(node*N_states)+ k] == 1)){
				//std::cout << "k" << k << std::endl;
				neighbors[node].push_back(k);
				int genk = find_gen(children,k);
				if (genk == i){
					
					
					current_edges[node].push_back(k);
					
				}
				else if (genk > i){
					
					forward_edges[node].push_back(k);
					
				}
				else{


					back_edges[node].push_back(k);

				}
			
			}
		}

	}
	
}
	//print_vec_vec(neighbors);
}


 
 ///call only after build_array
bool compare(Model * m, Model * m1, Graph_Structure * g, Graph_Structure * g1){
	
	g->reset(m->G.N);
	g1->reset(m1->G.N);
	
	g->build_array1(m->G);
	g1->build_array1(m1->G);
	if (g->E != g1->E){
		
		return true; 
		
	}
	
	
	g->build_array2(m->G);
	g1->build_array2(m1->G);
	
	if(g->children[1].size() != g1->children[1].size()){
					
			return true;
				
		}
		
		
	g->build_array3();
	g1->build_array3();	
		
	if (g->gen != g1->gen){
		
		return true;
		
	}
		
	for (int i = 2; i <= g->gen; i++){
			
		if(g->children[i].size() != g1->children[i].size()){
					
			return true;
				
		}
				
	}		
				
	g->build_array4();
	g1->build_array4();

	 //all generation counts are the same
	
		vector<vector<int>> Edges;
		vector<vector<int>> Edges1; 
		
		
		Edges.resize(g->gen);
		Edges1.resize(g1->gen);
		for (int i = 1; i <= g->gen; i++){ 
			
			for (int j = 0; j < g->children[i].size(); j++){
					
					
						int node = g->children[i][j];
						int node2 = g1->children[i][j];
						
						Edges[i-1].push_back(g->back_edges[node].size());
						Edges[i-1].push_back(g->forward_edges[node].size());
						Edges[i-1].push_back(g->current_edges[node].size());
						 Edges1[i-1].push_back(g1->back_edges[node2].size());
						Edges1[i-1].push_back(g1->forward_edges[node2].size());
						Edges1[i-1].push_back(g1->current_edges[node2].size()); 
						
			}			
		}

			//std::cout << "Edges size" << Edges[0].size() << std::endl;

		
		
		
			if(!col_perm(Edges,Edges1,g->gen,g->children,g1->children))
				
				return true;
		
		
				/* else if (g->child_connections[i] != g1->child_connections[i]){
					
					return true;
					
				}
			
				else if (g->child_count[i].size() != g1->child_count[i].size()){
					
					return true;
					
				}
				
				else if (count_perm(g->child_count[i], g1->child_count[i],g->N_states)){
							
						vector<int> rowsum1;
						vector<int> rowsum2;
						int tempsum1;
						int tempsum2;			
							
						for (int i = 0; i < g->N_states; i++){

							tempsum1 = 0;
							tempsum2 = 0;
							for (int j = 0; j <g->N_states; j++){
				
				
								tempsum1 = tempsum1 + g->adj[(i*g->N_states)+j];
								tempsum2 = tempsum2 + g1->adj[(i*g->N_states)+j];
						
								
				
							}
								//std::cout << "tempsum1" << tempsum1 << std::endl;
								//std::cout << "tempsum2" << tempsum2 << std::endl;
								rowsum1.push_back(tempsum1);
								rowsum2.push_back(tempsum2);
						}						
									
						if	(!count_perm(rowsum1,rowsum2,g->N_states)){
							return true;
						}
						else{

							return false;

						}					
							
				}						 */
		
		else {
			
			return false;
			
		}
		
	
	
	
	
}

