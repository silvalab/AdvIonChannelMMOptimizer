#include "cluster.hpp"
#include "model.hpp"



template <class T, class S>
int min_vec(std::vector<T> v,std::vector<S> distances ){
	 S min = distances[v[0]];
	// std::cout << "starting min" << min << std::endl;
	 int min_idx = 0;
	 for(int i = 0; i < v.size(); i++){
		if( distances[v[i]] <= min){ 
			min_idx = i;
			min = distances[v[i]];
		}
	 }
	// std::cout << min_idx << std::endl;
	 return min_idx;
}


template <class T,class S> 
void print_vec_vec_vec(std::vector<std::vector<std::vector<T>>> v, std::vector<S> distances){
	
	for (int i = 0; i < v.size(); i++){
		for (int j = 0; j < v[i].size(); j++){
			for(int k = 0; k < v[i][j].size(); k++){
					std::cout << distances[v[i][j][k]] << "\t";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}


template <class T,class S> 
int min_cluster_indx(std::vector<std::vector<std::vector<T>>> v, std::vector<S> distances){
	S min = distances[v[0][0][0]];
	int min_idx = 0;
	for (int i = 0; i < v.size(); i++){
		if (distances[v[i][0][0]] < min){
			min = distances[v[i][0][0]];
			min_idx = i;
		}
	}
	
	return min_idx;
}


template <class T>
void print_vec(std::vector<T> v){
	
	for (int i = 0; i < v.size(); i++){
		
		std::cout << v[i] << "\t";
		
	}
	std::cout << std::endl;
}

template <class S>
std::vector<int> min_vec_multiple(std::vector<int> cluster,std::vector<S> distances,int min_costs_to_add){
	
	std::vector<int> values_to_return (min_costs_to_add,0);
	//print_vec(cluster);
	//print_vec(distances);
	std::sort(cluster.begin(),cluster.end(), [&](int a, int b){ return distances[a] < distances[b];});
	//print_vec(cluster);
	for (int i = 0; i < min_costs_to_add; i++){
		values_to_return[i] = cluster[i];
	}
	
	return values_to_return;
}

template <class T>
void print_vec_vec(std::vector<std::vector<T>> v){
	
	for (int i = 0; i < v.size(); i++){
		
		for (int j = 0; j < v[i].size(); j++){
			
			std::cout << v[i][j] << "\t";
			
		}
		std::cout << std::endl;
	}
}


template <class T>
void print_list_list(std::list<std::list<T>> v){
	
	for (int i = 0; i < v.size(); i++){
		
		for (int j = 0; j < v[i].size(); j++){
			
			std::cout << v[i][j] << "\t";
			
		}
		std::cout << std::endl;
	}
}



template <class T>
void min_dist_clust(std::vector<std::vector<T>> v, int& i, int& j){
	 T min = v[0][1];
	for (int k = 0; k < v.size(); k++){
		if(v[k][1] <= min){
				i = k;
		}
	} 
	
	j = v[i][0];
	//std::cout << "i" << "\t" <<  i << "\t" << "j" << "\t" << j << std::endl;
}

template <class T> 
void new_D(std::vector<std::vector<T>>& D, int i, int j){
	
	//replace row i by min of row i and row j
	for (int k = 0; k < D.size(); k++){
		if(i != k){
			D[i][k] = std::min(D[i][k],D[j][k]);
		}
	}
	
	//infinity out row j 
	for (int k = 0; k < D.size(); k++){
		D[j][k] = inf ;
	}
	
	//infinity out column j
	for (int k = 0; k < D.size(); k++){
		D[k][j] = inf ;
	}
}


//void update_dmin

template <class T>
int min_idx(std::vector<T> v, int i){
	T min_dist = inf;
	int min_idx = 0;
	for (int k = 0; k < v.size(); k++){
		if(k != i){
			if(v[k] < min_dist ){
				min_dist = v[k];
				min_idx = k;
			}
		}
	}
	return min_idx;
}

template <class T>
void populate_dmin(std::vector<std::vector<T>>& dmin, std::vector<std::vector<T>> D){
	for(int i = 0; i < D.size(); i++){ //populate dmin
		dmin[i][0] = min_idx(D[i],i);
		dmin[i][1] = D[i][dmin[i][0]]; 
	}
}

int getCluster(std::vector<int> ind){
	
	for (int i = 0; i < ind.size(); i++){
		if( ind[i] == 0) return 1; // needs to run again one member has not been categoriized
	}
	return 0;
}

int ij_presence(std::vector<std::vector<std::vector<int>>> clusters,int i,int j){
	std::vector<int> new_ij(4, 0); //(i,j,both,neither)
	for(int k = 0; k < clusters.size(); k++){
				for(int r = 0; r < clusters[k].size(); r++){
					for (int s = 0; s < clusters[k][r].size(); s++){
						if( clusters[k][r][s] == i){ // j is new
							new_ij[0]++;				
						}
						else if(clusters[k][r][s] == j){ // i is  new	
							new_ij[1]++;
						}
					}
				}
			}
	if(new_ij[0] == 0 && new_ij[1] == 0){ //both i and j are new
		return 2;
	}
	else if(new_ij[0] == 1 && new_ij[1] == 1){ //neither i and j ara new
		return 3;
	}
	else if(new_ij[1] == 1){ // j, is present, i is new
		return 0;
	}
	else {
		return 1; // i is present, j is new
	}
}

void populate_clusters(std::vector<std::vector<std::vector<int>>>& clusters, int i, int j){
	for(int k = 0; k < clusters.size(); k++){
				for(int r = 0; r < clusters[k].size(); r++){
					for (int s = 0; s < clusters[k][r].size(); s++){
						if(ij_presence(clusters,i,j) == 0){ // i is new, j in dendrogram
							if(clusters[k][r][s] == j){
								clusters[k].push_back(std::vector<int> (1,i));
								return;
							}
						}
						else{ // j is new, i in dendrogram
							if(clusters[k][r][s] == i){
								clusters[k].push_back(std::vector<int> (1,j));
								return;
							}
							
						}
					}
				}
				
			}
}

void extractab(std::vector<std::vector<std::vector<int>>> clusters,int i,int j,int& a, int& b){
	int extracta = 1;
	for(int r = 0; r < clusters.size(); r++){
		for(int s = 0; s < clusters[r].size(); s++){
			for(int t = 0; t < clusters[r][s].size(); t++){
				if(extracta){
					if(clusters[r][s][t] == i || clusters[r][s][t] == j){
						a = r;
						extracta = 0;
					}
				}	
				else if(clusters[r][s][t] == i || clusters[r][s][t] == j){
					b = r;
					return;
				}
			}
		}
	}
}


void compress(std::vector<std::vector<std::vector<int>>>& clusters, int a, int b){
	
	// a guaranteed to  be less than b
	//merge the root
	for (int i = 0; i < clusters[b][0].size(); i++){
		clusters[a][0].push_back(clusters[b][0][i]);
	}
	for (int i = 1; i < clusters[b].size(); i++){
		
		clusters[a].push_back(clusters[b][i]);
	}
	
	//delete clusters[b]
	
	std::vector<std::vector<std::vector<int>>> clusters_temp;
	
	for (int i = 0; i < clusters.size(); i++){
		if(i != b){	
			clusters_temp.push_back(clusters[i]);
		}
	}	
	clusters.swap(clusters_temp);
}

void merge_clusters(std::vector<std::vector<std::vector<int>>>& clusters, int i, int j){
	int a,b;
	
	extractab(clusters,i,j,a,b);
	
	//std::cout << "a" << "\t" << a << "b" << "\t" << b << std::endl;
	
	//compress clusters between a,b 
	compress(clusters,a,b);
		
}
	
void single_linkage_cluster(std::vector<Model*>& models, std::vector<double>& distances, int HL, std::vector<ProtocolParameter> protos, int a_j){
	//std::vector<std::vector<Model*>>& models, int HL
	
	while (distances.size() >= HL){
		int N = distances.size();
		std::vector<std::vector<double>> D (N,std::vector<double>(N));
		std::vector<std::vector<double>> dmin (N, std::vector<double> (2));
		std::vector<std::vector<std::vector<int>>> clusters;
		//double inf = std::numeric_limits<double>::infinity();
		std::vector<int> ind (N,0);
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				if(i !=j){
					float d = std::abs(distances[i]-distances[j]);
					D[i][j] = d;
					D[j][i] = d;
				}
				else{
					D[i][j] = inf;
				}
			}
		}
	
		//print_vec_vec(D);
		//print_vec_vec(dmin);
	
	
		populate_dmin(dmin,D);
		//print_vec_vec(dmin);
	
		int i,j;
	
		min_dist_clust(dmin,i,j);
		if( j < i ){
			int itemp = i;
			i = j;
			j = itemp;
		}
		//std::cout << "i" << i << std::endl; 
		//std::cout << "j" << j << std::endl;  
	
		std::vector<int> clust;
		clust.push_back(i);
		clust.push_back(j);
		std::vector<std::vector<int>> cluster;
		cluster.push_back(clust);
		clusters.push_back(cluster);
		print_vec_vec_vec(clusters,distances);
		ind[i]++;
		ind[j]++;
		while(getCluster(ind)){
			new_D(D,i,j);
			//print_vec_vec(D);
			populate_dmin(dmin,D);
			//print_vec_vec(dmin);
	
			min_dist_clust(dmin,i,j); 
			if( j < i ){
				int itemp = i;
				i = j;
				j = itemp;
			}
			//std::cout << "i" << i << std::endl; 
			//std::cout << "j" << j << std::endl; 
			ind[i]++; ind[j]++;
			int new_ind = ij_presence(clusters,i,j);
			if(new_ind == 2){
				std::vector<int> clust;
				clust.push_back(i);
				clust.push_back(j);
				std::vector<std::vector<int>> cluster;
				cluster.push_back(clust);
				clusters.push_back(cluster);
			}
			else if(new_ind == 3){ // i and j are already in the dendrogram
				//std::cout << "i and j are in already in the dendrogram" << std::endl;
				merge_clusters(clusters,i,j);
			}
			else{ // either i or j is already in the dendrogram
				//find which
				populate_clusters(clusters,i,j);
			}
			//print_vec_vec_vec(clusters,distances);
		}
	
		//std::cout << "final cluster vec" << std::endl;
		//print_vec_vec_vec(clusters,distances);
		int numClusters = clusters.size();
		int min_costs_to_add = HL - numClusters;
		//std::cout << "min_coststoadd" << min_costs_to_add << std::endl;
		std::vector<double> distances_to_return;
		std::vector<Model*> models_to_return;
		std::vector<std::vector<int>> rep_clusters (clusters.size());
		int min_idx = min_cluster_indx(clusters,distances);
		//std::cout << "min_idx" << min_idx << std::endl;
		for(int i = 0; i < clusters.size(); i++){
			std::vector<int> cluster;
			//std::cout << "i" << i << std::endl;
			for(int j = 0; j < clusters[i].size(); j++){
				for(int k = 0; k < clusters[i][j].size(); k++){
					cluster.push_back(clusters[i][j][k]);
				}			
			}
			if( i == min_idx){ //min cluster
				if( min_costs_to_add < cluster.size()){
					std::vector<int> to_add = min_vec_multiple(cluster,distances,min_costs_to_add);
					//std::cout << "to_add" << std::endl;
					//print_vec(to_add);
					for (int s = 0; s < to_add.size(); s++){
						rep_clusters[i].push_back(to_add[s]); //take the first min_costs to _add of min cluster to reach HL
					}
				}
				else{ //just take the cluster
					for ( int s = 0; s < cluster.size(); s++){
						rep_clusters[i].push_back(cluster[s]);
					}
				}
					
			}
			else{
				rep_clusters[i].push_back(cluster[min_vec(cluster,distances)]); 
			}			//the index of each cluster to save vec length equal to number of clusters
		}
		//std::cout << "rep clusters" << std::endl;
		//print_vec_vec(rep_clusters);
		
		for (int i = 0; i < rep_clusters.size(); i++){
			for(int j = 0; j < rep_clusters[i].size(); j++){
				//std::cout << distances[rep_clusters[i][j]] << std::endl;
				distances_to_return.push_back(distances[rep_clusters[i][j]]);
				models_to_return.push_back(new Model(models[rep_clusters[i][j]]));
			}
		}
	
		for(int i = 0; i < models.size(); i++){
			delete (models[i]);
		}
		distances.swap(distances_to_return);
		models.swap(models_to_return);
	}
	
}

void single_linkage_cluster(std::vector<Model*>& models, std::vector<double>& distances, int HL, std::vector<ProtocolParameter> protos){
	//std::vector<std::vector<Model*>>& models, int HL
	
	while (distances.size() >= HL){
		int N = distances.size();
		std::vector<std::vector<double>> D (N,std::vector<double>(N));
		std::vector<std::vector<double>> dmin (N, std::vector<double> (2));
		std::vector<std::vector<std::vector<int>>> clusters;
		//double inf = std::numeric_limits<double>::infinity();
		std::vector<int> ind (N,0);
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				if(i !=j){
					float d = std::abs(distances[i]-distances[j]);
					D[i][j] = d;
					D[j][i] = d;
				}
				else{
					D[i][j] = inf;
				}
			}
		}
	
		//print_vec_vec(D);
		//print_vec_vec(dmin);
	
	
		populate_dmin(dmin,D);
		//print_vec_vec(dmin);
	
		int i,j;
	
		min_dist_clust(dmin,i,j);
		if( j < i ){
			int itemp = i;
			i = j;
			j = itemp;
		}
		//std::cout << "i" << i << std::endl; 
		//std::cout << "j" << j << std::endl;  
	
		std::vector<int> clust;
		clust.push_back(i);
		clust.push_back(j);
		std::vector<std::vector<int>> cluster;
		cluster.push_back(clust);
		clusters.push_back(cluster);
		//print_vec_vec_vec(clusters,distances);
		ind[i]++;
		ind[j]++;
		while(getCluster(ind)){
			new_D(D,i,j);
			//print_vec_vec(D);
			populate_dmin(dmin,D);
			//print_vec_vec(dmin);
	
			min_dist_clust(dmin,i,j); 
			if( j < i ){
				int itemp = i;
				i = j;
				j = itemp;
			}
			//std::cout << "i" << i << std::endl; 
			//std::cout << "j" << j << std::endl; 
			ind[i]++; ind[j]++;
			int new_ind = ij_presence(clusters,i,j);
			if(new_ind == 2){
			//	std::cout << " either i and j are new" << std::endl;
				std::vector<int> clust;
				clust.push_back(i);
				clust.push_back(j);
				std::vector<std::vector<int>> cluster;
				cluster.push_back(clust);
				clusters.push_back(cluster);
			}
			else if(new_ind == 3){ // i and j are already in the dendrogram
				//std::cout << "i and j are in already in the dendrogram" << std::endl;
				merge_clusters(clusters,i,j);
			}
			else{ // either i or j is already in the dendrogram
				//find which
				//std::cout << " either i and j are in already in the dendrogram" << std::endl;
				populate_clusters(clusters,i,j);
			}
			//print_vec_vec_vec(clusters,distances);
		}
	
		//std::cout << "final cluster vec" << std::endl;
		//print_vec_vec_vec(clusters,distances);
		int numClusters = clusters.size();
		int min_costs_to_add = HL - numClusters;
		//std::cout << "min_coststoadd" << min_costs_to_add << std::endl;
		std::vector<double> distances_to_return;
		std::vector<Model*> models_to_return;
		std::vector<std::vector<int>> rep_clusters (clusters.size());
		int min_idx = min_cluster_indx(clusters,distances);
		//std::cout << "min_idx" << min_idx << std::endl;
		for(int i = 0; i < clusters.size(); i++){
			std::vector<int> cluster;
			//std::cout << "i" << i << std::endl;
			for(int j = 0; j < clusters[i].size(); j++){
				for(int k = 0; k < clusters[i][j].size(); k++){
					cluster.push_back(clusters[i][j][k]);
				}			
			}
			if( i == min_idx){ //min cluster
				if( min_costs_to_add < cluster.size()){
					std::vector<int> to_add = min_vec_multiple(cluster,distances,min_costs_to_add);
					//std::cout << "to_add" << std::endl;
					//print_vec(to_add);
					for (int s = 0; s < to_add.size(); s++){
						rep_clusters[i].push_back(to_add[s]); //take the first min_costs to _add of min cluster to reach HL
					}
				}
				else{ //just take the cluster
					for ( int s = 0; s < cluster.size(); s++){
						rep_clusters[i].push_back(cluster[s]);
					}
				}
					
			}
			else{
				rep_clusters[i].push_back(cluster[min_vec(cluster,distances)]); 
			}			//the index of each cluster to save vec length equal to number of clusters
		}
		//std::cout << "rep clusters" << std::endl;
		//print_vec_vec(rep_clusters);
		
		for (int i = 0; i < rep_clusters.size(); i++){
			for(int j = 0; j < rep_clusters[i].size(); j++){
				//std::cout << distances[rep_clusters[i][j]] << std::endl;
				distances_to_return.push_back(distances[rep_clusters[i][j]]);
				models_to_return.push_back(new Model(models[rep_clusters[i][j]]));
			}
		}
	
		for(int i = 0; i < models.size(); i++){
			delete (models[i]);
		}
		distances.swap(distances_to_return);
		models.swap(models_to_return);
	}
	
}

/* int main() {
	//0.75, 173.25, 185.7, 90.85, 621.75, 555.20, 115.23
	std::vector<float> distances = {695.59,695.59,695.722, 695.215, 206.618, 675.286, 675.286, 59.7941, 233.784,514.024,7.719,28.9866,55.371,54.4328,4.9,112.92,445.268};
	auto clusters = single_linkage_cluster(distances, 10);
	
	print_vec(clusters);
	
	
	return 0;
} */