
template<typename T>
void print(std::vector<T> v,std::ofstream& out){

	for (auto i: v){
		out << i << "\t";
	}

	out << std::endl;
}

template<typename T>
void print(std::vector<std::vector<T>> v,std::ofstream& out){

	for (auto i: v){
		for(auto j : i){
			out << j << "\t";
		}
	}

	out << std::endl;
}

template<typename T>
void print(std::vector<std::vector<std::vector<T>>> v,std::ofstream& out){

	for (auto i: v){
		for(auto j : i){
			for(auto k : j){
				out << k << "\t";
			}
		}
	}

	out << std::endl;
}

template <typename T>
void disp(std::vector<T> v){

	for (auto i: v){
	
		std::cout << i << "\t";
	}

	std::cout << std::endl;
}

template <typename T>
void disp(std::vector<std::vector<T>> v){
	for (auto i: v){	
		for (auto j : i){
			std::cout << j << "\t";
		}
	}
	std::cout << std::endl;
}

template <typename T>
void disp(std::vector<std::vector<std::vector<T>>> v){
	for (auto i: v){	
		for (auto j : i){
			for(auto k : j){
				std::cout << k << "\t";
			}
		}
	}
	std::cout << std::endl;
}

template <class T>
int min(std::vector<T> v){
	T min = v[0];
	int min_idx = 0;
	for(int i = 1; i < v.size(); i++){
		if( v[i] < min){
			min = v[i];
			min_idx = i;
		}
	}
	return min_idx;
}

template <class T>
T max(std::vector<T> v){
	T max = v[0];
	for(int i = 1; i < v.size(); i++){
		if( v[i] > max){
			max = v[i];
		}
	}
	return max;
}

template<typename T>
T sum(std::vector<T> v){

	T sum = 0;

	for(int i = 0; i < v.size(); i++){
		sum += v[i];
	
	}

	return sum;
}

template <class T>
std::vector<T> transpose(std::vector<T>& matrix, int rows, int cols){
	std::vector<double> new_matrix (matrix.size());
	for(int i = 0; i < rows; i++){	
		for( int j = 0; j < cols; j++){		
		new_matrix[j*rows+i] = matrix[i*cols+j];		
		}
	}
	return new_matrix;
}

template<typename T>
int find(std::vector<T> v, T x){

	for(int i = 0; i < v.size(); i++){
		if(v[i] == x) return i;
	
	}

	return v.size();
}

template<typename T>
int contains(std::vector<T> v, T x){

	for(int i = 0; i < v.size(); i++){
		if(v[i] == x) return 1;
	
	}

	return 0;
}

