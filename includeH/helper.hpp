#include <vector>
#include <fstream>

template<typename T>
void print(std::vector<T> v,std::ofstream& out);

template<typename T>
void print(std::vector<std::vector<T>> v,std::ofstream& out);

template<typename T>
void print(std::vector<std::vector<std::vector<T>>> v,std::ofstream& out);

template <typename T>
void disp(std::vector<T> v);

template <typename T>
void disp(std::vector<std::vector<T>> v);

template <typename T>
void disp(std::vector<std::vector<std::vector<T>>> v);

template<typename T>
int min(std::vector<T> v);

template<typename T>
T max(std::vector<T> v);

template<typename T>
T sum(std::vector<T> v);

template<typename T>
std::vector<T> transpose(std::vector<T>& matrix, int rows, int cols);


#include "helper.tcc"