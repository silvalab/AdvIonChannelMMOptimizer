#ifndef helper_HPP_
#define helper_HPP_
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

template<typename T>
int find(std::vector<T> v, T x);

template<typename T>
int contains(std::vector<T> v, T x);

#include "helper.tcc"

#endif