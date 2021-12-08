#ifndef Matrix_HPP_
#define Matrix_HPP_

#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include <math.h>


//Writing a Matrix class to support matrix behavior (similar to MATLAB)


class Matrix {
	
	private:
	int rows;
	int columns;
	int format = -1; // 0: row major 1: column major
	public:
	std::vector<double> data;
	//int elements = rows*columns; //for nx1 vector only (for now)
	Matrix (int,int,int); //Constructor: initiate with some default
	Matrix();
	void populate_Matrix(int,int,int);
	int getrows();
	int getcolumns();
	int getformat();
	void setvalue(int);
	void addrow(int);
	void addcolumn(int);
	double matrix_sum();
	void operator+= (double);
	void operator-= (double);
	void operator*= (double);
	void operator/= (double);
	double& operator() (int,int);
		
	friend Matrix operator+ (Matrix& matrix1, Matrix& matrix2) {

		Matrix result(matrix1.rows,matrix1.columns,matrix1.format);

		if (matrix1.format == matrix2.format) {
			for (int i = 0; i < (matrix1.rows*matrix1.columns); i++) {
					result.data[i] = matrix1.data[i] + matrix2.data[i];
				}
			
			return result;
		}

		else { 
			std::cout << "Warning: You are trying to add in different formats!" << std::endl;
			return result;
		}
	}

	friend Matrix operator- (Matrix& matrix1, Matrix& matrix2) {

		Matrix result(matrix1.rows,matrix1.columns,matrix1.format);


		if (matrix1.format == matrix2.format) {
		for (int i = 0; i < (matrix1.rows*matrix1.columns); i++) {
			
			result.data[i] = matrix1.data[i] - matrix2.data[i];
		}
	
		return result;
		}

		else { std::cout << "Warning: You are trying to subtract in different formats!" << std::endl;
		return result;
		}
	}


	friend Matrix getTranspose(Matrix& mat) {

		Matrix result(mat.rows,mat.columns,mat.format);
	
		result = mat;
		result.rows = mat.columns;
		result.columns = mat.rows;

		if (mat.format == 0) {
			result.format = 1;

		}
		else {	
			result.format =0;
			}	
		return result;
	
	
	}
	friend Matrix power(Matrix& mat,double value) {

		Matrix result(mat.rows,mat.columns,mat.format);
	
		for (int i = 0; i < (mat.rows*mat.columns); i++) {
			
			result.data[i] = pow(mat.data[i],value);
		}

		
		return result;
	
	
	}
	
	friend Matrix dotProduct(Matrix& mat1, Matrix& mat2){

		if((mat1.format == mat2.format) && (mat1.columns == mat2.rows)){

			Matrix result(mat1.rows, mat2.columns, mat1.format);

			for (int i = 1; i <= mat1.rows; i++){

				for( int j = 1; j <= mat2.columns; j++) {
					for(int k = 1; k <= mat2.rows; k++){
						result(i,j) = mat1(i,k)*mat2(k,j) + result(i,j);

					}

				}
			}
			return result;

		}

		else {
			std::cout << "Error different formats!" << std::endl; 
			return mat1;
		}
	
	}


friend std::ostream& operator<< (std::ostream& out, Matrix& mat) {
		//std::ostringstream output;	// for writing to string
		out << "[";
		for (int i = 0; i < mat.rows; i++) {
			for(int j = 0; j < mat.columns; j++){
				out << mat.data[(i*mat.columns)+j] << " ";
			}
			out << std::endl;
		}
		out << "]" << std::endl;
		// transfer the contents of output to out
		return out;
	}

};
#endif