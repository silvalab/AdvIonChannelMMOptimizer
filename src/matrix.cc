#include "matrix.hpp"

// for functions in a template suported cass al declarations must have the syntax:
// template <clas T> ClassName:: ffname (parameters)
Matrix:: Matrix (int rows, int columns, int format){
	this-> rows = rows;
	this-> columns = columns;
	this-> format = format;
	data.resize(rows*columns,0);
	
}
Matrix::Matrix(){
    rows = 0;
	columns = 0;
	format = -1;
}

void Matrix::populate_Matrix(int rows,int columns,int format){
	
	this->rows = rows;
	this->columns = columns;
	this->format = format;
	data.resize(rows*columns,0);
	
}
int Matrix:: getrows () { return rows;}

int Matrix:: getcolumns() {return columns;}

int Matrix:: getformat() { return format;}
	 
void Matrix :: setvalue(int val) {
	fill(data.begin(), data.end(), val);
}

void Matrix::addrow(int location) {
	if (format== 0){
		std::vector<double>:: iterator it;
		it = data.begin();	
		it = data.insert(it+((location-1)*columns),0);
		data.insert(it,columns-1,0);
		rows = rows + 1;
	}


	else{

		std::vector<double> :: iterator it;
		it = data.begin();
		it = data.insert(it+location-1,0);

		for (int i = 0; i < columns-1; i++) {
			it = data.insert(it+columns-1,0);
		}		
		rows= rows+1;
	}
} 

void Matrix::addcolumn(int location) {
	if (format == 0){// row major
		std::vector<double> :: iterator it;
		it = data.begin();
		it = data.insert(it+(location-1),0);
	
		for (int i = 0; i < rows-1; i++) {
			it = data.insert(it+columns+1,0);
		}

		columns= columns+1;

	}


	else{
		std::vector<double> :: iterator it;
		it = data.begin();
		it = data.insert(it+((location-1)*rows),0);
		data.insert(it,rows-1,0);
		columns= columns+1;
	}
} 


double Matrix::matrix_sum(){
	
	double sum = 0;
	for(int i = 0; i < data.size(); i++){
		sum += data[i];
		
	}	
	return sum;
}


void Matrix:: operator+= (double value) {

	for (auto &q: data) q+= value; }

void Matrix:: operator-= (double value) {

	for (auto &q: data) q-= value; }

void Matrix:: operator*= (double value) {

	for (auto &q: data) q*= value; }

void Matrix:: operator/= (double value) {

	for (auto &q: data) q/= value; }


double& Matrix::operator() (int a, int b) {
	if (format== 0){
		return data[(columns*a) + b];
	}
	/*  else{
		return 0;
	}  */
}


/* int main () {
	Matrix<double> matrix1 (5,3,0);
	Matrix<double> matrix2 (5,3,1);
	Matrix<double> result (3,5,1);
	Matrix<double> transpose (3,5,1);
	Matrix<double> product (3,3,1);
	//cout << "Original data: [";
	//for (int i =0; i < matrix1.getrows(); i++) cout << matrix1[i]
	//	<< " ";
	//cout << "]" << endl;
//	matrix1 += 10;
//	cout << "[";
//	for (int i =0; i <matrix1.getrows(); i++)
//	cout << matrix1[i] << " ";
//	cout << "]" << endl;

	//cout << "Data:" << matrix1 << endl;
	//cout << matrix1.getrows() << endl;	
	//cout << matrix1.getcolumns() << endl;
	//cout << matrix1.getformat() << endl;
	matrix1.setvalue(1);
	matrix2.setvalue(2);
	matrix1(1,2) = 6;
	matrix2(1,2) = 3;
	//result = matrix1-matrix2;
	//cout << result << endl;
	//cout << "Original Matrix1" << matrix1 << endl;
	transpose = getTranspose(matrix1);
	cout << "Transpose" << transpose << endl;
	cout << "Original Matrix2" << matrix2 << endl;
	product = dotProduct(transpose,matrix2);
	cout << product << endl;
	return 0;
	} */
