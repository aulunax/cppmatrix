#include "Matrix.h"
#include "MatrixMisc.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <thread>
#include <future>



template <typename T>
Matrix<T>::Matrix(const std::string &content) : Matrix()
{
MEASURE_TIME_START
	std::istringstream contentStream(content);
	std::string rowStr;
	int m=0,n=0;
	n = Mat::countOccurrences(content,';') + 1;

    std::vector<std::future<void>> futures;

	int i = 0;
	while (std::getline(contentStream, rowStr, ';')) {

		m = Mat::countOccurrences(rowStr,',') + 1;

		if (size.m == 0) {
			reserve(n, m);
		}
		else if (m != size.m) { 	// check for correct number of commas
			size = {0,0};
			rawData = Data(); 
			std::cerr << "Invalid matrix format" << std::endl;
			return;
		}
		futures.push_back(std::async(std::launch::async, [this, rowStr, i]() {
			std::istringstream rowStream(rowStr);
			std::string cell;
			int j = 0;
			while (std::getline(rowStream, cell, ',')) {
				std::istringstream numberStream(cell);
				T value;
				numberStream >> value;
				rawData(i,j) = value;
				j++;
			}
		}));
		i++;
	}

	for (auto& future : futures) {
        future.wait();
    }

MEASURE_TIME_END
}

template <typename T>
Matrix<T>::Matrix(int n, int m)
{
    reserve(n,m);
}

template <typename T>
Matrix<T>::Matrix(int n, int m, T value)
{
	fill(n, m, value);
}

template<typename T>
Matrix<T>::Matrix(Dimensions size) {
	reserve(size.n, size.m);
}

template <typename T>
Matrix<T>::Matrix(Dimensions size, T value)
{
	fill(size.n, size.m, value);
}

template<typename T>
Matrix<T>::Matrix(const Matrix& other)
{
	reserve(other.size.n, other.size.m);
    if (this != &other) {
		rawData = other.rawData;
    }
}


template<typename T>
Matrix<T>::Matrix(Matrix&& other)
{
	size = other.size;
    if (this != &other) {
		other.size = Dimensions{0,0};
        rawData = std::move(other.rawData);
    }
}

template <typename T>
void Matrix<T>::print()
{
	if (size == Dimensions{0,0}) {
		std::cout << "Uninitialized Matrix\n\n";
		return;
	}

	for (int i = 0; i < size.n; i++) {
		for (int j = 0; j < size.m; j++) {
			std::cout << rawData(i,j) << ' ';
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}

template <typename T>
void Matrix<T>::reserve(int n, int m)
{
MEASURE_TIME_START
	size = {n,m};
	rawData = Data(n, m);
MEASURE_TIME_END
}

template <typename T>
void Matrix<T>::fill(int n, int m, T value)
{
MEASURE_TIME_START
	size = {n,m};
	rawData = Data(n, m, value);
MEASURE_TIME_END
}

template <typename T>
Matrix<T> Matrix<T>::inv() const
{
	if (size.n != size.m)
		throw MatrixSizeDisparityException();

MEASURE_TIME_START
	bool isDiagonal = true;
	for (int i = 0; i < size.n; i++) {
		for (int j = 0; j < size.m; j++) {
			if (i != j && rawData(i,j) != 0) {
				isDiagonal = false;
				break;
			}
			if (i == j && rawData(i,j) == 0) {
				throw MatrixNotInvertible();
			}		
		}
		if (!isDiagonal)
			break;
	}

	Matrix<T> result(size, 0);
	if (isDiagonal) {
		for (int i = 0; i < size.n; i++) {
			result(i,i) = 1/rawData(i,i);
		}
	}
MEASURE_TIME_END
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::diag(int k) const
{
	if (size.m != size.n) {
		throw NotSquareMatrixException();
	}
MEASURE_TIME_START
	Matrix<T> d(size, 0); 
	int absK = abs(k);
	int startRow = k < 0 ? absK : 0;
	int startCol = k > 0 ? absK : 0;
	for (int j = 0; j+absK < size.m ; j++) {
		d(j+startRow,j+startCol) = rawData(j+startRow,j+startCol);
	}
MEASURE_TIME_END
	return d;
	// Matrix<T> args = Matrix<T>(1, 1, k);
	// Matrix<T> other;
	// return returnTypeFunctionWrapper("diag(int k)", &Matrix::threadedMatrixOperation, this, args, other, &Matrix::getDiagonal, Dimensions{0,0});

	//return threadedMatrixOperation(args, other, &getDiagonal);
}

template <typename T>
void Matrix<T>::getDiagonal(const Data& args, const Data &data1, const Data &data2, Data &result, int startRow, int endRow)
{
	int diagonal = args(0,0);
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data1.numCols(); j++) {
			if (i+diagonal == j) 
				result(i,j) = data1(i,j);
			else
				result(i,j) = 0;
		}
	}
}

template <typename T>
Matrix<T> Matrix<T>::tril(int k) const
{
	if (size.m != size.n) {
		throw NotSquareMatrixException();
	}

    Matrix<T> args = Matrix<T>(1,2);
	args(0,0) = k;
	args(0,1) = -1;

	Matrix<T> other;
	return returnTypeFunctionWrapper("tril(int k)", &Matrix::threadedMatrixOperation, this, args, other, &Matrix::getTriangle, Dimensions{0,0});

	//return threadedMatrixOperation(args, other, &getTriangle);
}

template <typename T>
Matrix<T> Matrix<T>::triu(int k) const
{
	if (size.m != size.n) {
		throw NotSquareMatrixException();
	}

    Matrix<T> args = Matrix<T>(1,2);
	args(0,0) = k;
	args(0,1) = 1;

	Matrix<T> other;
	return returnTypeFunctionWrapper("triu(int k) ", &Matrix::threadedMatrixOperation, this, args, other, &Matrix::getTriangle, Dimensions{0,0});

	//return threadedMatrixOperation(args, other, &getTriangle);
}

template <typename T>
void Matrix<T>::getTransposed(const Data &args, const Data &data1, const Data &data2, Data &result, int startRow, int endRow)
{
	for (int i = startRow; i < endRow; ++i) {
        for (int j = 0; j < data1.numCols(); ++j) {
            result(j,i) = data1(i,j);
        }
    }
}

template <typename T>
Matrix<T> Matrix<T>::transpose() const
{
	Matrix<T> args;
	Matrix<T> other;
	return returnTypeFunctionWrapper("transpose()", &Matrix::threadedMatrixOperation, this, args, other, &Matrix::getTransposed,  Dimensions{size.m, size.n});

	//return threadedMatrixOperation(args, other, &getTransposed, Dimensions{size.m, size.n});
}

template <typename T>
void Matrix<T>::getTriangle(const Data& args, const Data &data1, const Data &data2, Data &result, int startRow, int endRow)
{
	int direction = args(0,1);
	int k = args(0,0);
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data1.numCols(); j++) {
			if (direction == 1 ? i+k <= j : i+k>=j) 
				result(i,j) = data1(i,j);
			else
				result(i,j) = 0;
		}
	}
}



template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix& other) {
	reserve(other.size.n, other.size.m);
    if (this != &other) {
		rawData = other.rawData;
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(Matrix&& other) noexcept {
	size = other.size;
    if (this != &other) {
		other.size = Dimensions{0,0};
        rawData = std::move(other.rawData);
    }
    return *this;
}


template <typename T>
void Matrix<T>::fillWithValue(const Data &args, const Data &data1, const Data &data2, Data &result, int startRow, int endRow)
{
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < result.numCols(); j++) {
			result(i,j) = data2(0,0);
		}
	}
}

template <typename T>
void Matrix<T>::addMatrices(const Data& args, const Data &data1, const Data &data2, Data &result, int startRow, int endRow)
{
    for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data1.numCols(); j++) {
			result(i,j) = data1(i,j) + data2(i,j);
		}
	}
}

template<typename T>
void Matrix<T>::substractMatrices(const Data& args, const Data& data1, const Data& data2, Data& result, int startRow, int endRow) {
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data1.numCols(); j++) {
			result(i,j) = data1(i,j) - data2(i,j);
		}
	}
}

template<typename T>
void Matrix<T>::multiplyMatrices(const Data& args, const Data& data1, const Data& data2, Data& result, int startRow, int endRow) {
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data2.numCols(); j++) {
			for (int k = 0; k < data1.numCols(); k++) {
				result(i,j) += data1(i,k) * data2(k,j);
			}
		}
	}
}

template<typename T>
void Matrix<T>::multiplyByConstant(const Data& args, const Data& data1, const Data& data2, Data& result, int startRow, int endRow) {
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data1.numCols(); j++) {
			result(i,j) = data1(i,j) * args(0,0);
		}
	}
}

template<typename T>
void Matrix<T>::divideByConstant(const Data& args, const Data& data1, const Data& data2, Data& result, int startRow, int endRow) {
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data1.numCols(); j++) {
			result(i,j) = data1(i,j) / args(0,0);
		}
	}
}

template<typename T>
Matrix<T> Matrix<T>::threadedMatrixOperation(const Matrix& args, const Matrix& other, ThreadedFunction operation, Dimensions resultSize) const {
	Matrix<T> temp;
	if (resultSize != Dimensions{0,0}) {
		temp = Matrix<T>(resultSize);
	}
	else if (size.m == other.size.n)
		temp = Matrix<T>(size.n, other.size.m);
	else
		temp = Matrix<T>(size);


	int numOfThreads = (int)std::thread::hardware_concurrency();

	// if matrix has only few rows, don't multithread
	if (size.n > numOfThreads*2 && size.m > 1) {

		int rowsPerThread = size.n / numOfThreads;

		std::vector<std::thread> threads;
		for (int i = 0; i < numOfThreads - 1; ++i) {
			threads.emplace_back(operation, std::ref(args.rawData), std::ref(rawData), std::ref(other.rawData), std::ref(temp.rawData), i * rowsPerThread, (i + 1) * rowsPerThread);
		}
		threads.emplace_back(operation, std::ref(args.rawData), std::ref(rawData), std::ref(other.rawData), std::ref(temp.rawData), (numOfThreads - 1) * rowsPerThread, size.n);


		for (auto& thread : threads) {
			thread.join();
		}
	}
	else {
		operation(std::ref(args.rawData), std::ref(rawData), std::ref(other.rawData), std::ref(temp.rawData), 0, size.n);
	}

	return temp;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix& other) const
{
	if (size != other.size) {
		throw MatrixSizeDisparityException();
	}
	Matrix<T> args;

	return returnTypeFunctionWrapper("operator+(const Matrix& other)", &Matrix::threadedMatrixOperation, this, args, other, &Matrix::addMatrices, Dimensions{0,0});
	//return threadedMatrixOperation(args, other, &addMatrices);
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix& other) const
{
	if (size != other.size) {
		throw MatrixSizeDisparityException();
	}
	Matrix<T> args;
	return returnTypeFunctionWrapper("operator-(const Matrix& other)", &Matrix::threadedMatrixOperation, this, args, other, &Matrix::substractMatrices, Dimensions{0,0});
	//return threadedMatrixOperation(args, other, &substractMatrices);
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix& other) const
{
	if (size.m != other.size.n) {
		throw MatrixSizeDisparityException();
	}
	Matrix<T> args;

	return returnTypeFunctionWrapper("operator*(const Matrix& other)", &Matrix::threadedMatrixOperation, this, args, other, &Matrix::multiplyMatrices, Dimensions{0,0});
	//return threadedMatrixOperation(args, other, &multiplyMatrices);
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T& value) const
{
	Matrix<T> args = Matrix<T>(1, 1, value);
	Matrix<T> other;
	return returnTypeFunctionWrapper("operator*(const T& value)", &Matrix::threadedMatrixOperation, this, args, other, &Matrix::multiplyByConstant, Dimensions{0,0});
	//return threadedMatrixOperation(args, other, &multiplyByConstant);
}

template<typename T>
Matrix<T> Matrix<T>::operator/(const T& value) const
{
	if (value == 0) {
		throw std::runtime_error("Division by zero error");
	}
	Matrix<T> args = Matrix<T>(1, 1, value);
	Matrix<T> other;
	return returnTypeFunctionWrapper("operator/(const T& value)", &Matrix::threadedMatrixOperation, this, args, other, &Matrix::divideByConstant, Dimensions{0,0});
	//return threadedMatrixOperation(args, other, &divideByConstant);
}

template <typename T>
Matrix<T> Matrix<T>::operator-() const
{
    return (*this)*-1;
}

template <typename T>
bool Matrix<T>::HasUniqueSolution() const
{
	int ranks = 0;
	for (int i = 0; i < size.n; i++) {
		bool zeros = true;
		for (int j = 0; j < size.m; j++) {
			if (rawData(i,j) != 0) {
				zeros = false;
				break;
			}
		}
		if (!zeros)
			ranks += 1;
	}
	if (ranks != size.n)
    	return false;
	return true;
}

template <typename T>
bool Matrix<T>::isTril() const
{
MEASURE_TIME_START
	for (int i = 0; i < size.n; i++) {
		for (int j = i+1; j < size.m; j++) {
			if (rawData(i,j) != 0) {
				MEASURE_TIME_END
				return false;
			}
		}
	}
MEASURE_TIME_END
	return true;
}

template <typename T>
bool Matrix<T>::isTriu() const
{
MEASURE_TIME_START
   for (int i = 0; i < size.n; i++) {
		for (int j = 0; j < i; j++) {
			if (rawData(i,j) != 0) {
				MEASURE_TIME_END
				return false;
			}
		}
	}
MEASURE_TIME_END
	return true;
}

template <typename T>
Matrix<T> Matrix<T>::backSubstitution(const Matrix &b) const
{
	Matrix<T> x(size.n, 1);

    for (int i = size.n - 1; i >= 0; i--) {
        x(i,0) = b(i,0);

        // Subtract the contributions of already solved variables
        for (int j = i + 1; j < size.n; j++) {
            x(i,0) -= rawData(i,j) * x(j,0);
        }

        // Divide by the diagonal element
        x(i,0) /= rawData(i,i);
    }
    return x;
}

template <typename T>
Matrix<T> Matrix<T>::forwardSubstitution(const Matrix &b) const
{
	Matrix<T> x(size.n, 1);

	for (int i = 0; i < size.n; i++) {
		x(i,0) = b(i,0);

		// Subtract the contributions of already solved variables
		for (int j = 0; j < i; j++) {
			x(i,0) -= rawData(i,j) * x(j,0);
		}

		// Divide by the diagonal element
		x(i,0) /= rawData(i,i);
	}
    return x;
}

template <typename T>
Matrix<T> Matrix<T>::operator/(const Matrix &other) const
{
MEASURE_TIME_START

	// transpose xA = b into  A^T * x^T = b^t
	// and perform left division
	Matrix<T> xTransposed = (other.transpose() | (*this).transpose());

	// tranpose result to get x from x^T
	return xTransposed.transpose();
}

template <typename T>
Matrix<T> Matrix<T>::operator|(const Matrix &other) const
{
#ifdef MATRIX_DEBUG
	std::cout<<"operator| starts\n";
#endif
MEASURE_TIME_START
	// check if it is a triangular matrix
	bool tril = isTril();
	bool triu = isTriu();


	Matrix<T> x(size.n, 1);

	// case 1: is lower triangular matrix
	if (tril) {
		x = forwardSubstitution(other);
	}
	// case 2: is upper triangular matrix
	else if (triu) {
		x = backSubstitution(other);
	}
	// case 3: is not triangular matrix
	// https://en.wikipedia.org/wiki/LU_decomposition#MATLAB_code_example
	else {
		Matrix<T> U = (*this);
		Matrix<T> L = Matrix<T>(size,1).diag();
		for(int i = 1; i < size.n; i++) {
			for(int j = 0; j < i; j++) {
				L(i,j) = U(i,j) / U(j,j);
				for(int k = 0; k < size.n; k++) {
					U(i,k) = U(i,k) - L(i,j)*U(j,k);
				}
			}
		}
		
		Matrix<T> y = L.forwardSubstitution(other);
		if (!U.HasUniqueSolution())
			throw MatrixEquationNoUniqueSolutionException();
		x = U.backSubstitution(y);
	}

MEASURE_TIME_END
	return x;
}

template <typename T>
T &Matrix<T>::operator[](Dimensions indecies)
{
	return rawData(indecies.n,indecies.m);
}

template class Matrix<int>;
template class Matrix<double>;
