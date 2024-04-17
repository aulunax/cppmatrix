#include "Matrix.h"
#include "MatrixMisc.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <thread>
#include <future>

#ifdef MATRIX_DEBUG
	#include <chrono>

	#define Duration(a) std::chrono::duration_cast<std::chrono::microseconds>(a)
#endif

template <typename T>
Matrix<T>::Matrix(const std::string &content) : Matrix()
{
	#ifdef MATRIX_DEBUG
		auto startTime = std::chrono::high_resolution_clock::now();
	#endif
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
				rawData[i][j] = value;
				j++;
			}
		}));
		i++;
	}

	for (auto& future : futures) {
        future.wait();
    }

	#ifdef MATRIX_DEBUG
		auto endTime = std::chrono::high_resolution_clock::now();
    	auto duration = Duration(endTime - startTime);
		std::cout << "Matrix(const std::string &content) finished in " << duration.count()/1000.0 << " miliseconds.\n\n";
	#endif
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
		int i = 0;
		for (const auto& row : other.rawData) {
       		rawData[i] = row;
			i++;
    	}
    }
}


template<typename T>
Matrix<T>::Matrix(Matrix&& other)
{
	size = other.size;
    if (this != &other) {
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


	for (const auto& row : rawData) {
		for (const auto& elem : row) {
			std::cout << elem << ' ';
		}
		std::cout << '\n';
	}
	std::cout << '\n';
}

template <typename T>
void Matrix<T>::reserve(int n, int m)
{
	size = {n,m};
	rawData = Data(n, std::vector<T>(m));
}

template <typename T>
void Matrix<T>::fill(int n, int m, T value)
{
#ifdef MATRIX_DEBUG
	//std::cout << "Starting Matrix<T>::fill(int n, int m, T value)\n";
	auto startTime = std::chrono::high_resolution_clock::now();
#endif

	size = {n,m};
	rawData = Data(n, std::vector<T>(m, value));
	
#ifdef MATRIX_DEBUG
	auto endTime = std::chrono::high_resolution_clock::now();
	auto duration = Duration(endTime - startTime);
	//std::cout << "Matrix<T>::fill(int n, int m, T value) finished in " << duration.count()/1000.0 << " miliseconds.\n\n";
#endif
}


template <typename T>
Matrix<T> Matrix<T>::diag(int k) const
{
#ifdef MATRIX_DEBUG
	std::cout << "Starting Matrix<T>::diag(int k) \n";
#endif
	if (size.m != size.n) {
		throw NotSquareMatrixException();
	}

	Matrix<T> args = Matrix<T>(1, 1, k);
	Matrix<T> other;

	return threadedMatrixOperation(args, other, &getDiagonal);
}

template <typename T>
void Matrix<T>::getDiagonal(const Data& args, const Data &data1, const Data &data2, Data &result, int startRow, int endRow)
{
	int diagonal = args[0][0];
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data1[0].size(); j++) {
			if (i+diagonal == j) 
				result[i][j] = data1[i][j];
			else
				result[i][j] = 0;
		}
	}
}

template <typename T>
Matrix<T> Matrix<T>::tril(int k) const
{
#ifdef MATRIX_DEBUG
	std::cout << "Starting Matrix<T>::tril(int k) \n";
#endif
	if (size.m != size.n) {
		throw NotSquareMatrixException();
	}

    Matrix<T> args = Matrix<T>(1,2);
	args[{0,0}] = k;
	args[{0,1}] = -1;

	Matrix<T> other;

	return threadedMatrixOperation(args, other, &getTriangle);
}

template <typename T>
Matrix<T> Matrix<T>::triu(int k) const
{
#ifdef MATRIX_DEBUG
	std::cout << "Starting Matrix<T>::triu(int k) \n";
#endif
	if (size.m != size.n) {
		throw NotSquareMatrixException();
	}

    Matrix<T> args = Matrix<T>(1,2);
	args[{0,0}] = k;
	args[{0,1}] = 1;

	Matrix<T> other;

	return threadedMatrixOperation(args, other, &getTriangle);
}

template <typename T>
void Matrix<T>::getTransposed(const Data &args, const Data &data1, const Data &data2, Data &result, int startRow, int endRow)
{
	for (int i = startRow; i < endRow; ++i) {
        for (int j = 0; j < data1[0].size(); ++j) {
            result[j][i] = data1[i][j];
        }
    }
}

template <typename T>
Matrix<T> Matrix<T>::transpose() const
{
#ifdef MATRIX_DEBUG
	std::cout << "Starting Matrix<T>::transpose() \n";
	//auto startTime = std::chrono::high_resolution_clock::now();
#endif
	Matrix<T> args;
	Matrix<T> other;

	return threadedMatrixOperation(args, other, &getTransposed, Dimensions{size.m, size.n});
}

template <typename T>
void Matrix<T>::getTriangle(const Data& args, const Data &data1, const Data &data2, Data &result, int startRow, int endRow)
{
	int direction = args[0][1];
	int k = args[0][0];
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data1[0].size(); j++) {
			if (direction == 1 ? i+k <= j : i+k>=j) 
				result[i][j] = data1[i][j];
			else
				result[i][j] = 0;
		}
	}
}



template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix& other) {
	reserve(other.size.n, other.size.m);
    if (this != &other) {
		int i = 0;
		for (const auto& row : other.rawData) {
       		rawData[i] = row;
			i++;
    	}
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(Matrix&& other) noexcept {
	size = other.size;
    if (this != &other) {
        rawData = std::move(other.rawData);
    }
    return *this;
}

template <typename T>
void Matrix<T>::fillWithValue(const Data& args, const Data &data1, const Data &data2, Data &result, int startRow, int endRow)
{
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < result[0].size(); j++) {
			result[i][j] = data2[0][0];
		}
	}
}

template <typename T>
void Matrix<T>::addMatrices(const Data& args, const Data &data1, const Data &data2, Data &result, int startRow, int endRow)
{
    for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data1[0].size(); j++) {
			result[i][j] = data1[i][j] + data2[i][j];
		}
	}
}

template<typename T>
void Matrix<T>::substractMatrices(const Data& args, const Data& data1, const Data& data2, Data& result, int startRow, int endRow) {
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data1[0].size(); j++) {
			result[i][j] = data1[i][j] - data2[i][j];
		}
	}
}

template<typename T>
void Matrix<T>::multiplyMatrices(const Data& args, const Data& data1, const Data& data2, Data& result, int startRow, int endRow) {
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data2[0].size(); j++) {
			for (int k = 0; k < data1[0].size(); k++) {
				result[i][j] += data1[i][k] * data2[k][j];
			}
		}
	}
}

template<typename T>
void Matrix<T>::multiplyByConstant(const Data& args, const Data& data1, const Data& data2, Data& result, int startRow, int endRow) {
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data1[0].size(); j++) {
			result[i][j] = data1[i][j] * args[0][0];
		}
	}
}

template<typename T>
void Matrix<T>::divideByConstant(const Data& args, const Data& data1, const Data& data2, Data& result, int startRow, int endRow) {
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data1[0].size(); j++) {
			result[i][j] = data1[i][j] / args[0][0];
		}
	}
}

template<typename T>
Matrix<T> Matrix<T>::threadedMatrixOperation(const Matrix& args, const Matrix& other, ThreadedFunction operation, Dimensions resultSize) const {

#ifdef MATRIX_DEBUG
	auto startTime = std::chrono::high_resolution_clock::now();
#endif

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
	if (size.n > numOfThreads*2) {

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
#ifdef MATRIX_DEBUG
	auto endTime = std::chrono::high_resolution_clock::now();
	auto duration = Duration(endTime - startTime);
	std::cout << "threadedMatrixOperation finished in " << duration.count()/1000.0 << " miliseconds.\n\n";
#endif
	return temp;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix& other) const
{
#ifdef MATRIX_DEBUG
	std::cout << "Starting Matrix<T>::operator+(const Matrix& other)\n";
#endif
	if (size != other.size) {
		throw MatrixSizeDisparityException();
	}
	Matrix<T> args;

	return threadedMatrixOperation(args, other, &addMatrices);
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix& other) const
{
#ifdef MATRIX_DEBUG
	std::cout << "Starting Matrix<T>::operator-(const Matrix& other)\n";
#endif
	if (size != other.size) {
		throw MatrixSizeDisparityException();
	}
	Matrix<T> args;

	return threadedMatrixOperation(args, other, &substractMatrices);
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix& other) const
{
#ifdef MATRIX_DEBUG
	std::cout << "Starting Matrix<T>::operator*(const Matrix& other)\n";
#endif
	if (size.m != other.size.n) {
		throw MatrixSizeDisparityException();
	}
	Matrix<T> args;

	return threadedMatrixOperation(args, other, &multiplyMatrices);
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T& value) const
{
#ifdef MATRIX_DEBUG
	std::cout << "Starting Matrix<T>::operator*(const T& value)\n";
#endif
	Matrix<T> args = Matrix<T>(1, 1, value);
	Matrix<T> other;

	return threadedMatrixOperation(args, other, &multiplyByConstant);
}

template<typename T>
Matrix<T> Matrix<T>::operator/(const T& value) const
{
#ifdef MATRIX_DEBUG
	std::cout << "Starting Matrix<T>::operator/(const T& value)\n";
#endif
	if (value == 0) {
		throw std::runtime_error("Division by zero error");
	}
	Matrix<T> args = Matrix<T>(1, 1, value);
	Matrix<T> other;

	return threadedMatrixOperation(args, other, &divideByConstant);
}

template <typename T>
Matrix<T> Matrix<T>::operator-() const
{
    return (*this)*-1;
}

template <typename T>
Matrix<T> Matrix<T>::gaussianElimination(const Matrix &A, const Matrix &b)
{
	if (A.size.n != b.size.n || b.size.m != 1) {
		throw MatrixSizeDisparityException();
	}

	Matrix<T> Ab(A.size.n, A.size.m+1);
	for (int i = 0; i < A.size.n; i++) {
		for (int j = 0; j < A.size.m; j++) {
			Ab.rawData[i][j] = A.rawData[i][j];	
		}
		Ab.rawData[i][A.size.m] = b.rawData[i][0];	
	}

#ifdef MATRIX_DEBUG
	std::cout << "Concat result:\n";
	Ab.print();
#endif

	int h = 0;
	int k = 0;

	while (h < Ab.size.n && k < Ab.size.m) {
		int valMax = abs(Ab[{h,k}]);
		int iMax = h;
		for (int i = h; i < Ab.size.n; i++) {
			if (abs(Ab[{i, k}]) > iMax) {
				valMax = abs(Ab[{i, k}]);
				iMax = i;
			}
		}

		if (Ab[{iMax, k}] == 0)
        	k += 1;
		else {
			std::swap(Ab.rawData[h], Ab.rawData[iMax]);
			for (int i = h+1; i < Ab.size.n; i++) {
				T f = Ab[{i, k}] / Ab[{h, k}];
				Ab[{i, k}] = 0;
				for (int j = k+1; j < Ab.size.m; j++) {
					Ab[{i, j}] = Ab[{i, j}] - Ab[{h, j}] * f;
				}
			}
			h += 1;
			k += 1;
		}
#ifdef MATRIX_DEBUG
	std::cout << "Gaussian elimination step:\n";
	Ab.print();
#endif
	}
    return Ab;
}

template <typename T>
Matrix<T> Matrix<T>::operator/(const Matrix &other) const
{
#ifdef MATRIX_DEBUG
	std::cout << "Starting Matrix<T>::operator/(const Matrix &other)\n\n";
#endif

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
	std::cout << "Starting Matrix<T>::operator|(const Matrix &other)\n\n";
	auto startTime = std::chrono::high_resolution_clock::now();
#endif

	Matrix<T> gaussianResult = gaussianElimination((*this), other);

#ifdef MATRIX_DEBUG
	std::cout << "Gaussian elimination result:\n\n";
	gaussianResult.print();
#endif

	// check if there is unique solution by calculating the rank of the matrix
	// https://en.wikipedia.org/wiki/Rank_(linear_algebra)
	int ranks = 0;
	for (int i = 0; i < size.n; i++) {
		bool zeros = true;
		for (int j = 0; j < gaussianResult.size.m - 1; j++) {
			if (gaussianResult[{i, j}] != 0) {
				zeros = false;
				break;
			}
		}
		if (!zeros)
			ranks += 1;
	}
	if (ranks != gaussianResult.size.m - 1) {
		throw MatrixEquationNoUniqueSolutionException();
	}

	Matrix<T> x(size.n, 1);

	// get result by back substitution
	// https://algowiki-project.org/en/Backward_substitution
	for (int i = size.n - 1; i >= 0; i--) {
		x[{i, 0}] = gaussianResult[{i, size.n}];
		for (int j = i + 1; j < size.n; ++j) {
			x[{i, 0}] -= gaussianResult[{i, j}] * x[{j, 0}];
		}
		x[{i, 0}] /= gaussianResult[{i, i}];
	}

#ifdef MATRIX_DEBUG
	auto endTime = std::chrono::high_resolution_clock::now();
	auto duration = Duration(endTime - startTime);
	std::cout << "Matrix<T>::operator|(const Matrix &other) finished in " << duration.count()/1000.0 << " miliseconds.\n\n";
#endif
	return x;
}

template <typename T>
T &Matrix<T>::operator[](Dimensions indecies)
{
	return rawData[indecies.n][indecies.m];
}

template <typename T>
bool Matrix<T>::operator==(const Matrix &other)
{
	if (size != other.size)
		return false;
    for (int i = 0; i < size.n; i++) {
		for (int j = 0; j < size.m; j++) {
			if (rawData[i][j] != other.rawData[i][j])
				return false;
		}
	}
	return true;
}

template <typename T>
bool Matrix<T>::operator!=(const Matrix &other)
{
    return !(*this == other);
}

template class Matrix<int>;
template class Matrix<float>;
template class Matrix<double>;
