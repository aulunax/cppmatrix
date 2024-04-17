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
			rawData = std::vector<std::vector<T>>(); 
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
	rawData = std::vector<std::vector<T>>(n, std::vector<T>(m));
}

template <typename T>
void Matrix<T>::fill(int n, int m, T value)
{
#ifdef MATRIX_DEBUG
	std::cout << "Starting Matrix<T>::fill(int n, int m, T value)\n";
	auto startTime = std::chrono::high_resolution_clock::now();
#endif

	size = {n,m};
	rawData = std::vector<std::vector<T>>(n, std::vector<T>(m, value));
	
#ifdef MATRIX_DEBUG
	auto endTime = std::chrono::high_resolution_clock::now();
	auto duration = Duration(endTime - startTime);
	std::cout << "Matrix<T>::fill(int n, int m, T value) finished in " << duration.count()/1000.0 << " miliseconds.\n\n";
#endif
}


template <typename T>
Matrix<T> Matrix<T>::diag(int k) const
{
#ifdef MATRIX_DEBUG
	std::cout << "Starting Matrix<T>::diag(int k) \n";
#endif
	Matrix<T> args = Matrix<T>(1, 1, k);

	return threadedMatrixOperation(args, &getDiagonal);
}

template <typename T>
void Matrix<T>::getDiagonal(const std::vector<std::vector<T>> &data1, const std::vector<std::vector<T>> &data2, std::vector<std::vector<T>> &result, int startRow, int endRow)
{
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data1[0].size(); j++) {
			if (i+data2[0][0] == j) 
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
    Matrix<T> args = Matrix<T>(1,2);
	args[{0,0}] = k;
	args[{0,1}] = -1;

	return threadedMatrixOperation(args, &getTriangle);
}

template <typename T>
Matrix<T> Matrix<T>::triu(int k) const
{
#ifdef MATRIX_DEBUG
	std::cout << "Starting Matrix<T>::triu(int k) \n";
#endif
    Matrix<T> args = Matrix<T>(1,2);
	args[{0,0}] = k;
	args[{0,1}] = 1;

	return threadedMatrixOperation(args, &getTriangle);
}

template <typename T>
void Matrix<T>::getTriangle(const std::vector<std::vector<T>> &data1, const std::vector<std::vector<T>> &data2, std::vector<std::vector<T>> &result, int startRow, int endRow)
{
	int direction = data2[0][1];
	int k = data2[0][0];
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
void Matrix<T>::fillWithValue(const std::vector<std::vector<T>> &data1, const std::vector<std::vector<T>> &data2, std::vector<std::vector<T>> &result, int startRow, int endRow)
{
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < result[0].size(); j++) {
			result[i][j] = data2[0][0];
		}
	}
}

template <typename T>
void Matrix<T>::addMatrices(const std::vector<std::vector<T>> &data1, const std::vector<std::vector<T>> &data2, std::vector<std::vector<T>> &result, int startRow, int endRow)
{
    for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data1[0].size(); j++) {
			result[i][j] = data1[i][j] + data2[i][j];
		}
	}
}

template<typename T>
void Matrix<T>::substractMatrices(const std::vector<std::vector<T>>& data1, const std::vector<std::vector<T>>& data2, std::vector<std::vector<T>>& result, int startRow, int endRow) {
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data1[0].size(); j++) {
			result[i][j] = data1[i][j] - data2[i][j];
		}
	}
}

template<typename T>
void Matrix<T>::multiplyMatrices(const std::vector<std::vector<T>>& data1, const std::vector<std::vector<T>>& data2, std::vector<std::vector<T>>& result, int startRow, int endRow) {
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data2[0].size(); j++) {
			for (int k = 0; k < data1[0].size(); k++) {
				result[i][j] += data1[i][k] * data2[k][j];
			}
		}
	}
}

template<typename T>
void Matrix<T>::multiplyByConstant(const std::vector<std::vector<T>>& data1, const std::vector<std::vector<T>>& data2, std::vector<std::vector<T>>& result, int startRow, int endRow) {
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data1[0].size(); j++) {
			// assume data2 is (1,1) matrix
			result[i][j] = data1[i][j] * data2[0][0];
		}
	}
}

template<typename T>
void Matrix<T>::divideByConstant(const std::vector<std::vector<T>>& data1, const std::vector<std::vector<T>>& data2, std::vector<std::vector<T>>& result, int startRow, int endRow) {
	for (int i = startRow; i < endRow; i++) {
		for (int j = 0; j < data1[0].size(); j++) {
			// assume data2 is (1,1) matrix
			result[i][j] = data1[i][j] / data2[0][0];
		}
	}
}

template<typename T>
Matrix<T> Matrix<T>::threadedMatrixOperation(const Matrix& other, std::function<void(const std::vector<std::vector<T>>&, const std::vector<std::vector<T>>&, std::vector<std::vector<T>>&, int, int)> operation) const {

#ifdef MATRIX_DEBUG
	auto startTime = std::chrono::high_resolution_clock::now();
#endif
	int numOfThreads = std::min(size.n, (int)std::thread::hardware_concurrency());
	int rowsPerThread = size.n / numOfThreads;

	Matrix<T> temp = Matrix<T>(size);

	std::vector<std::thread> threads;
	for (int i = 0; i < numOfThreads - 1; ++i) {
		threads.emplace_back(operation, std::ref(rawData), std::ref(other.rawData), std::ref(temp.rawData), i * rowsPerThread, (i + 1) * rowsPerThread);
	}
	threads.emplace_back(operation, std::ref(rawData), std::ref(other.rawData), std::ref(temp.rawData), (numOfThreads - 1) * rowsPerThread, size.n);


	for (auto& thread : threads) {
		thread.join();
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

	return threadedMatrixOperation(other, &addMatrices);
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

	return threadedMatrixOperation(other, &substractMatrices);
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

	return threadedMatrixOperation(other, &multiplyMatrices);
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T& value) const
{
#ifdef MATRIX_DEBUG
	std::cout << "Starting Matrix<T>::operator*(const T& value)\n";
#endif
	Matrix<T> other = Matrix<T>(1, 1, value);

	return threadedMatrixOperation(other, &multiplyByConstant);
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
	Matrix<T> other = Matrix<T>(1, 1, value);

	return threadedMatrixOperation(other, &divideByConstant);
}

template <typename T>
T &Matrix<T>::operator[](Dimensions indecies)
{
	return rawData[indecies.n][indecies.m];
}

template class Matrix<int>;
template class Matrix<float>;
template class Matrix<double>;
