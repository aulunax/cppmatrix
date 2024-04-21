#pragma once
#include "MatrixStructs.h"
#include <vector>
#include <string>
#include <stdexcept>
#include <functional>

#ifdef MATRIX_DEBUG
#define Duration(a) std::chrono::duration_cast<std::chrono::microseconds>(a)
#endif




class MatrixSizeDisparityException : public std::exception {
public:
    virtual const char* what() const throw() {
        return "Matrices aren't of the same size";
    }
};

class NotSquareMatrixException : public std::exception {
public:
    virtual const char* what() const throw() {
        return "Matrix needs to be square";
    }
};

class MatrixEquationNoUniqueSolutionException : public std::exception {
public:
	virtual const char* what() const throw() {
		return "Matrix equation doesn't have unique solution";
	}
};

class MatrixNotInvertible : public std::exception {
public:
	virtual const char* what() const throw() {
		return "Matrix is not invertible";
	}
};


template<typename T>
class Matrix
{
	typedef std::vector<std::vector<T>> Data;
	typedef std::function<void(const Data&, const Data&, const Data&, Data&, int, int)> ThreadedFunction;

	Dimensions size;
	Data rawData;

#ifdef MATRIX_DEBUG
	template<typename Func, typename... Args>
	auto returnTypeFunctionWrapper(const std::string& functionName, Func&& func, Args&&... args) const {
		auto startTime = std::chrono::high_resolution_clock::now();
		auto result = std::invoke(std::forward<Func>(func), std::forward<Args>(args)...);
		auto endTime = std::chrono::high_resolution_clock::now();
		auto duration = Duration(endTime - startTime);
		std::cout << "Function: " << functionName << "\nTime taken: " << duration.count()/1000.0 << " milliseconds" << std::endl;
		return result;
	}
	#else
	template<typename Func, typename... Args>
	auto returnTypeFunctionWrapper(const std::string& functionName, Func&& func, Args&&... args) const {
		return std::invoke(std::forward<Func>(func), std::forward<Args>(args)...);
	}
#endif


	// unused
    static void fillWithValue(const Data& args, const Data &data1, const Data &data2, Data &result, int startRow, int endRow);

	// functions that can be threaded using threadedMatrixOperation(...)
	static void addMatrices(const Data& args, const Data& matrix1, const Data& matrix2, Data& result, int startRow, int endRow);
	static void substractMatrices(const Data& args, const Data& data1, const Data& data2, Data& result, int startRow, int endRow);
	static void multiplyMatrices(const Data& args, const Data& data1, const Data& data2, Data& result, int startRow, int endRow);
    static void multiplyByConstant(const Data& args, const Data &data1, const Data &data2, Data &result, int startRow, int endRow);
    static void divideByConstant(const Data& args, const Data &data1, const Data &data2, Data &result, int startRow, int endRow);
	// args contains k value of diag function
    static void getDiagonal(const Data& args, const Data &data1, const Data &data2, Data &result, int startRow, int endRow);
	// args[0][0] contains k value of diag function
	// args[0][1] contains direction of triangle
	// direction: -1 means down, 1 means up 
    static void getTriangle(const Data& args, const Data &data1, const Data &data2, Data &result, int startRow, int endRow);
	static void getTransposed(const Data& args, const Data &data1, const Data &data2, Data &result, int startRow, int endRow);

	// Splits the matrix into even chunks and performs given 'operation' on each chunk using multithreading
	// this - first operand
	// other - second operand
	// The multithreaded function has to have following properties:
	// - performs 'operation' on all rows from startRow to endRow
	// - all cells of the result matrix must be filled with some value
	// - data1 contains rawData of 'this' matrix 
	// - data2 contains rawData of 'other' matrix for operations with 2 operands
	// - args contains additional arguments
    Matrix threadedMatrixOperation(const Matrix& args, const Matrix& other, ThreadedFunction operation, Dimensions resultSize = {0,0}) const;

	// https://en.wikipedia.org/wiki/Gaussian_elimination#Pseudocode
	// returns copy of concatenated matrix A and b in row echelon form
	// unused (replaced by LU factorization)
	static Matrix gaussianElimination(const Matrix& A, const Matrix& b);

	// checks if rank of matrix is equal to number of its columns
	bool HasUniqueSolution() const;

	bool isTril() const;
	bool isTriu() const;

	// returns solution to (this*)*b = x
	Matrix backSubstitution(const Matrix &b) const;
	Matrix forwardSubstitution(const Matrix &b) const;
public:
	Matrix() : size({0,0}) {};
	Matrix(const std::string& content);
	Matrix(int n, int m);
	Matrix(int n, int m, T value);
	Matrix(Dimensions size);
	Matrix(Dimensions size, T value);

	Matrix(const Matrix& other);
	Matrix(Matrix&& other);

	Dimensions getSize() const {return size;};
	void print();
	void reserve(int n, int m);
	void fill(int n, int m, T value);

	// only implemented for diagonal matrices
	Matrix inv() const;

	// k means which diagonal should resulting matrix contain (0 means main diagonal)
	Matrix diag(int k=0) const;
	// k means from which diagonal should resulting matrix start (0 means main diagonal)
	Matrix tril(int k=0) const;
	Matrix triu(int k=0) const;

	Matrix transpose() const;

	Matrix& operator=(const Matrix& other);
	Matrix& operator=(Matrix&& other) noexcept;

	Matrix operator+(const Matrix& other) const;
	Matrix operator-(const Matrix& other) const;
	Matrix operator*(const Matrix& other) const;
	Matrix operator*(const T& value) const;
	Matrix operator/(const T& value) const;
	Matrix operator-() const;

	// solves x*(*this) = other
	Matrix operator/(const Matrix& other) const;
	// solves (*this)*x = other
	// handles 3 cases:
	// upper triangular, lower triangular and miscellaneous matrix
	Matrix operator|(const Matrix& other) const;

	T& operator[](Dimensions indecies);

	friend bool operator==(const Matrix<T>& a, const Matrix<T>& b) {
		if (a.size != b.size)
		return false;
		for (int i = 0; i < a.size.n; i++) {
			for (int j = 0; j < a.size.m; j++) {
				if (a.rawData[i][j] != b.rawData[i][j])
					return false;
		}
	}
	return true;
	}
	friend bool operator!=(const Matrix<T>& a, const Matrix<T>& b) {
		return !(a == b);
	}

};

