#pragma once
#include "MatrixStructs.h"
#include <vector>
#include <string>
#include <stdexcept>
#include <functional>


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

template<typename T>
class Matrix
{
	typedef std::vector<std::vector<T>> Data;
	typedef std::function<void(const Data&, const Data&, const Data&, Data&, int, int)> ThreadedFunction;

	Dimensions size;
	Data rawData;

	// unused
    static void fillWithValue(const Data& args, const Data &data1, const Data &data2, Data &result, int startRow, int endRow);

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
	static Matrix gaussianElimination(const Matrix& A, const Matrix& b);

public:
	Matrix() : size({0,0}) {};
	Matrix(const std::string& content);
	Matrix(int n, int m);
	Matrix(int n, int m, T value);
	Matrix(Dimensions size);
	Matrix(Dimensions size, T value);

	Matrix(const Matrix& other);
	Matrix(Matrix&& other);

	void print();
	void reserve(int n, int m);
	void fill(int n, int m, T value);

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
	Matrix operator|(const Matrix& other) const;


	T& operator[](Dimensions indecies);
	bool operator==(const Matrix& other);
	bool operator!=(const Matrix& other);
};

