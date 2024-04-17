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

template<typename T>
class Matrix
{
	Dimensions size;
	std::vector<std::vector<T>> rawData;

	// unused
    static void fillWithValue(const std::vector<std::vector<T>> &data1, const std::vector<std::vector<T>> &data2, std::vector<std::vector<T>> &result, int startRow, int endRow);

	static void addMatrices(const std::vector<std::vector<T>>& matrix1, const std::vector<std::vector<T>>& matrix2, std::vector<std::vector<T>>& result, int startRow, int endRow);
	static void substractMatrices(const std::vector<std::vector<T>>& data1, const std::vector<std::vector<T>>& data2, std::vector<std::vector<T>>& result, int startRow, int endRow);
	static void multiplyMatrices(const std::vector<std::vector<T>>& data1, const std::vector<std::vector<T>>& data2, std::vector<std::vector<T>>& result, int startRow, int endRow);
    static void multiplyByConstant(const std::vector<std::vector<T>> &data1, const std::vector<std::vector<T>> &data2, std::vector<std::vector<T>> &result, int startRow, int endRow);
    static void divideByConstant(const std::vector<std::vector<T>> &data1, const std::vector<std::vector<T>> &data2, std::vector<std::vector<T>> &result, int startRow, int endRow);

	// data2 contains k value of diag function
    static void getDiagonal(const std::vector<std::vector<T>> &data1, const std::vector<std::vector<T>> &data2, std::vector<std::vector<T>> &result, int startRow, int endRow);
	// data2[0][0] contains k value of diag function
	// data2[0][1] contains direction of triangle
	// direction: -1 means down, 1 means up 
    static void getTriangle(const std::vector<std::vector<T>> &data1, const std::vector<std::vector<T>> &data2, std::vector<std::vector<T>> &result, int startRow, int endRow);

	// Splits the matrix into even chunks and performs given 'operation' on each chunk using multithreading
	// this - first operand
	// other - second operand
	// The multithreaded function has to have following properties:
	// - performs 'operation' on all rows from startRow to endRow
	// - all cells of the result matrix must be filled with some value
	// - data1 contains rawData of 'this' matrix 
	// - data2 contains rawData of 'other' matrix for operations with 2 operands
	// - data2 contains additional arguments for operations with 1 operand
    Matrix threadedMatrixOperation(const Matrix& other, std::function<void(const std::vector<std::vector<T>>&, const std::vector<std::vector<T>>&, std::vector<std::vector<T>>&, int, int)> operation) const;

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


	Matrix& operator=(const Matrix& other);
	Matrix& operator=(Matrix&& other) noexcept;

	Matrix operator+(const Matrix& other) const;
	Matrix operator-(const Matrix& other) const;
	Matrix operator*(const Matrix& other) const;
	Matrix operator*(const T& value) const;
	Matrix operator/(const T& value) const;

	T& operator[](Dimensions indecies);

	


	
};

