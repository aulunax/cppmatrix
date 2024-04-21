#include <cuda_runtime.h>
#include <vector>
#include "Matrix.h"


void mulMatrix(const Matrix<double>& a, const Matrix<double>& b, Matrix<double>& c);

__global__ void matrixMul(double* a, double* b, double* c, int rowsA, int colsA, int colsB);
