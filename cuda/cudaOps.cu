#include "cudaOps.h"
#include <device_launch_parameters.h>


__global__ void matrixMul(double* a, double* b, double* c, int rowsA, int colsA, int colsB) {
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (row < rowsA && col < colsB) {
        double sum = 0.0;
        for (int k = 0; k < colsA; ++k) {
            sum += a[row * colsA + k] * b[k * colsB + col];
        }
        c[row * colsB + col] = sum;
    }
}

// Host function to perform matrix multiplication
void mulMatrix(const Matrix<double>& a, const Matrix<double>& b, Matrix<double>& c) {
    int rowsA = a.getSize().n;
    int colsA = a.getSize().m;
    int colsB = b.getSize().m;

    // Allocate device memory
    double* dev_a, *dev_b, *dev_c;
    cudaMalloc((void**)&dev_a, rowsA * colsA * sizeof(double));
    cudaMalloc((void**)&dev_b, colsA * colsB * sizeof(double));
    cudaMalloc((void**)&dev_c, rowsA * colsB * sizeof(double));

    // Copy input matrices from host to device
    cudaMemcpy(dev_a, a.getDataVec().data(), rowsA * colsA * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, b.getDataVec().data(), colsA * colsB * sizeof(double), cudaMemcpyHostToDevice);

    // Define grid and block dimensions
    dim3 blockDim(16, 16);
    dim3 gridDim((colsB + blockDim.x - 1) / blockDim.x, (rowsA + blockDim.y - 1) / blockDim.y);

    // Launch kernel
    matrixMul<<<gridDim, blockDim>>>(dev_a, dev_b, dev_c, rowsA, colsA, colsB);

    // Copy result matrix from device to host
    cudaMemcpy(c.getDataVec().data(), dev_c, rowsA * colsB * sizeof(double), cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_c);
}