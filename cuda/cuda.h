#include <cuda_runtime.h>

template<typename T>
__global__ void cudaKernel(T* data, int size);