import copy
from py_matrix import Matrix, MatrixSizeDisparityException
from jacobi import jacobi_iter_method
from gauss_seidel import gauss_seidel_iter_method
from common import create_matrix, norm
from math import sin

N = 3000 # 948
a1, a2, a3 = 7, -1, -1

A = create_matrix(N,a1,a2,a3)

b = Matrix(N,1,1)
for i in range(N):
    b[(i,0)] = sin((i+1)*4)

# x, err_norms, iterations, time_taken = jacobi_iter_method(A, b, 100)
# print("iterations: " + str(iterations))
# print("time_taken: " + str(time_taken))
# print("err_norm: " + str(err_norms))

x, err_norms, iterations, time_taken = gauss_seidel_iter_method(A, b, 100)
print("iterations: " + str(iterations))
print("time_taken: " + str(time_taken))
print("err_norm: " + str(err_norms))
