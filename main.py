from py_matrix import Matrix
from python.direct import direct_method
from python.jacobi import jacobi_iter_method
from python.gauss_seidel import gauss_seidel_iter_method
from python.common import create_matrix
from math import sin


for num in [100, 200, 500]:
    print("\nMatrix size: " + str(num) + "\n")
    N = num # 948
    a1, a2, a3 = 7, -1, -1

    A = create_matrix(N,a1,a2,a3)

    b = Matrix(N,1,1)
    for i in range(N):
        b[(i,0)] = sin((i+1)*4)

    x, err_norms, iterations, time_taken = jacobi_iter_method(A, b, 100)
    print("Jacobi method:")
    print("iterations: " + str(iterations))
    print("time_taken: " + str(time_taken))
    print("err_norm: " + str(err_norms))

    x, err_norms, iterations, time_taken = gauss_seidel_iter_method(A, b, 100)
    print("Gauss-Seidel method:")
    print("iterations: " + str(iterations))
    print("time_taken: " + str(time_taken))
    print("err_norm: " + str(err_norms))

    x, time_taken = direct_method(A, b)
    print("Direct method:")
    print("time_taken: " + str(time_taken))
