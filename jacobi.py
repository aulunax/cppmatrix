from py_matrix import Matrix
import time
from common import norm

def jacobi_iter_method(A: Matrix, b: Matrix, max_iterations: int):
    D = A.diag(0)
    L = A.tril(-1)
    U = A.triu(1)

    D_inv = D.inv()

    M = -D_inv*(L+U)
    bm = D_inv*b

    x = Matrix(A.getSize()[0], 1, 1)

    start_time = time.time()
    iterations = 0
    err_norms = []
    for _ in range(max_iterations):
        iterations += 1
        x = M*x+bm

        err_norm = norm(A*x-b)
        err_norms.append(err_norm)
        if err_norm < 1e-9:
            break


    end_time = time.time()
    time_taken = end_time - start_time

    return x, err_norms, iterations, time_taken