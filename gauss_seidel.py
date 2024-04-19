from py_matrix import Matrix
from common import norm
import time

def gauss_seidel_iter_method(A: Matrix, b: Matrix, max_iterations: int):
    start_time = time.time()

    D = A.diag(0)
    L = A.tril(-1)
    U = A.triu(1)

    DL = D+L

    M = -DL
    bm = DL | b

    x = Matrix(A.getSize()[0], 1, 1)

    iterations = 0
    err_norms = []
    for _ in range(max_iterations):
        iterations += 1
        x = (M | (U*x)) + bm

        err_norm = norm(A*x-b)
        err_norms.append(err_norm)
        if err_norm < 1e-9:
            break


    end_time = time.time()
    time_taken = end_time - start_time

    return x, err_norms, iterations, time_taken