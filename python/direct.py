from py_matrix import Matrix
import time

from python.common import norm

def direct_method(A: Matrix, b: Matrix):
    start_time = time.time()

    x = A | b

    end_time = time.time()
    time_taken = end_time - start_time

    err_norms = norm(A*x-b)

    return x, err_norms, time_taken