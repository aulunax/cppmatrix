from py_matrix import Matrix
import time

def direct_method(A: Matrix, b: Matrix):
    start_time = time.time()

    x = A | b

    end_time = time.time()
    time_taken = end_time - start_time

    return x, time_taken