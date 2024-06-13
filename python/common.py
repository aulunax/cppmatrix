from py_matrix import Matrix
import math as mmm
def norm(A: Matrix): 
    norm = 0
    for index in range(A.getSize()[0]):
        norm += A[(index,0)]**2
    return mmm.sqrt(norm)

def create_matrix(N, a1, a2, a3):
    mat = Matrix(N,N,a1)
    mat = mat.diag(0)
    diag2 = Matrix(N,N,a2).diag(-1) + Matrix(N,N,a2).diag(1)
    diag3 = Matrix(N,N,a3).diag(-2) + Matrix(N,N,a3).diag(2)
    mat = mat + diag2 + diag3

    return mat
