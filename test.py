from py_matrix import Matrix, MatrixSizeDisparityException
import copy

matrix = Matrix("1,2,3;4,5,6;7,8,9")
matrix.print()
transposed_matrix = matrix.transpose()
transposed_matrix.print()

multi = transposed_matrix * matrix
multi.print()

test1 = Matrix(1000,1000,1)
test2 = Matrix(1000,1000,1)
check = Matrix(1000,1000,1000)

result = test1 * test2

try:
    xd = result * multi
except MatrixSizeDisparityException:
    print("bad sizes")

print(result[(1,1)])
print(test1 == test2)
print(result == check)

A = Matrix("2,1;3,-2")
b = Matrix("5;7")

x = A | b
x.print()

c = copy.deepcopy(x)
c[(0,0)] = 0

c.print()
x.print()

empty = Matrix()
empty.print()

size = c.getSize()
print(size)