import py_matrix as mat

matrix = mat.Matrix("1,2,3;4,5,6;7,8,9")
matrix.print()
transposed_matrix = matrix.transpose()
transposed_matrix.print()

multi = transposed_matrix * matrix
multi.print()

test1 = mat.Matrix(2000,1000,1)
test2 = mat.Matrix(1000,2000,1)
check = mat.Matrix(2000,2000,1000)

result = test1 * test2

try:
    xd = result * multi
except mat.MatrixSizeDisparityException:
    print("bad sizes")

print(result[mat.Dimensions(1,1)])
print(test1 == test2)
print(result == check)

A = mat.Matrix("2,1;3,-2")
b = mat.Matrix("5;7")

x = A | b
x.print()