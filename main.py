import os
from py_matrix import Matrix
from python.direct import direct_method
from python.jacobi import jacobi_iter_method
from python.gauss_seidel import gauss_seidel_iter_method
from python.common import create_matrix
from math import sin
import matplotlib.pyplot as plt
import numpy as np

# ZADANIE A,B
a1, a2, a3 = 7, -1, -1
N = 948
A = create_matrix(N,a1,a2,a3)
b = Matrix(N,1,1)
for i in range(N):
    b[(i,0)] = sin((i+1)*4)

x, jacobi_err_norms, jacobi_iterations, time_taken = jacobi_iter_method(A, b, 100)
print("Jacobi method:")
print("iterations: " + str(jacobi_iterations))
print("time_taken: " + str(time_taken))
print("err_norm: " + str(jacobi_err_norms))

x, gauss_seidel_err_norms, gauss_seidel_iterations, time_taken = gauss_seidel_iter_method(A, b, 100)
print("Gauss-Seidel method:")
print("iterations: " + str(gauss_seidel_iterations))
print("time_taken: " + str(time_taken))
print("err_norm: " + str(gauss_seidel_err_norms))

plt.figure()
plt.plot(range(1, jacobi_iterations + 1), jacobi_err_norms, label='Jacobi Method')
plt.plot(range(1, gauss_seidel_iterations + 1), gauss_seidel_err_norms, label='Gauss-Seidel Method')
plt.yscale('log')
plt.xlabel('Iterations')
plt.ylabel('Error Norm')
plt.title('Error Norm Comparison - Jacobi vs Gauss-Seidel')
plt.legend()
plt.grid(True)

# Save the plot as an image
if not os.path.exists('img'):
    os.makedirs('img')
plt.savefig('img/zadanieB.png')

# ZADANIE C, D

a1, a2, a3 = 3, -1, -1
N = 948
A = create_matrix(N,a1,a2,a3)
b = Matrix(N,1,1)
for i in range(N):
    b[(i,0)] = sin((i+1)*4)

x, jacobi_err_norms, jacobi_iterations, time_taken = jacobi_iter_method(A, b, 100)
print("Jacobi method:")
print("iterations: " + str(jacobi_iterations))
print("time_taken: " + str(time_taken))
print("err_norm: " + str(jacobi_err_norms))

x, gauss_seidel_err_norms, gauss_seidel_iterations, time_taken = gauss_seidel_iter_method(A, b, 100)
print("Gauss-Seidel method:")
print("iterations: " + str(gauss_seidel_iterations))
print("time_taken: " + str(time_taken))
print("err_norm: " + str(gauss_seidel_err_norms))

x, err_norms, time_taken = direct_method(A, b)
print("Direct method:")
print("time_taken: " + str(time_taken))
print("err_norm: " + str(err_norms))

plt.figure()
plt.plot(range(1, jacobi_iterations + 1), jacobi_err_norms, label='Jacobi Method')
plt.plot(range(1, gauss_seidel_iterations + 1), gauss_seidel_err_norms, label='Gauss-Seidel Method')
plt.yscale('log')
plt.xlabel('Iterations')
plt.ylabel('Error Norm')
plt.title('Error Norm Comparison - Jacobi vs Gauss-Seidel')
plt.legend()
plt.grid(True)

# Save the plot as an image
if not os.path.exists('img'):
    os.makedirs('img')
plt.savefig('img/zadanieC.png')


# ZADANIE E

jacobi_times = []
gauss_seidel_times = []
direct_times = []

sizes = [100, 500, 1000, 2000, 3000, 4000, 5000]

for num in sizes:
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
    jacobi_times.append(time_taken)

    x, err_norms, iterations, time_taken = gauss_seidel_iter_method(A, b, 100)
    print("Gauss-Seidel method:")
    print("iterations: " + str(iterations))
    print("time_taken: " + str(time_taken))
    print("err_norm: " + str(err_norms))
    gauss_seidel_times.append(time_taken)


    x, err_norms, time_taken = direct_method(A, b)
    print("Direct method:")
    print("time_taken: " + str(time_taken))
    print("err_norm: " + str(err_norms))
    direct_times.append(time_taken)


plt.figure(figsize=(10, 6))
plt.plot(sizes, jacobi_times, marker='o', label='Jacobi Method')
plt.plot(sizes, gauss_seidel_times, marker='o', label='Gauss-Seidel Method')
plt.plot(sizes, direct_times, marker='o', label='Direct Method')
plt.yscale('log')
plt.xlabel('Matrix Size')
plt.ylabel('Time Taken (seconds)')
plt.title('Comparison of Time Taken by Each Method')
plt.xticks(sizes)
plt.grid(True)
plt.legend()

if not os.path.exists('img'):
    os.makedirs('img')
plt.savefig('img/zadanieElog.png')