#!/usr/bin/python
import math
import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt


def Sparse_CRS_MVMult(N, row, col, val, v):
    result = [0] * N
    for i in range(0, N):
        for k in range(row[i], row[i + 1]):
            result[i] = result[i] + val[k] * v[col[k]]
    return result


N = 100
N1 = N - 2
N2 = N1 ** 2
dx = dy = 1.0 / (N - 1)
bet = dx / dy
alf = -2 * (1 + bet ** 2)

row = []
row_m = []
col = []
col_m = []
val = []
val_m = []
k = 0
k_m = 0

for i in range(0, N2):
    row.append(k)
    row_m.append(k_m)

    for j in range(0, N2):
        if (i == j):

            val.append(-4)
            col.append(j)
            k += 1
            val_m.append(1 / alf)
            col_m.append(j)
            k_m += 1

        elif (i == j + 1):

            if (i != 0 and i % (N1) == 0):

                val.append(0)
                col.append(j)
                k += 1

            else:

                val.append(1)
                col.append(j)
                k += 1


        elif (i + 1 == j):

            if (j != 0 and j % (N1) == 0):

                val.append(0)
                col.append(j)
                k += 1

            else:

                val.append(1)
                col.append(j)
                k += 1


        elif (i == j + N1):

            val.append(pow(bet, 2))
            col.append(j)
            k += 1

        elif (i + N1 == j):

            val.append(pow(bet, 2))
            col.append(j)
            k += 1

row.append(k)
row_m.append(k_m)

x = []
y = []
u = []
for i in range(0, N):
    u.append([0] * N)

for i in range(0, N):
    x.append(i * dx)
    y.append(i * dx)
    u[i][0] = math.sin(math.pi * x[i])
    u[i][N - 1] = math.sin(math.pi * x[i]) * math.exp(-math.pi)
    u[0][i] = 0
    u[N - 1][i] = 0

k = 0

B = []

for i in range(1, N - 1):
    for j in range(1, N - 1):
        B.append(-u[i + 1][j] - u[i - 1][j] - pow(bet, 2) * u[i][j + 1] - pow(bet, 2) * u[i][j - 1])
        k += 1

s = [0] * N2
s = Sparse_CRS_MVMult(N2, row, col, val, s)

r = [0] * N2
for i in range(0, N2):
    r[i] = B[i] - s[i]

z = Sparse_CRS_MVMult(N2, row_m, col_m, val_m, r)

p = z

rho = np.dot(r, z)

niter = 0

tol = pow(10,-16)

maxiter = 10000
U = [0] * N2
normB = np.linalg.norm(B)
while ((np.linalg.norm(r) / normB) > tol):  # Test break condition

    a = Sparse_CRS_MVMult(N2, row, col, val, p);
    alpha = rho / np.dot(a, p);
    U = scipy.linalg.blas.daxpy(p, U, N2, alpha, 0, 1, 0, 1)
    r = scipy.linalg.blas.daxpy(a, r, N2, -alpha, 0, 1, 0, 1)
    z = Sparse_CRS_MVMult(N2, row_m, col_m, val_m, r)
    rho_new = np.dot(r, z)

    for i in range(0, N2):
        p[i] = z[i] + (rho_new / rho) * p[i]

    rho = rho_new
    niter += 1
    if (niter == maxiter):  # if max. number of iterations  is reached, break
        print("ERROR!!!! \n Max Iterations Reached!!!!\n")
    	break

k = 0
for i in range(1, N - 1):
    for j in range(1, N - 1):
        u[i][j] = U[k]
    	k += 1

plt.contourf(u)
plt.show()
