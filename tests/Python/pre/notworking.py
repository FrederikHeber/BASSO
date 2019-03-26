#!/usr/bin/env python3

# WARNING: minieigen before pyBasso!!
# Otherwise it will be possible to construct the Eigen matrix properly
# MatrixX.__init__() will say: did not match C++ signature: ...
import minieigen

import pyBasso
import numpy as np

X = pyBasso.create_LpSpace(10, 2., 2.)
Y = pyBasso.create_LpSpace(12, 2., 2.)

x = X.createElement()
x[0] = 1.

matrix = np.asarray(np.zeros((12,10)))

M = minieigen.MatrixX(matrix)
M[5,0] = 1.
print(M)

A = pyBasso.create_LinearMapping(X,Y, M, False)
print(X.dim)
print(Y.dim)
At = A.adjoint
print(A.sourcespace.dim)
print(A.targetspace.dim)
print(At.sourcespace.dim)
print(At.targetspace.dim)

print(x.space.dim)
Ax = A * x
print(Ax)
print(Ax.space.dim)
print(At*Ax)
