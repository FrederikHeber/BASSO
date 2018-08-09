#!/usr/bin/env python3

# WARNING: minieigen before pyBasso!!
# Otherwise it will be impossible to construct the Eigen matrix properly
# MatrixX.__init__() will say: did not match C++ signature: ...
import minieigen

import pyBasso
import numpy as np

np.random.seed(426)

X = pyBasso.create_LpSpace(10, 2., 2.)
Y = pyBasso.create_LpSpace(11, 2., 2.)

x = X.createElement()
x[0] = 1.

matrix = np.asarray(np.random.uniform(low=-1., high=1., size=(Y.dim,X.dim)))

righthandside = np.asarray(np.random.uniform(low=-1., high=1., size=(Y.dim)))

M = minieigen.MatrixX(matrix)
#print(M)
A = pyBasso.create_LinearMapping(X,Y, M, False)

rhs = minieigen.VectorX(righthandside)
print(rhs)
y = Y.createElement()
y.set(rhs)
print(y)

inverse_problem = pyBasso.InverseProblem(A, X, Y, y);

# create opts
opts = pyBasso.Options()
opts.algorithm_name = "SESOP"
opts.delta = 0.001
opts.C = 1
opts.maxiter = 10
opts.stopping_criteria = "RelativeResiduum || MaxIterationCount"
opts.orthogonalization_type = 1 # MetricProjection
opts.setValues()
opts.setVerbosity(0)

# create zero start value and project onto range
xstart = Y.dualspace.createElement()
inverse_problem.project(opts, xstart)

print("Projected rhs: "+str(inverse_problem.y))
print("Norm difference: "+str((inverse_problem.y-y).norm))
