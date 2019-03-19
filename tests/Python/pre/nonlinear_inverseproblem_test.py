#!/usr/bin/env python3

# WARNING: minieigen before pyBasso!!
# Otherwise it will be impossible to construct the Eigen matrix properly
# MatrixX.__init__() will say: did not match C++ signature: ...
import minieigen

import pyBasso
import numpy as np

np.random.seed(426)

X = pyBasso.create_LpSpace(10, 2., 2.)
Y = pyBasso.create_LpSpace(10, 2., 2.)

matrix = np.asarray(np.random.uniform(low=-1., high=1., size=(10,10)))
M = minieigen.MatrixX(matrix)

righthandside = np.asarray(np.random.uniform(low=-1., high=1., size=(10)))
rhs = minieigen.VectorX(righthandside)

#A = pyBasso.create_LinearMapping(X,Y, matrix, False)

def mapping(_argument):
	#print("map: "+str(_argument))
	return M * _argument

def derivative(_argument):
	#print("adjoint: "+str(_argument))
	return M.transpose() * _argument

A = pyBasso.create_NonLinearMapping(X,Y, mapping, derivative, False)

print(rhs)
y = Y.dualspace.createElement()
y.set(rhs)
print(y)

# make it such that y has a valid minimum-norm solution
J_p = X.dualspace.getDualityMapping()
Asys = A.adjoint(y)
print(Asys)
precursor = J_p(Asys)
print(precursor)
xmin = precursor*(1./precursor.norm)
print(xmin)
y = A(xmin)
print(y)
inverse_problem = pyBasso.InverseProblem(A, X, Y, y);

# create opts
opts = pyBasso.CommandLineOptions()
opts.algorithm_name = "SESOP"
opts.delta = 0.001
opts.C = 1
opts.maxiter = 1000
opts.stopping_criteria = "RelativeResiduum || MaxIterationCount"
opts.orthogonalization_type = 1 # MetricProjection
opts.verbose=2
opts.setVerbosity()
opts.setValues()

# create zero start value and solve
xstart = X.createElement()
inverse_problem.solve(opts, xstart)

print("Expected solution: "+str(xmin))
print("Approximate solution: "+str(inverse_problem.x))
print("Norm of difference: "+str((inverse_problem.x-xmin).norm))
