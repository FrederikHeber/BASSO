#!/usr/bin/env python3

import pyBasso
import numpy as np

X = pyBasso.create_LpSpace(10, 2., 2.)
print(X.dim)

x = X.createElement()

norm = X.getNorm()
print(norm(x))
print(norm.p)
print(X.space.dim)
print(X.dualspace.dim)
print(X.getDualityMapping().power)
