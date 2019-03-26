#!/usr/bin/env python3

import pyBasso
import numpy as np

X = pyBasso.create_LpSpace(10, 2., 2.)

x = X.createElement()
x.setZero()

print(x.norm)
print(x.space.dim)
print(x.isZero())
print(x.isZero(1e-1))

x[1] = 1.
print(x[1])

print(x.isZero())
print(x.isZero(1e+1))

y = X.createElement()

y[0] = 1.

z = x-y

print(z.isNonnegative())

print((z-z.getSignVector()).norm)

print(z.norm)

print(z.getAbsVector())

