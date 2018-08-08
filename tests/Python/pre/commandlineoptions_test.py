#!/usr/bin/env python3

import pyBasso

opts = pyBasso.CommandLineOptions()
opts.algorithm_name = "Landweber"
print(opts.algorithm_name)
print("N="+str(opts.N))
print(opts.delta)

# TODO: Options::store looks at vm not at actual values.
# Override CommandLineOptions and reimplement writeValue<> class.
print(opts)
