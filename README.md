libBasso - BAnach Sequential Subspace Optimization library
==========================================================

Basso is a library to minimize a constrained functional in a Banach space
setting. It is an extension of the sequential subspace optimization
method by Narkiss by multiple search directions.

For more details we refer to the extensive userguide.

Installation instructions:
--------------------------

Basso depends on the following packages:
 * cmake - make utility for checking prerequisites
 * GSL - GNU Scientific library for function minimization
 * NLopt - (optional) alternative (to GSL) function minimization library
 * Eigen - C++ library for linear algebra
 * ARPack - library for linear algebra routines
 * CppUnit - C++ library for unit testing
 * boost - almost standard library (we require 1.54 due to logging)
 * poco - poco is used (optionally) to write iteration information to
   a sqlite3 database
 * libpng - library for writing png files (used by a helper program)

... on Ubuntu
-------------
Eigen and CppUnit packages and reside in standard paths, hence we have to do nothing

It is advised to compile the source in a distinct build directory, e.g. in a subfolder "build64". This is referred to as out-of-source build. To prepare the compilation, call cmake then as follows:

    cmake [options] ..

These *options* might be one of the following. Take caution to prepend them
before the path containing the **CMakeLists.txt** file!

### using boost


If additionally boost libraries do not reside in standard folders, such as
/usr/include, /usr/lib ... but in ~/packages/boost the following options might help:

    .. -DBoost_NO_SYSTEM_PATHS=ON -DBOOST_ROOT=<path to boost>/boost/

Note that boost libraries should have been compiled with "layout=tagged". Otherwise, FindBoost.cmake can not find the libraries and headers.

To allow for debugging, you might futher add "-DBoost_DEBUG=ON".

### using Eigen

If Eigen resides in some specific folder, add it as follows

    .. -DEIGEN3_INCLUDE_DIR=<path to eigen>/include/eigen3

### Debug or Release builds

In general, to compile in debugging information add

    -DCMAKE_BUILD_TYPE=Debug

to compile without use

    -DCMAKE_BUILD_TYPE=Release


Documentation
-------------

To generate the doxygen documenation, you use

    make doc

The generated API html can the be found under **src/Documentation/doc/html/index.html**.

Moreover,

   cd doc; make

generates the userguide in **doc/basso.pdf** and **doc/basso.html**.

Timing Measurements
-------------------

Add

    -DUSE_TIMINGS=true

to activate the timing measurements of each part of the algorithm.


Parallelization
---------------

The example MatrixFactorizer has been implemented for parallel execution with
both MPI and OpenMP. However, only one of the two can be activated during
compilation.

For **MPI** add

    -DUSE_MPI=True

and for **OpenMP** add

    -DUSE_OPENMP=True

to the cmake command line. Note that only one of them may be True at the same
time (but both can be False in which case factorization is done sequentially.)

Special purpose diff
--------------------

Result files are compared against stored files in the regression tests using
diff by default. `diff` compares literally and does not know about numerical
(im)precision. Hence, a different diff, such as ndiff, may be specified via

    -DUSE_NDIFF=ON -DNDIFF_PATH=<path to ndiff's bin folder>
    -DNDIFF_OPTIONS:LIST="-relerr;1e-4;-q"

to activate comparison of outputfiles with ndiff. The semicolons ";" are
required to ensure that we have a list of options (otherwise they are presented
to ndiff as single string which will not work).

Note: There is another `ndiff` program used to compare XML output files from
nmap. Hence, the path NDIFF_PATH to ndiff's bin folder should always be given.


Various Tests
-------------

Various compile switches exist to enable more testing of the algorithm.

TRUESOLUTION			Calculates true solution in inner MatrixFactorizer loop to calculate Bregman distance and check convergence (see src/MatrixFactorizer.cpp).
FULLMATRIXNORM			Calculates operator norm via SVD instead of matrix norm (takes O(N^3) !) (cmake define or compile switch)


2015-03-18 Frederik Heber
