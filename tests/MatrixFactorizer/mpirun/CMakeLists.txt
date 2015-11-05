# test 1: Matrix factorization using mpirun with rank 1 
add_test (MatrixFactorizer_mpirun
	${MPIEXEC} -${MPIEXEC_NUMPROC_FLAG} 3
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer
		--verbose 0
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01 
		--maxiter 50
		--max-loops 10
		--number-directions 1
		--sparse-dim 1
		--data ${PROJECT_SOURCE_DIR}/tests/MatrixFactorizer/mpirun/pre/data.m 
		--solution-first-factor first-factor.m
		--solution-second-factor second-factor.m)
set_tests_properties (MatrixFactorizer_mpirun
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_mpirun
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of .*5.9")

# test 1b: checking results
add_test (MatrixFactorizer_mpirun_check_first
	${DIFFCMD} ${NDIFF_OPTIONS} ${PROJECT_SOURCE_DIR}/tests/MatrixFactorizer/mpirun/post/first-factor.m first-factor.m)
set_tests_properties (MatrixFactorizer_mpirun_check_first
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_mpirun_check_first
	PROPERTIES DEPENDS MatrixFactorizer_mpirun)

# test 1c: checking results
add_test (MatrixFactorizer_mpirun_check_second
	${DIFFCMD} ${NDIFF_OPTIONS} ${PROJECT_SOURCE_DIR}/tests/MatrixFactorizer/mpirun/post/second-factor.m second-factor.m)
set_tests_properties (MatrixFactorizer_mpirun_check_second
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_mpirun_check_second
	PROPERTIES DEPENDS MatrixFactorizer_mpirun)
