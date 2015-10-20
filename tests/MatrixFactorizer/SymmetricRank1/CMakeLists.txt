# test 1: Matrix factorization with rank 1 
add_test (MatrixFactorizer_SymmetricRank1 
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01 
		--maxiter 50
		--number-directions 1
		--sparse-dim 1
		--data ${PROJECT_SOURCE_DIR}/tests/MatrixFactorizer/SymmetricRank1/pre/data.m 
		--solution-first-factor first-factor.m
		--solution-second-factor second-factor.m)
set_tests_properties (MatrixFactorizer_SymmetricRank1
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of .*2.32")

# test 1b: checking results
add_test (MatrixFactorizer_SymmetricRank1_check_first
	${DIFFCMD} ${NDIFF_OPTIONS} ${PROJECT_SOURCE_DIR}/tests/MatrixFactorizer/SymmetricRank1/post/first-factor.m first-factor.m)
set_tests_properties (MatrixFactorizer_SymmetricRank1_check_first
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_SymmetricRank1_check_first
	PROPERTIES DEPENDS MatrixFactorizer_SymmetricRank1)

# test 1c: checking results
add_test (MatrixFactorizer_SymmetricRank1_check_second
	${DIFFCMD} ${NDIFF_OPTIONS} ${PROJECT_SOURCE_DIR}/tests/MatrixFactorizer/SymmetricRank1/post/second-factor.m second-factor.m)
set_tests_properties (MatrixFactorizer_SymmetricRank1_check_second
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_SymmetricRank1_check_second
	PROPERTIES DEPENDS MatrixFactorizer_SymmetricRank1)

# test 2: Matrix factorization with rank 1 but more than one search direction
add_test (MatrixFactorizer_SymmetricRank1_fail 
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01 
		--number-directions 2
		--sparse-dim 1
		--data ${PROJECT_SOURCE_DIR}/tests/MatrixFactorizer/SymmetricRank1/pre/data.m)
set_tests_properties (MatrixFactorizer_SymmetricRank1_fail
	PROPERTIES WILL_FAIL true)
