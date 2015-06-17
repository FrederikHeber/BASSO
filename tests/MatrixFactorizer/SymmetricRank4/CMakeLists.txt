# does simple test call work?
add_test (MatrixFactorizer_SymmetricRank4
	${PROJECT_BINARY_DIR}/src/MatrixFactorizer 
		--normx 2 
		--normy 2 
		--powery 2 
		--delta 0.01 
		--number-directions 1
		--sparse-dim 1
		--data ${PROJECT_SOURCE_DIR}/tests/MatrixFactorizer/SymmetricRank4/pre/data.m 
		--solution-first-factor first-factor.m
		--solution-second-factor second-factor.m)
#set_tests_properties (MatrixFactorizer_SymmetricRank4
#	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of .*5.9")

# test 1b: checking results
add_test (MatrixFactorizer_SymmetricRank4_check_first
	diff ${PROJECT_SOURCE_DIR}/tests/MatrixFactorizer/SymmetricRank4/post/first-factor.m first-factor.m)
set_tests_properties (MatrixFactorizer_SymmetricRank4_check_first
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_SymmetricRank4_check_first
	PROPERTIES DEPENDS MatrixFactorizer_SymmetricRank4)

# test 1c: checking results
add_test (MatrixFactorizer_SymmetricRank4_check_second
	diff ${PROJECT_SOURCE_DIR}/tests/MatrixFactorizer/SymmetricRank4/post/second-factor.m second-factor.m)
set_tests_properties (MatrixFactorizer_SymmetricRank4_check_second
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_SymmetricRank4_check_second
	PROPERTIES DEPENDS MatrixFactorizer_SymmetricRank4)
