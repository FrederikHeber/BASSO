# does simple test call work?
add_test (MatrixFactorizer_SymmetricRank2
	${PROJECT_BINARY_DIR}/src/MatrixFactorizerBase/MatrixFactorizer 
		--normx 2 
		--normy 2 
		--powery 2 
		--delta 0.01 
		--number-directions 1
		--sparse-dim 1
		--data ${PROJECT_SOURCE_DIR}/tests/MatrixFactorizer/SymmetricRank2/pre/data.m 
		--solution-first-factor first-factor.m
		--solution-second-factor second-factor.m)
set_tests_properties (MatrixFactorizer_SymmetricRank2
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of .*3.92")

# test 1b: checking results
add_test (MatrixFactorizer_SymmetricRank2_check_first
	diff ${PROJECT_SOURCE_DIR}/tests/MatrixFactorizer/SymmetricRank2/post/first-factor.m first-factor.m)
set_tests_properties (MatrixFactorizer_SymmetricRank2_check_first
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_SymmetricRank2_check_first
	PROPERTIES DEPENDS MatrixFactorizer_SymmetricRank2)

# test 1c: checking results
add_test (MatrixFactorizer_SymmetricRank2_check_second
	diff ${PROJECT_SOURCE_DIR}/tests/MatrixFactorizer/SymmetricRank2/post/second-factor.m second-factor.m)
set_tests_properties (MatrixFactorizer_SymmetricRank2_check_second
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_SymmetricRank2_check_second
	PROPERTIES DEPENDS MatrixFactorizer_SymmetricRank2)
