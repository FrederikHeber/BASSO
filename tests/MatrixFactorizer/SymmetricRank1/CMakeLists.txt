# does simple test call work?
add_test (MatrixFactorizer_SymmetricRank1 
	${PROJECT_BINARY_DIR}/src/MatrixFactorizer 
		--normx 2 
		--normy 2 
		--powery 2 
		--delta 0.01 
		--number-directions 1
		--sparse-dim 1
		--data ${PROJECT_SOURCE_DIR}/tests/MatrixFactorizer/SymmetricRank1/pre/data.m 
		--solution-first-factor first-factor.m
		--solution-second-factor first-factor.m)
#set_tests_properties (MatrixFactorizer_SymmetricRank1
#	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of X.*2")
