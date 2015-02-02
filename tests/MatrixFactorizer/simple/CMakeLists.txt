# does simple test call fail when matrix does not exist?
add_test (MatrixFactorizerSimple_NonExistantMatrix 
	${PROJECT_BINARY_DIR}/src/MatrixFactorizer 
		--normx 2 
		--normy 2 
		--powery 2 
		--delta 0.01
		--sparse-dim 1
		--data ${PROJECT_SOURCE_DIR}/tests/MatrixFactorizer/simple/pre/nomatrix.m 
		--solution-first-factor first-factor.m
		--solution-second-factor second-factor.m)
set_tests_properties (MatrixFactorizerSimple_NonExistantMatrix
	PROPERTIES WILL_FAIL true)
	