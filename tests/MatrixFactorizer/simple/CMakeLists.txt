# does simple test call fail when matrix does not exist?
add_test (MatrixFactorizerSimple_NonExistantMatrix 
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01
		--sparse-dim 1
		--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/nomatrix.m 
		--solution-first-factor first-factor.m
		--solution-second-factor second-factor.m)
set_tests_properties (MatrixFactorizerSimple_NonExistantMatrix
	PROPERTIES WILL_FAIL true)
	
# does simple test succeed with given stopping criteria
add_test (MatrixFactorizerSimple_StoppingCriteria 
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01
		--sparse-dim 1
		--number-directions 1
		--max-loops 10
		--factorization-stopping-criteria "MaxIterationCount || Residuum"
		--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/matrix.m 
		--solution-first-factor first-factor.m
		--solution-second-factor second-factor.m)
set_tests_properties (MatrixFactorizerSimple_StoppingCriteria
	PROPERTIES WILL_FAIL false)
	
# does simple test succeed with given stopping criteria
add_test (MatrixFactorizerSimple_EmptyStoppingCriteria 
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01
		--sparse-dim 1
		--number-directions 1
		--factorization-stopping-criteria ""
		--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/matrix.m 
		--solution-first-factor first-factor.m
		--solution-second-factor second-factor.m)
set_tests_properties (MatrixFactorizerSimple_EmptyStoppingCriteria
	PROPERTIES WILL_FAIL true)
	
	