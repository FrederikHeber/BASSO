# test 1a: solving random matrix with Landweber
add_test (Basso_Landweber_RandomMatrix
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm Landweber
		--normx 2 
		--normy 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--C 1 
		--maxiter 100 
		--iteration-file Landweber.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/Landweber/RandomMatrix/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/Landweber/RandomMatrix/pre/randomvector.m
		--solution solution.m)
set_tests_properties (Basso_Landweber_RandomMatrix
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 100 with relative residual of 0.87.*")

# test 1b: checking results
add_test (Basso_Landweber_RandomMatrix_check
	diff ${PROJECT_SOURCE_DIR}/tests/Basso/Landweber/RandomMatrix/post/solution.m  solution.m)
set_tests_properties (Basso_Landweber_RandomMatrix_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_Landweber_RandomMatrix_check
	PROPERTIES DEPENDS Basso_Landweber_RandomMatrix )
	
	