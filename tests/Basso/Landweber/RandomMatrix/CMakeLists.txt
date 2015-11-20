# test 1a: solving random matrix with Landweber
add_test (Basso_Landweber_RandomMatrix
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm Landweber
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--C 1 
		--maxiter 100 
		--iteration-file Landweber.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--solution solution.m)
set_tests_properties (Basso_Landweber_RandomMatrix
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 100 with relative residual of 0.87.*")

# test 1b: checking results
add_test (Basso_Landweber_RandomMatrix_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution.m  solution.m)
set_tests_properties (Basso_Landweber_RandomMatrix_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_Landweber_RandomMatrix_check
	PROPERTIES DEPENDS Basso_Landweber_RandomMatrix )
	
	