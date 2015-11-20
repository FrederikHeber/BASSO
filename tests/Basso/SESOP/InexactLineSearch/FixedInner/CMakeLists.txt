# test 1a: with fixed inner iterations
add_test (Basso_SESOP_RandomMatrix_InexactLinesearch_fixedinner
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
		--update-algorithm 0
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--iteration-file SESOP.db
		--max-inner-iterations 10
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--solution solution.m)
set_tests_properties (Basso_SESOP_RandomMatrix_InexactLinesearch_fixedinner
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 18 with relative residual of 0.0098.*")

# test 1b: checking results
add_test (Basso_SESOP_RandomMatrix_InexactLinesearch_fixedinner_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution.m solution.m)
set_tests_properties (Basso_SESOP_RandomMatrix_InexactLinesearch_fixedinner_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_InexactLinesearch_fixedinner_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_InexactLinesearch_fixedinner)
	