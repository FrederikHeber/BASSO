# test 1a: with fixed inner iterations
add_test (Basso_SESOP_RandomMatrix_InexactLinesearch_fixedinner
	${PROJECT_BINARY_DIR}/src/BassoBase/Basso 
		--algorithm SESOP
		--number-directions 2
		--update-algorithm 0
		--normx 2 
		--normy 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--iteration-file SESOP.db
		--max-inner-iterations 10
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/InexactLineSearch/FixedInner/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/InexactLineSearch/FixedInner/pre/randomvector.m
		--solution solution.m)
set_tests_properties (Basso_SESOP_RandomMatrix_InexactLinesearch_fixedinner
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 18 with relative residual of 0.0098.*")

# test 1b: checking results
add_test (Basso_SESOP_RandomMatrix_InexactLinesearch_fixedinner_check
	diff ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/InexactLineSearch/FixedInner/post/solution.m solution.m)
set_tests_properties (Basso_SESOP_RandomMatrix_InexactLinesearch_fixedinner_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_InexactLinesearch_fixedinner_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_InexactLinesearch_fixedinner)
	