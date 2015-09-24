# test 1: with Wolfe conditions
add_test (Basso_SESOP_RandomMatrix_InexactLinesearch_Wolfe
	${PROJECT_BINARY_DIR}/Examples/BassoBase/Basso 
		--algorithm SESOP
		--number-directions 2
		--update-algorithm 0
		--inexact-linesearch 1
		--wolfe-constants 1e-2 0.6
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--iteration-file SESOP.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/InexactLineSearch/WolfeConditions/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/InexactLineSearch/WolfeConditions/pre/randomvector.m
		--solution solution.m)
set_tests_properties (Basso_SESOP_RandomMatrix_InexactLinesearch_Wolfe
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 33 with relative residual of 0.008.*")

# test 1b: checking results
add_test (Basso_SESOP_RandomMatrix_InexactLinesearch_Wolfe_check
	diff ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/InexactLineSearch/WolfeConditions/post/solution.m solution.m)
set_tests_properties (Basso_SESOP_RandomMatrix_InexactLinesearch_Wolfe_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_InexactLinesearch_Wolfe_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_InexactLinesearch_Wolfe)
	