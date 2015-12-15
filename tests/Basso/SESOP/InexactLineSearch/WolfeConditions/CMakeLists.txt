if (NLOPT_FOUND)
	# test 1: with Wolfe conditions
	add_test (Basso_SESOP_RandomMatrix_InexactLinesearch_Wolfe
		${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
			--algorithm SESOP
			--number-directions 2
			--update-algorithm 0
			--minimization-library nlopt
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
			--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
			--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
			--solution solution.m)
	set_tests_properties (Basso_SESOP_RandomMatrix_InexactLinesearch_Wolfe
		PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 31 with relative residual of 0.007.*")
	
	# test 1b: checking results
	add_test (Basso_SESOP_RandomMatrix_InexactLinesearch_Wolfe_check
		${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution.m solution.m)
	set_tests_properties (Basso_SESOP_RandomMatrix_InexactLinesearch_Wolfe_check
		PROPERTIES WILL_FAIL false)
	set_tests_properties (Basso_SESOP_RandomMatrix_InexactLinesearch_Wolfe_check
		PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_InexactLinesearch_Wolfe)
endif (NLOPT_FOUND)
