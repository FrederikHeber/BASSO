# test 1: 
add_test (Basso_SESOP_RandomMatrix_NemirovskyDirections
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso
		--algorithm SESOP
		--number-directions 2
		--searchspace "Nemirovsky"
		--update-algorithm 0
		--type-space-x "lp" 
		--px 2
		--type-space-y "lp" 
		--py 0
		--powerx 2
		--powery 2
		--delta 0.01
		--maxiter 100
		--iteration-file SESOP1.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--solution solution.m)
set_tests_properties (Basso_SESOP_RandomMatrix_NemirovskyDirections
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 24 with relative residual of 0.0095.*")

# test 1b: checking results
add_test (Basso_SESOP_RandomMatrix_NemirovskyDirections_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution.m  solution.m)
set_tests_properties (Basso_SESOP_RandomMatrix_NemirovskyDirections_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_NemirovskyDirections_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_NemirovskyDirections)

# test 2: Nemirovsky with not two search directions
add_test (Basso_SESOP_RandomMatrix_NemirovskyDirections_fail
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso
		--algorithm SESOP
		--number-directions 4
		--searchspace "Nemirovsky"
		--update-algorithm 0
		--type-space-x "lp" 
		--px 2
		--type-space-y "lp" 
		--py 0
		--powerx 2
		--powery 2
		--delta 0.01
		--maxiter 100
		--iteration-file SESOP2.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m)
set_tests_properties (Basso_SESOP_RandomMatrix_NemirovskyDirections_fail
	PROPERTIES WILL_FAIL "true")

