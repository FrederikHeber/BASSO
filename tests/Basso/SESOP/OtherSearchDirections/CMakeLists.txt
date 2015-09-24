# test 1: 
add_test (Basso_SESOP_RandomMatrix_NemirovskyDirections
	${PROJECT_BINARY_DIR}/Examples/BassoBase/Basso
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
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/OtherSearchDirections/pre/randommatrix.m
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/OtherSearchDirections/pre/randomvector.m
		--solution solution.m)
set_tests_properties (Basso_SESOP_RandomMatrix_NemirovskyDirections
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 24 with relative residual of 0.0095.*")

# test 1b: checking results
add_test (Basso_SESOP_RandomMatrix_NemirovskyDirections_check
	diff ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/OtherSearchDirections/post/solution.m  solution.m)
set_tests_properties (Basso_SESOP_RandomMatrix_NemirovskyDirections_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_NemirovskyDirections_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_NemirovskyDirections)

# test 2: Nemirovsky with not two search directions
add_test (Basso_SESOP_RandomMatrix_NemirovskyDirections_fail
	${PROJECT_BINARY_DIR}/Examples/BassoBase/Basso
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
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/OtherSearchDirections/pre/randommatrix.m
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/OtherSearchDirections/pre/randomvector.m)
set_tests_properties (Basso_SESOP_RandomMatrix_NemirovskyDirections_fail
	PROPERTIES WILL_FAIL "true")

