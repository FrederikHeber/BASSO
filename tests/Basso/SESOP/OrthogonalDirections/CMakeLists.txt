# test 1a: 
add_test (Basso_SESOP_RandomMatrix_p2_r2_orthogonalized_N1
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 1
		--update-algorithm 0
		--orthogonal-directions 1
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--iteration-file SESOP_p2_r2_N1.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--solution solution_p2_r2_N1.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_r2_orthogonalized_N1
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 4 with relative residual of 0.0093.*")

# test 1b: checking results
add_test (Basso_SESOP_RandomMatrix_p2_r2_orthogonalized_N1_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p2_r2.m solution_p2_r2_N1.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_r2_orthogonalized_N1_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_r2_orthogonalized_N1_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p2_r2_orthogonalized_N1)

# test 1c: 
add_test (Basso_SESOP_RandomMatrix_p2_r2_orthogonalized_N2
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
		--update-algorithm 0
		--orthogonal-directions 1
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--iteration-file SESOP_p2_r2_N2.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--solution solution_p2_r2_N2.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_r2_orthogonalized_N2
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 4 with relative residual of 0.0093.*")

# test 1d: checking results
add_test (Basso_SESOP_RandomMatrix_p2_r2_orthogonalized_N2_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p2_r2.m solution_p2_r2_N2.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_r2_orthogonalized_N2_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_r2_orthogonalized_N2_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p2_r2_orthogonalized_N2)

# test 1e: 
add_test (Basso_SESOP_RandomMatrix_p2_r2_orthogonalized_N4
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 4
		--update-algorithm 0
		--orthogonal-directions 1
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--iteration-file SESOP_p2_r2_N4.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--solution solution_p2_r2_N4.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_r2_orthogonalized_N4
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 4 with relative residual of 0.0093.*")

# test 1f: checking results
add_test (Basso_SESOP_RandomMatrix_p2_r2_orthogonalized_N4_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p2_r2.m solution_p2_r2_N4.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_r2_orthogonalized_N4_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_r2_orthogonalized_N4_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p2_r2_orthogonalized_N4)

# test 2a: 
add_test (Basso_SESOP_RandomMatrix_p2_r11_orthogonalized
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
		--orthogonal-directions 1
		--update-algorithm 0
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 11
		--powerx 2
		--powery 11
		--delta 0.01 
		--maxiter 100 
		--iteration-file SESOP_p2_r11.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--solution solution_p2_r11.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_r11_orthogonalized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 10 with relative residual of 0.034.*")

# test 2b: checking results
add_test (Basso_SESOP_RandomMatrix_p2_r11_orthogonalized_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p2_r11.m solution_p2_r11.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_r11_orthogonalized_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_r11_orthogonalized_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p2_r11_orthogonalized)

# test 3a: 
add_test (Basso_SESOP_RandomMatrix_p2_rinf_orthogonalized
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
		--orthogonal-directions 1
		--update-algorithm 0
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--iteration-file SESOP_p2_rinf.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--solution solution_p2_rinf.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_rinf_orthogonalized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 20 with relative residual of 0.005.*")

# test 3b: checking results
add_test (Basso_SESOP_RandomMatrix_p2_rinf_orthogonalized_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p2_rinf.m solution_p2_rinf.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_rinf_orthogonalized_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_rinf_orthogonalized_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p2_rinf_orthogonalized)

# test 4a: 
add_test (Basso_SESOP_RandomMatrix_p1_1_r2_orthogonalized
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
		--orthogonal-directions 1
		--update-algorithm 0
		--type-space-x "lp" 
		--px 1.1 
		--type-space-y "lp" 
		--py 2 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--iteration-file SESOP_p1_1_r2.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--solution solution_p1_1_r2.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p1_1_r2_orthogonalized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 18 with relative residual of 0.0079.*")

# test 4b: checking results
add_test (Basso_SESOP_RandomMatrix_p1_1_r2_orthogonalized_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p1_1_r2.m solution_p1_1_r2.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p1_1_r2_orthogonalized_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p1_1_r2_orthogonalized_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p1_1_r2_orthogonalized)

# test 5a: 
add_test (Basso_SESOP_RandomMatrix_p1_1_r11_orthogonalized
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
		--orthogonal-directions 1
		--update-algorithm 0
		--type-space-x "lp" 
		--px 1.1
		--type-space-y "lp" 
		--py 11
		--powerx 2
		--powery 11
		--delta 0.01 
		--maxiter 100 
		--iteration-file SESOP_p1_1_r11.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--solution solution_p1_1_r11.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p1_1_r11_orthogonalized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 21 with relative residual of 0.042.*")

# test 5b: checking results
add_test (Basso_SESOP_RandomMatrix_p1_1_r11_orthogonalized_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p1_1_r11.m solution_p1_1_r11.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p1_1_r11_orthogonalized_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p1_1_r11_orthogonalized_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p1_1_r11_orthogonalized)

# test 6a: 
add_test (Basso_SESOP_RandomMatrix_p1_1_rinf_orthogonalized
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
		--orthogonal-directions 1
		--update-algorithm 0
		--type-space-x "lp" 
		--px 1.1
		--type-space-y "lp" 
		--py 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--iteration-file SESOP_p1_1_rinf.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--solution solution_p1_1_rinf.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p1_1_rinf_orthogonalized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 67 with relative residual of 0.0094.*")

# test 6b: checking results
add_test (Basso_SESOP_RandomMatrix_p1_1_rinf_orthogonalized_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p1_1_rinf.m solution_p1_1_rinf.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p1_1_rinf_orthogonalized_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p1_1_rinf_orthogonalized_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p1_1_rinf_orthogonalized)

# test 7a: 
add_test (Basso_SESOP_RandomMatrix_p6_r2_orthogonalized
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
		--orthogonal-directions 1
		--update-algorithm 0
		--type-space-x "lp" 
		--px 6 
		--type-space-y "lp" 
		--py 2 
		--powerx 6
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--iteration-file SESOP_p6_r2.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--solution solution_p6_r2.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p6_r2_orthogonalized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 22 with relative residual of 0.0083.*")

# test 7b: checking results
add_test (Basso_SESOP_RandomMatrix_p6_r2_orthogonalized_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p6_r2.m solution_p6_r2.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p6_r2_orthogonalized_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p6_r2_orthogonalized_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p6_r2_orthogonalized)

# test 8a: 
add_test (Basso_SESOP_RandomMatrix_p6_r11_orthogonalized
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
		--orthogonal-directions 1
		--update-algorithm 0
		--type-space-x "lp" 
		--px 6 
		--type-space-y "lp" 
		--py 11
		--powerx 6
		--powery 11
		--delta 0.01 
		--maxiter 100 
		--iteration-file SESOP_p6_r11.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--solution solution_p6_r11.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p6_r11_orthogonalized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 31 with relative residual of 0.05.*")

# test 8b: checking results
add_test (Basso_SESOP_RandomMatrix_p6_r11_orthogonalized_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p6_r11.m solution_p6_r11.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p6_r11_orthogonalized_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p6_r11_orthogonalized_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p6_r11_orthogonalized)

# test 9a: 
add_test (Basso_SESOP_RandomMatrix_p6_rinf_orthogonalized
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
		--orthogonal-directions 1
		--update-algorithm 0
		--type-space-x "lp" 
		--px 6 
		--type-space-y "lp" 
		--py 0 
		--powerx 6
		--powery 2
		--delta 0.01 
		--maxiter 200 
		--iteration-file SESOP_p6_rinf.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--solution solution_p6_rinf.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p6_rinf_orthogonalized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 79 with relative residual of 0.0094.*")

# test 9b: checking results
add_test (Basso_SESOP_RandomMatrix_p6_rinf_orthogonalized_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p6_rinf.m solution_p6_rinf.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p6_rinf_orthogonalized_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p6_rinf_orthogonalized_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p6_rinf_orthogonalized)
	