# test 1a: 
add_test (Basso_SESOP_RandomMatrix_p2_r2
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
		--update-algorithm 0
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--iteration-file SESOP_p2_r2.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--solution solution_p2_r2.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_r2
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 5 with relative residual of 0.0066.*")

# test 1b: checking results
add_test (Basso_SESOP_RandomMatrix_p2_r2_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p2_r2.m solution_p2_r2.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_r2_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_r2_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p2_r2)

# test 2a: 
add_test (Basso_SESOP_RandomMatrix_p2_r11
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
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
set_tests_properties (Basso_SESOP_RandomMatrix_p2_r11
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 12 with relative residual of 0.0386.*")

# test 2b: checking results
add_test (Basso_SESOP_RandomMatrix_p2_r11_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p2_r11.m solution_p2_r11.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_r11_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_r11_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p2_r11)

# test 3a: 
add_test (Basso_SESOP_RandomMatrix_p2_rinf
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
		--iteration-file SESOP_p2_rinf.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--solution solution_p2_rinf.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_rinf
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 18 with relative residual of 0.009.*")

# test 3b: checking results
add_test (Basso_SESOP_RandomMatrix_p2_rinf_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p2_rinf.m solution_p2_rinf.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_rinf_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p2_rinf_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p2_rinf)

# test 4a: 
add_test (Basso_SESOP_RandomMatrix_p1_1_r2
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
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
set_tests_properties (Basso_SESOP_RandomMatrix_p1_1_r2
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 24 with relative residual of 0.0079.*")

# test 4b: checking results
add_test (Basso_SESOP_RandomMatrix_p1_1_r2_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p1_1_r2.m solution_p1_1_r2.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p1_1_r2_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p1_1_r2_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p1_1_r2)

# test 5a: 
add_test (Basso_SESOP_RandomMatrix_p1_1_r11
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
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
set_tests_properties (Basso_SESOP_RandomMatrix_p1_1_r11
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 16 with relative residual of 0.044.*")

# test 5b: checking results
add_test (Basso_SESOP_RandomMatrix_p1_1_r11_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p1_1_r11.m solution_p1_1_r11.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p1_1_r11_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p1_1_r11_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p1_1_r11)

# test 6a: 
add_test (Basso_SESOP_RandomMatrix_p1_1_rinf
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
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
set_tests_properties (Basso_SESOP_RandomMatrix_p1_1_rinf
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 73 with relative residual of 0.0091.*")

# test 6b: checking results
add_test (Basso_SESOP_RandomMatrix_p1_1_rinf_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p1_1_rinf.m solution_p1_1_rinf.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p1_1_rinf_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p1_1_rinf_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p1_1_rinf)

# test 7a: 
add_test (Basso_SESOP_RandomMatrix_p6_r2
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
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
set_tests_properties (Basso_SESOP_RandomMatrix_p6_r2
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 60 with relative residual of 0.0094.*")

# test 7b: checking results
add_test (Basso_SESOP_RandomMatrix_p6_r2_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p6_r2.m solution_p6_r2.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p6_r2_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p6_r2_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p6_r2)

# test 8a: 
add_test (Basso_SESOP_RandomMatrix_p6_r11
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
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
set_tests_properties (Basso_SESOP_RandomMatrix_p6_r11
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 31 with relative residual of 0.051.*")

# test 8b: checking results
add_test (Basso_SESOP_RandomMatrix_p6_r11_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p6_r11.m solution_p6_r11.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p6_r11_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p6_r11_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p6_r11)

# test 9a: 
add_test (Basso_SESOP_RandomMatrix_p6_rinf
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
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
set_tests_properties (Basso_SESOP_RandomMatrix_p6_rinf
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 129 with relative residual of 0.0084.*")

# test 9b: checking results
add_test (Basso_SESOP_RandomMatrix_p6_rinf_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p6_rinf.m solution_p6_rinf.m)
set_tests_properties (Basso_SESOP_RandomMatrix_p6_rinf_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_p6_rinf_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_p6_rinf)
