# test 1a: n=100, p=1.1, r=2
add_test (Basso_RESESOP_RandomMatrix_p1_1_r2
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm RESESOP
		--number-directions 2
		--type-space-x "lp" 
		--px 1.1
		--type-space-y "lp" 
		--py 2
		--powerx 2
		--powery 2
		--tau 1.1
		--delta 0.025
		--maxiter 500 
		--iteration-file RESESOP_p1_1_r2.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--compare-against ${CMAKE_CURRENT_SOURCE_DIR}/pre/solution.m
		--solution solution_p1_1_r2.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_r2
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 42 with relative residual of 0.013.*")

# test 1b: checking results
add_test (Basso_RESESOP_RandomMatrix_p1_1_r2_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p1_1_r2.m  solution_p1_1_r2.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_r2_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_r2_check
	PROPERTIES DEPENDS Basso_RESESOP_RandomMatrix_p1_1_r2 )

# test 2a: n=100, p=1.1, r=11
add_test (Basso_RESESOP_RandomMatrix_p1_1_r11
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm RESESOP
		--number-directions 2
		--type-space-x "lp" 
		--px 1.1
		--type-space-y "lp" 
		--py 11
		--powerx 2
		--powery 11
		--tau 1.1
		--delta 0.025
		--maxiter 500 
		--iteration-file RESESOP_p1_1_r11.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--compare-against ${CMAKE_CURRENT_SOURCE_DIR}/pre/solution.m
		--solution solution_p1_1_r11.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_r11
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 17 with relative residual of 0.05.*")

# test 2b: checking results
add_test (Basso_RESESOP_RandomMatrix_p1_1_r11_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p1_1_r11.m  solution_p1_1_r11.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_r11_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_r11_check
	PROPERTIES DEPENDS Basso_RESESOP_RandomMatrix_p1_1_r11 )

# test 3a: n=100, p=1.1, r=inf
add_test (Basso_RESESOP_RandomMatrix_p1_1_rinf
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm RESESOP
		--number-directions 2
		--type-space-x "lp" 
		--px 1.1
		--type-space-y "lp" 
		--py 0 
		--powerx 2
		--powery 2
		--tau 1.1
		--delta 0.025
		--maxiter 500 
		--iteration-file RESESOP_p1_1_rinf.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--compare-against ${CMAKE_CURRENT_SOURCE_DIR}/pre/solution.m
		--solution solution_p1_1_rinf.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_rinf
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 20 with relative residual of 0.072.*")

# test 3b: checking results
add_test (Basso_RESESOP_RandomMatrix_p1_1_rinf_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p1_1_rinf.m  solution_p1_1_rinf.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_rinf_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_rinf_check
	PROPERTIES DEPENDS Basso_RESESOP_RandomMatrix_p1_1_rinf )

# test 4a: n=100, p=2, r=2
add_test (Basso_RESESOP_RandomMatrix_p2_r2
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm RESESOP
		--number-directions 2
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2
		--powerx 2
		--powery 2
		--tau 1.1
		--delta 0.025
		--maxiter 500 
		--iteration-file RESESOP_p2_r2.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--compare-against ${CMAKE_CURRENT_SOURCE_DIR}/pre/solution.m
		--solution solution_p2_r2.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_r2
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 94 with relative residual of 0.013.*")

# test 4b: checking results
add_test (Basso_RESESOP_RandomMatrix_p2_r2_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p2_r2.m  solution_p2_r2.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_r2_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_r2_check
	PROPERTIES DEPENDS Basso_RESESOP_RandomMatrix_p2_r2 )

# test 5a: n=100, p=2, r=11
add_test (Basso_RESESOP_RandomMatrix_p2_r11
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm RESESOP
		--number-directions 2
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 11
		--powerx 2
		--powery 11
		--tau 1.1
		--delta 0.025
		--maxiter 500 
		--iteration-file RESESOP_p2_r11.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--compare-against ${CMAKE_CURRENT_SOURCE_DIR}/pre/solution.m
		--solution solution_p2_r11.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_r11
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 56 with relative residual of 0.056.*")

# test 5b: checking results
add_test (Basso_RESESOP_RandomMatrix_p2_r11_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p2_r11.m  solution_p2_r11.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_r11_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_r11_check
	PROPERTIES DEPENDS Basso_RESESOP_RandomMatrix_p2_r11 )

# test 6a: n=100, p=2, r=inf
add_test (Basso_RESESOP_RandomMatrix_p2_rinf
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm RESESOP
		--number-directions 2
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 0 
		--powerx 2
		--powery 2
		--tau 1.1
		--delta 0.025
		--maxiter 500 
		--iteration-file RESESOP_p2_rinf.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--compare-against ${CMAKE_CURRENT_SOURCE_DIR}/pre/solution.m
		--solution solution_p2_rinf.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_rinf
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 44 with relative residual of 0.076.*")

# test 6b: checking results
add_test (Basso_RESESOP_RandomMatrix_p2_rinf_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p2_rinf.m  solution_p2_rinf.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_rinf_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_rinf_check
	PROPERTIES DEPENDS Basso_RESESOP_RandomMatrix_p2_rinf )


# test 7a: n=100, p=6, r=2
add_test (Basso_RESESOP_RandomMatrix_p6_r2
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm RESESOP
		--number-directions 2
		--type-space-x "lp" 
		--px 6
		--type-space-y "lp" 
		--py 2
		--powerx 6
		--powery 2
		--tau 1.1
		--delta 0.025
		--maxiter 500 
		--iteration-file RESESOP_p6_r2.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--compare-against ${CMAKE_CURRENT_SOURCE_DIR}/pre/solution.m
		--solution solution_p6_r2.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p6_r2
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 500 with relative residual of 0.163.*")

# test 7b: checking results
add_test (Basso_RESESOP_RandomMatrix_p6_r2_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p6_r2.m  solution_p6_r2.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p6_r2_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_RESESOP_RandomMatrix_p6_r2_check
	PROPERTIES DEPENDS Basso_RESESOP_RandomMatrix_p6_r2)

# test 8a: n=100, p=6, r=11
add_test (Basso_RESESOP_RandomMatrix_p6_r11
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm RESESOP
		--number-directions 2
		--type-space-x "lp" 
		--px 6
		--type-space-y "lp" 
		--py 11
		--powerx 6
		--powery 11
		--tau 1.1
		--delta 0.025
		--maxiter 500 
		--iteration-file RESESOP_p6_r11.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--compare-against ${CMAKE_CURRENT_SOURCE_DIR}/pre/solution.m
		--solution solution_p6_r11.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p6_r11
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 500 with relative residual of 0.30.*")

# test 8b: checking results
add_test (Basso_RESESOP_RandomMatrix_p6_r11_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p6_r11.m  solution_p6_r11.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p6_r11_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_RESESOP_RandomMatrix_p6_r11_check
	PROPERTIES DEPENDS Basso_RESESOP_RandomMatrix_p6_r11)

# test 9a: n=100, p=6, r=inf
add_test (Basso_RESESOP_RandomMatrix_p6_rinf
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm RESESOP
		--number-directions 2
		--type-space-x "lp" 
		--px 6
		--type-space-y "lp" 
		--py 0 
		--powerx 6
		--powery 2
		--tau 1.1
		--delta 0.025
		--maxiter 500 
		--iteration-file RESESOP_p6_rinf.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m
		--compare-against ${CMAKE_CURRENT_SOURCE_DIR}/pre/solution.m
		--solution solution_p6_rinf.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p6_rinf
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 46 with relative residual of 0.074.*")

# test 9b: checking results
add_test (Basso_RESESOP_RandomMatrix_p6_rinf_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution_p6_rinf.m  solution_p6_rinf.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p6_rinf_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_RESESOP_RandomMatrix_p6_rinf_check
	PROPERTIES DEPENDS Basso_RESESOP_RandomMatrix_p6_rinf)
	