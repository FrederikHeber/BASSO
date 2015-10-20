# test 1a: n=100, p=1.1, r=2
add_test (Basso_RESESOP_RandomMatrix_p1_1_r2_regularized
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
		--delta 0.1496300733072937
		--maxiter 500 
		--iteration-file RESESOP_p1_1_r2.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/solution.m
		--solution solution_p1_1_r2.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_r2_regularized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 23 with relative residual of 0.082.*")

# test 1b: checking results
add_test (Basso_RESESOP_RandomMatrix_p1_1_r2_regularized_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/post/solution_p1_1_r2.m  solution_p1_1_r2.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_r2_regularized_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_r2_regularized_check
	PROPERTIES DEPENDS Basso_RESESOP_RandomMatrix_p1_1_r2_regularized)

# test 2a: n=100, p=1.1, r=11
add_test (Basso_RESESOP_RandomMatrix_p1_1_r11_regularized
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
		--delta 0.02967849813273046
		--maxiter 500 
		--iteration-file RESESOP_p1_1_r11.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/solution.m
		--solution solution_p1_1_r11.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_r11_regularized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 30 with relative residual of 0.074.*")

# test 2b: checking results
add_test (Basso_RESESOP_RandomMatrix_p1_1_r11_regularized_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/post/solution_p1_1_r11.m  solution_p1_1_r11.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_r11_regularized_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_r11_regularized_check
	PROPERTIES DEPENDS Basso_RESESOP_RandomMatrix_p1_1_r11_regularized)

# test 3a: n=100, p=1.1, r=inf
add_test (Basso_RESESOP_RandomMatrix_p1_1_rinf_regularized
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
		--delta 0.02495206594467164
		--maxiter 500 
		--iteration-file RESESOP_p1_1_rinf.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/solution.m
		--solution solution_p1_1_rinf.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_rinf_regularized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 21 with relative residual of 0.074.*")

# test 3b: checking results
add_test (Basso_RESESOP_RandomMatrix_p1_1_rinf_regularized_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/post/solution_p1_1_rinf.m  solution_p1_1_rinf.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_rinf_regularized_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_rinf_regularized_check
	PROPERTIES DEPENDS Basso_RESESOP_RandomMatrix_p1_1_rinf_regularized)
	
# test 4a: n=100, p=2, r=2
add_test (Basso_RESESOP_RandomMatrix_p2_r2_regularized
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
		--delta 0.1496300733072937
		--maxiter 500 
		--iteration-file RESESOP_p2_r2.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/solution.m
		--solution solution_p2_r2.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_r2_regularized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 29 with relative residual of 0.082.*")

# test 4b: checking results
add_test (Basso_RESESOP_RandomMatrix_p2_r2_regularized_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/post/solution_p2_r2.m  solution_p2_r2.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_r2_regularized_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_r2_regularized_check
	PROPERTIES DEPENDS Basso_RESESOP_RandomMatrix_p2_r2_regularized)

# test 5a: n=100, p=2, r=11
add_test (Basso_RESESOP_RandomMatrix_p2_r11_regularized
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
		--delta 0.02967849813273046
		--maxiter 500 
		--iteration-file RESESOP_p2_r11.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/solution.m
		--solution solution_p2_r11.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_r11_regularized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 50 with relative residual of 0.074.*")

# test 5b: checking results
add_test (Basso_RESESOP_RandomMatrix_p2_r11_regularized_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/post/solution_p2_r11.m  solution_p2_r11.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_r11_regularized_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_r11_regularized_check
	PROPERTIES DEPENDS Basso_RESESOP_RandomMatrix_p2_r11_regularized)

# test 6a: n=100, p=2, r=inf
add_test (Basso_RESESOP_RandomMatrix_p2_rinf_regularized
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
		--delta 0.02495206594467164
		--maxiter 500 
		--iteration-file RESESOP_p2_rinf.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/solution.m
		--solution solution_p2_rinf.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_rinf_regularized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 59 with relative residual of 0.077.*")

# test 6b: checking results
add_test (Basso_RESESOP_RandomMatrix_p2_rinf_regularized_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/post/solution_p2_rinf.m  solution_p2_rinf.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_rinf_regularized_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_rinf_regularized_check
	PROPERTIES DEPENDS Basso_RESESOP_RandomMatrix_p2_rinf_regularized)

# test 7a: n=100, p=6, r=2
add_test (Basso_RESESOP_RandomMatrix_p6_r2_regularized
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
		--delta 0.1496300733072937
		--maxiter 500 
		--iteration-file RESESOP_p6_r2.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/solution.m
		--solution solution_p6_r2.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p6_r2_regularized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 500 with relative residual of 0.089.*")

# test 7b: checking results
add_test (Basso_RESESOP_RandomMatrix_p6_r2_regularized_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/post/solution_p6_r2.m  solution_p6_r2.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p6_r2_regularized_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_RESESOP_RandomMatrix_p6_r2_regularized_check
	PROPERTIES DEPENDS Basso_RESESOP_RandomMatrix_p6_r2_regularized)

# test 8a: n=100, p=6, r=11
add_test (Basso_RESESOP_RandomMatrix_p6_r11_regularized
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
		--delta 0.02967849813273046
		--maxiter 500 
		--iteration-file RESESOP_p6_r11.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/solution.m
		--solution solution_p6_r11.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p6_r11_regularized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 500 with relative residual of 0.14.*")

# test 8b: checking results
add_test (Basso_RESESOP_RandomMatrix_p6_r11_regularized_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/post/solution_p6_r11.m  solution_p6_r11.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p6_r11_regularized_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_RESESOP_RandomMatrix_p6_r11_regularized_check
	PROPERTIES DEPENDS Basso_RESESOP_RandomMatrix_p6_r11_regularized)

# test 9a: n=100, p=6, r=inf
add_test (Basso_RESESOP_RandomMatrix_p6_rinf_regularized
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
		--delta 0.02495206594467164
		--maxiter 500 
		--iteration-file RESESOP_p6_rinf.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/pre/solution.m
		--solution solution_p6_rinf.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p6_rinf_regularized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 500 with relative residual of 0.23.*")

# test 9b: checking results
add_test (Basso_RESESOP_RandomMatrix_p6_rinf_regularized_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${PROJECT_SOURCE_DIR}/tests/Basso/RESESOP/noisy/post/solution_p6_rinf.m  solution_p6_rinf.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p6_rinf_regularized_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_RESESOP_RandomMatrix_p6_rinf_regularized_check
	PROPERTIES DEPENDS Basso_RESESOP_RandomMatrix_p6_rinf_regularized)
	
