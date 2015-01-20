# test 1: n=100, p=1.1, r=2
add_test (Basso_RESESOP_RandomMatrix_p1_1_r2_regularized
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm RESESOP
		--number-directions 2
		--normx 1.1
		--normy 2
		--powerx 2
		--powery 2
		--tau 1.1
		--delta 0.1496300733072937
		--maxiter 500 
		--matrix ${PROJECT_SOURCE_DIR}/tests/RESESOP/noisy/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/RESESOP/noisy/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/RESESOP/noisy/pre/solution.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_r2_regularized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 23 with relative residual of 0.082.*")

# test 2: n=100, p=1.1, r=11
add_test (Basso_RESESOP_RandomMatrix_p1_1_r11_regularized
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm RESESOP
		--number-directions 2
		--normx 1.1
		--normy 11
		--powerx 2
		--powery 11
		--tau 1.1
		--delta 0.02967849813273046
		--maxiter 500 
		--matrix ${PROJECT_SOURCE_DIR}/tests/RESESOP/noisy/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/RESESOP/noisy/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/RESESOP/noisy/pre/solution.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_r11_regularized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 30 with relative residual of 0.074.*")

# test 3: n=100, p=1.1, r=inf
add_test (Basso_RESESOP_RandomMatrix_p1_1_rinf_regularized
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm RESESOP
		--number-directions 2
		--normx 1.1
		--normy 0 
		--powerx 2
		--powery 2
		--tau 1.1
		--delta 0.02495206594467164
		--maxiter 500 
		--matrix ${PROJECT_SOURCE_DIR}/tests/RESESOP/noisy/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/RESESOP/noisy/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/RESESOP/noisy/pre/solution.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_rinf_regularized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 21 with relative residual of 0.074.*")
	
# test 4: n=100, p=2, r=2
add_test (Basso_RESESOP_RandomMatrix_p2_r2_regularized
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm RESESOP
		--number-directions 2
		--normx 2 
		--normy 2
		--powerx 2
		--powery 2
		--tau 1.1
		--delta 0.1496300733072937
		--maxiter 500 
		--matrix ${PROJECT_SOURCE_DIR}/tests/RESESOP/noisy/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/RESESOP/noisy/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/RESESOP/noisy/pre/solution.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_r2_regularized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 29 with relative residual of 0.082.*")

# test 5: n=100, p=2, r=11
add_test (Basso_RESESOP_RandomMatrix_p2_r11_regularized
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm RESESOP
		--number-directions 2
		--normx 2 
		--normy 11
		--powerx 2
		--powery 11
		--tau 1.1
		--delta 0.02967849813273046
		--maxiter 500 
		--matrix ${PROJECT_SOURCE_DIR}/tests/RESESOP/noisy/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/RESESOP/noisy/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/RESESOP/noisy/pre/solution.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_r11_regularized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 50 with relative residual of 0.074.*")

# test 6: n=100, p=2, r=inf
add_test (Basso_RESESOP_RandomMatrix_p2_rinf_regularized
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm RESESOP
		--number-directions 2
		--normx 2 
		--normy 0 
		--powerx 2
		--powery 2
		--tau 1.1
		--delta 0.02495206594467164
		--maxiter 500 
		--matrix ${PROJECT_SOURCE_DIR}/tests/RESESOP/noisy/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/RESESOP/noisy/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/RESESOP/noisy/pre/solution.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_rinf_regularized
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 59 with relative residual of 0.077.*")

