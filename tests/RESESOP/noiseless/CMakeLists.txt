# test 1: n=100, p=1.1, r=2
add_test (Basso_RESESOP_RandomMatrix_p1_1_r2
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm RESESOP
		--number-directions 2
		--normx 1.1
		--normy 2
		--powerx 2
		--powery 2
		--tau 1.1
		--delta 0.025
		--maxiter 500 
		--matrix ${PROJECT_SOURCE_DIR}/tests/RESESOP/noiseless/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/RESESOP/noiseless/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/RESESOP/noiseless/pre/solution.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_r2
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 42 with relative residual of 0.013.*")

# test 2: n=100, p=1.1, r=11
add_test (Basso_RESESOP_RandomMatrix_p1_1_r11
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm RESESOP
		--number-directions 2
		--normx 1.1
		--normy 11
		--powerx 2
		--powery 11
		--tau 1.1
		--delta 0.025
		--maxiter 500 
		--matrix ${PROJECT_SOURCE_DIR}/tests/RESESOP/noiseless/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/RESESOP/noiseless/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/RESESOP/noiseless/pre/solution.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_r11
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 17 with relative residual of 0.05.*")

# test 3: n=100, p=1.1, r=inf
add_test (Basso_RESESOP_RandomMatrix_p1_1_rinf
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm RESESOP
		--number-directions 2
		--normx 1.1
		--normy 0 
		--powerx 2
		--powery 2
		--tau 1.1
		--delta 0.025
		--maxiter 500 
		--matrix ${PROJECT_SOURCE_DIR}/tests/RESESOP/noiseless/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/RESESOP/noiseless/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/RESESOP/noiseless/pre/solution.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p1_1_rinf
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 20 with relative residual of 0.072.*")

# test 4: n=100, p=2, r=2
add_test (Basso_RESESOP_RandomMatrix_p2_r2
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm RESESOP
		--number-directions 2
		--normx 2 
		--normy 2
		--powerx 2
		--powery 2
		--tau 1.1
		--delta 0.025
		--maxiter 500 
		--matrix ${PROJECT_SOURCE_DIR}/tests/RESESOP/noiseless/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/RESESOP/noiseless/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/RESESOP/noiseless/pre/solution.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_r2
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 100 with relative residual of 0.013.*")

# test 5: n=100, p=2, r=11
add_test (Basso_RESESOP_RandomMatrix_p2_r11
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm RESESOP
		--number-directions 2
		--normx 2 
		--normy 11
		--powerx 2
		--powery 11
		--tau 1.1
		--delta 0.025
		--maxiter 500 
		--matrix ${PROJECT_SOURCE_DIR}/tests/RESESOP/noiseless/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/RESESOP/noiseless/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/RESESOP/noiseless/pre/solution.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_r11
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 56 with relative residual of 0.056.*")

# test 6: n=100, p=2, r=inf
add_test (Basso_RESESOP_RandomMatrix_p2_rinf
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm RESESOP
		--number-directions 2
		--normx 2 
		--normy 0 
		--powerx 2
		--powery 2
		--tau 1.1
		--delta 0.025
		--maxiter 500 
		--matrix ${PROJECT_SOURCE_DIR}/tests/RESESOP/noiseless/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/RESESOP/noiseless/pre/randomvector.m
		--compare-against ${PROJECT_SOURCE_DIR}/tests/RESESOP/noiseless/pre/solution.m)
set_tests_properties (Basso_RESESOP_RandomMatrix_p2_rinf
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 44 with relative residual of 0.076.*")

