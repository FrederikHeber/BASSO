# test 1: 
add_test (Basso_SSO_RandomMatrix
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm SSO
		--number-directions 2
		--update-algorithm 0
		--normx 2 
		--normy 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--matrix ${PROJECT_SOURCE_DIR}/tests/SSO/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/SSO/pre/randomvector.m)
set_tests_properties (Basso_SSO_RandomMatrix
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 18 with relative residual of 0.009.*")

# test 2: check angles
add_test (Basso_SSO_RandomMatrix_Angles
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm SSO
		--number-directions 4
		--verbose 1
		--normx 1.5
		--normy 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 5
		--matrix ${PROJECT_SOURCE_DIR}/tests/SSO/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/SSO/pre/randomvector.m)
set_tests_properties (Basso_SSO_RandomMatrix_Angles
	PROPERTIES PASS_REGULAR_EXPRESSION "Bregman Angles #0 is 0.1")
set_tests_properties (Basso_SSO_RandomMatrix_Angles
	PROPERTIES PASS_REGULAR_EXPRESSION "Bregman Angles #1 is 0.1.*")
set_tests_properties (Basso_SSO_RandomMatrix_Angles
	PROPERTIES PASS_REGULAR_EXPRESSION "Bregman Angles #2 is 0.04.*")
set_tests_properties (Basso_SSO_RandomMatrix_Angles
	PROPERTIES PASS_REGULAR_EXPRESSION "Bregman Angles #3 is 0.05.*")
set_tests_properties (Basso_SSO_RandomMatrix_Angles
	PROPERTIES PASS_REGULAR_EXPRESSION "Angles #0 is 0.22.*")
set_tests_properties (Basso_SSO_RandomMatrix_Angles
	PROPERTIES PASS_REGULAR_EXPRESSION "Angles #1 is 0.38.*")
set_tests_properties (Basso_SSO_RandomMatrix_Angles
	PROPERTIES PASS_REGULAR_EXPRESSION "Angles #2 is 0.16.*")
set_tests_properties (Basso_SSO_RandomMatrix_Angles
	PROPERTIES PASS_REGULAR_EXPRESSION "Angles #3 is 0.23.*")

# test 3: updateToMostParallel
add_test (Basso_SSO_RandomMatrix_UpdateToMostParallel
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm SSO
		--number-directions 2
		--update-algorithm 1
		--normx 2 
		--normy 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--matrix ${PROJECT_SOURCE_DIR}/tests/SSO/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/SSO/pre/randomvector.m)
set_tests_properties (Basso_SSO_RandomMatrix_UpdateToMostParallel
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 24 with relative residual of 0.009.*")
