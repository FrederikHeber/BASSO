# test 1: 
add_test (Basso_SSO_RandomMatrix
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm SSO
		--number-directions 2
		--normx 2 
		--normy 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--matrix ${PROJECT_SOURCE_DIR}/tests/SSO/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/SSO/pre/randomvector.m)
set_tests_properties (Basso_SSO_RandomMatrix
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 83 with residual error of 0.01.*")
