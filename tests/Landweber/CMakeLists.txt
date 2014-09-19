# test 1: 
add_test (Basso_Landweber_StraightLine1
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm landweber
		--normx 2 
		--normy 2 
		--delta 0.0001 
		--C 1 
		--maxiter 50 
		--matrix ${PROJECT_SOURCE_DIR}/tests/Landweber/pre/matrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Landweber/pre/vector.m)
set_tests_properties (Basso_Landweber_StraightLine1
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution .* is 8.0.*e-01 4.0.*e-01")

# test 2: 
add_test (Basso_Landweber_StraightLine2
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm landweber
		--normx 4 
		--normy 20 
		--delta 0.01 
		--C 6 
		--maxiter 50 
		--matrix ${PROJECT_SOURCE_DIR}/tests/Landweber/pre/matrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Landweber/pre/vector.m)
set_tests_properties (Basso_Landweber_StraightLine2
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution .* is 7.1.*e-01 5.6.*e-01")

# test 3: 
add_test (Basso_Landweber_StraightLine3
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm landweber
		--normx 1.5 
		--normy 10 
		--delta 0.0001 
		--C 2 
		--maxiter 50 
		--matrix ${PROJECT_SOURCE_DIR}/tests/Landweber/pre/matrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Landweber/pre/vector.m)
set_tests_properties (Basso_Landweber_StraightLine3
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution .* is 8.8.*e-01 2.2.*e-01")

# test 4: 
add_test (Basso_Landweber_RandomMatrix
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm landweber
		--normx 2 
		--normy 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--C 1 
		--maxiter 100 
		--matrix ${PROJECT_SOURCE_DIR}/tests/Landweber/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Landweber/pre/randomvector.m)
set_tests_properties (Basso_Landweber_RandomMatrix
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 100 with residual error of 1.59.*")
