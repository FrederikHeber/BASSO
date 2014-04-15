# test 1: 
add_test (Basso_Schoepfer2006_StraightLine1 
	${PROJECT_BINARY_DIR}/src/Basso 
		--normx 2 
		--normy 2 
		--delta 0.01 
		--C 1 
		--maxiter 50 
		--matrix ${PROJECT_SOURCE_DIR}/tests/Schoepfer2006_StraightLine/pre/matrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Schoepfer2006_StraightLine/pre/vector.m)
set_tests_properties (Basso_Schoepfer2006_StraightLine1
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution is 8.0.*e-01 4.0.*e-01")

# test 2: 
add_test (Basso_Schoepfer2006_StraightLine2 
	${PROJECT_BINARY_DIR}/src/Basso 
		--normx 4 
		--normy 4 
		--delta 0.01
		--C 6
		--maxiter 10000
		--matrix ${PROJECT_SOURCE_DIR}/tests/Schoepfer2006_StraightLine/pre/matrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Schoepfer2006_StraightLine/pre/vector.m)
set_tests_properties (Basso_Schoepfer2006_StraightLine2
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution is 7.1.*e-01 5.7.*e-01")

# test 3: 
add_test (Basso_Schoepfer2006_StraightLine3 
	${PROJECT_BINARY_DIR}/src/Basso 
		--normx 4 
		--normy 20 
		--delta 0.01
		--C 6
		--maxiter 10000
		--matrix ${PROJECT_SOURCE_DIR}/tests/Schoepfer2006_StraightLine/pre/matrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Schoepfer2006_StraightLine/pre/vector.m)
set_tests_properties (Basso_Schoepfer2006_StraightLine3
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution is 7.1.*e-01 5.7.*e-01")

# test 4: 
add_test (Basso_Schoepfer2006_StraightLine4 
	${PROJECT_BINARY_DIR}/src/Basso 
		--normx 6 
		--normy 40 
		--delta 0.01
		--C 10
		--maxiter 10000
		--matrix ${PROJECT_SOURCE_DIR}/tests/Schoepfer2006_StraightLine/pre/matrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Schoepfer2006_StraightLine/pre/vector.m)
set_tests_properties (Basso_Schoepfer2006_StraightLine4
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution is 7.0.*e-01 6.0.*e-01")

# test 5: 
add_test (Basso_Schoepfer2006_StraightLine5 
	${PROJECT_BINARY_DIR}/src/Basso 
		--normx 1.5 
		--normy 2 
		--delta 0.01
		--C 2
		--maxiter 1000
		--matrix ${PROJECT_SOURCE_DIR}/tests/Schoepfer2006_StraightLine/pre/matrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Schoepfer2006_StraightLine/pre/vector.m)
set_tests_properties (Basso_Schoepfer2006_StraightLine5
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution is 8.9.*e-01 2.2.*e-01")

# test 6: 
add_test (Basso_Schoepfer2006_StraightLine6 
	${PROJECT_BINARY_DIR}/src/Basso 
		--normx 1.5 
		--normy 10 
		--delta 0.01
		--C 2
		--maxiter 1000
		--matrix ${PROJECT_SOURCE_DIR}/tests/Schoepfer2006_StraightLine/pre/matrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Schoepfer2006_StraightLine/pre/vector.m)
set_tests_properties (Basso_Schoepfer2006_StraightLine6
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution is 8.9.*e-01 2.2.*e-01")

# test 7: 
add_test(Basso_Schoepfer2006_StraightLine7 
	${PROJECT_BINARY_DIR}/src/Basso 
		--normx 1.1 
		--normy 5 
		--delta 0.01
		--C 8
		--maxiter 5000
		--matrix ${PROJECT_SOURCE_DIR}/tests/Schoepfer2006_StraightLine/pre/matrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Schoepfer2006_StraightLine/pre/vector.m)
set_tests_properties (Basso_Schoepfer2006_StraightLine7
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution is 1.0.*e-01 0.0.*e-01")
	