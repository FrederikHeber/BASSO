# test 1: 
add_test (Basso_SESOP_RandomMatrix
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm SESOP
		--number-directions 2
		--update-algorithm 0
		--normx 2 
		--normy 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/pre/randomvector.m)
set_tests_properties (Basso_SESOP_RandomMatrix
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 18 with relative residual of 0.009.*")

# test 2: check angles
add_test (Basso_SESOP_RandomMatrix_Angles
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm SESOP
		--number-directions 4
		--calculateAngles 1
		--verbose 1
		--normx 1.5
		--normy 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 5
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/pre/randomvector.m)
set_tests_properties (Basso_SESOP_RandomMatrix_Angles
	PROPERTIES PASS_REGULAR_EXPRESSION "Bregman Angles #0 is 0.1.*")
set_tests_properties (Basso_SESOP_RandomMatrix_Angles
	PROPERTIES PASS_REGULAR_EXPRESSION "Bregman Angles #1 is 0.15.*")
set_tests_properties (Basso_SESOP_RandomMatrix_Angles
	PROPERTIES PASS_REGULAR_EXPRESSION "Bregman Angles #2 is 0.04.*")
set_tests_properties (Basso_SESOP_RandomMatrix_Angles
	PROPERTIES PASS_REGULAR_EXPRESSION "Bregman Angles #3 is 0.05.*")
set_tests_properties (Basso_SESOP_RandomMatrix_Angles
	PROPERTIES PASS_REGULAR_EXPRESSION "Angles #0 is 0.11.*")
set_tests_properties (Basso_SESOP_RandomMatrix_Angles
	PROPERTIES PASS_REGULAR_EXPRESSION "Angles #1 is 0.15.*")
set_tests_properties (Basso_SESOP_RandomMatrix_Angles
	PROPERTIES PASS_REGULAR_EXPRESSION "Angles #2 is 0.04.*")
set_tests_properties (Basso_SESOP_RandomMatrix_Angles
	PROPERTIES PASS_REGULAR_EXPRESSION "Angles #3 is 0.05.*")

# test 3: updateToMostParallel
add_test (Basso_SESOP_RandomMatrix_UpdateToMostParallel
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm SESOP
		--number-directions 2
		--update-algorithm 1
		--enforceRandomMapping 0
		--normx 2 
		--normy 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/pre/randomvector.m)
set_tests_properties (Basso_SESOP_RandomMatrix_UpdateToMostParallel
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 19 with relative residual of 0.0089.*")

# test 4: updateToMostOrthogonal
add_test (Basso_SESOP_RandomMatrix_UpdateToMostOrthogonal
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm SESOP
		--number-directions 2
		--update-algorithm 2
		--enforceRandomMapping 0
		--normx 2 
		--normy 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/pre/randomvector.m)
set_tests_properties (Basso_SESOP_RandomMatrix_UpdateToMostOrthogonal
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 25 with relative residual of 0.0047.*")

# test 5: updateToMostParallel_randomMapping
add_test (Basso_SESOP_RandomMatrix_UpdateToMostParallel_randomMapping
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm SESOP
		--number-directions 2
		--update-algorithm 1
		--enforceRandomMapping 1
		--normx 2 
		--normy 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/pre/randomvector.m)
set_tests_properties (Basso_SESOP_RandomMatrix_UpdateToMostParallel_randomMapping
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 21 with relative residual of 0.0083.*")

# test 6: updateToMostOrthogonal_randomMapping
add_test (Basso_SESOP_RandomMatrix_UpdateToMostOrthogonal_randomMapping
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm SESOP
		--number-directions 2
		--update-algorithm 2
		--enforceRandomMapping 1
		--normx 2 
		--normy 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/pre/randomvector.m)
set_tests_properties (Basso_SESOP_RandomMatrix_UpdateToMostOrthogonal_randomMapping
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 23 with relative residual of 0.0094.*")

# test 7: 
add_test (Basso_SESOP_RandomMatrix_InexactLinesearch
	${PROJECT_BINARY_DIR}/src/Basso 
		--algorithm SESOP
		--number-directions 2
		--update-algorithm 0
		--inexact-linesearch 1
		--wolfe-constants 1e-2 0.6
		--normx 2 
		--normy 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/pre/randomvector.m)
set_tests_properties (Basso_SESOP_RandomMatrix_InexactLinesearch
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 33 with relative residual of 0.008.*")

# test 8: 
add_test (Basso_SESOP_RandomMatrix_NemirovskyDirections
	${PROJECT_BINARY_DIR}/src/Basso
		--algorithm SESOP
		--number-directions 2
		--searchspace "Nemirovsky"
		--update-algorithm 0
		--normx 2
		--normy 0
		--powerx 2
		--powery 2
		--delta 0.01
		--maxiter 100
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/pre/randommatrix.m
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/pre/randomvector.m)
set_tests_properties (Basso_SESOP_RandomMatrix_NemirovskyDirections
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 24 with relative residual of 0.0095.*")

# test 9: Nemirovsky with not two search directions
add_test (Basso_SESOP_RandomMatrix_NemirovskyDirections_fail
	${PROJECT_BINARY_DIR}/src/Basso
		--algorithm SESOP
		--number-directions 4
		--searchspace "Nemirovsky"
		--update-algorithm 0
		--normx 2
		--normy 0
		--powerx 2
		--powery 2
		--delta 0.01
		--maxiter 100
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/pre/randommatrix.m
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/pre/randomvector.m)
set_tests_properties (Basso_SESOP_RandomMatrix_NemirovskyDirections_fail
	PROPERTIES WILL_FAIL "true")

