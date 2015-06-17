# test 1: check angles
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
		--iteration-file SESOP.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/Angles/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/Angles/pre/randomvector.m)
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
