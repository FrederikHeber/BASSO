# test 1: check angles
add_test (Basso_SESOP_RandomMatrix_Angles
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 4
		--calculateAngles 1
		--verbose 1
		--type-space-x "lp" 
		--px 1.5
		--type-space-y "lp" 
		--py 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 5
		--iteration-file SESOP.db
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/randommatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/randomvector.m)
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
