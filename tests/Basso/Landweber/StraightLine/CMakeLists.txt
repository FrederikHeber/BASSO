# test 1a: 
add_test (Basso_Landweber_StraightLine1
	${PROJECT_BINARY_DIR}/src/BassoBase/Basso 
		--algorithm Landweber
		--normx 2 
		--normy 2 
		--delta 0.0001 
		--C 1 
		--maxiter 50 
		--iteration-file Landweber1.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/Landweber/StraightLine/pre/matrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/Landweber/StraightLine/pre/vector.m
		--solution straightline1.m)
set_tests_properties (Basso_Landweber_StraightLine1
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution .* is 8.0.*e-01 4.0.*e-01")

# test 1b: checking results
add_test (Basso_Landweber_StraightLine1_check
	diff ${PROJECT_SOURCE_DIR}/tests/Basso/Landweber/StraightLine/post/straightline1.m  straightline1.m)
set_tests_properties (Basso_Landweber_StraightLine1_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_Landweber_StraightLine1_check
	PROPERTIES DEPENDS Basso_Landweber_StraightLine1 )

# test 2: 
add_test (Basso_Landweber_StraightLine2
	${PROJECT_BINARY_DIR}/src/BassoBase/Basso 
		--algorithm Landweber
		--normx 4 
		--normy 20 
		--delta 0.01 
		--C 6 
		--maxiter 50 
		--iteration-file Landweber2.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/Landweber/StraightLine/pre/matrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/Landweber/StraightLine/pre/vector.m
		--solution straightline2.m)
set_tests_properties (Basso_Landweber_StraightLine2
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution .* is 7.1.*e-01 5.6.*e-01")

# test 2b: checking results
add_test (Basso_Landweber_StraightLine2_check
	diff ${PROJECT_SOURCE_DIR}/tests/Basso/Landweber/StraightLine/post/straightline2.m  straightline2.m)
set_tests_properties (Basso_Landweber_StraightLine2_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_Landweber_StraightLine2_check
	PROPERTIES DEPENDS Basso_Landweber_StraightLine2 )

# test 3: 
add_test (Basso_Landweber_StraightLine3
	${PROJECT_BINARY_DIR}/src/BassoBase/Basso 
		--algorithm Landweber
		--normx 1.5 
		--normy 10 
		--delta 0.0001 
		--C 2 
		--maxiter 50 
		--iteration-file Landweber3.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/Landweber/StraightLine/pre/matrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/Landweber/StraightLine/pre/vector.m
		--solution straightline3.m)
set_tests_properties (Basso_Landweber_StraightLine3
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution .* is 8.8.*e-01 2.2.*e-01")

# test 3b: checking results
add_test (Basso_Landweber_StraightLine3_check
	diff ${PROJECT_SOURCE_DIR}/tests/Basso/Landweber/StraightLine/post/straightline3.m  straightline3.m)
set_tests_properties (Basso_Landweber_StraightLine3_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_Landweber_StraightLine3_check
	PROPERTIES DEPENDS Basso_Landweber_StraightLine3 )

	