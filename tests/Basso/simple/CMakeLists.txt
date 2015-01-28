# does simple test call work?
add_test (BassoSimple 
	${PROJECT_BINARY_DIR}/src/Basso 
		--normx 2 
		--normy 2 
		--powery 2 
		--delta 0.01 
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/simple/pre/matrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/simple/pre/vector.m)
set_tests_properties (BassoSimple
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of X.*2")
set_tests_properties (BassoSimple
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of Y.*2")
set_tests_properties (BassoSimple
	PROPERTIES PASS_REGULAR_EXPRESSION ".*noise.*0.01")
set_tests_properties (BassoSimple
	PROPERTIES PASS_REGULAR_EXPRESSION ".*Ax = y.*1 0")
