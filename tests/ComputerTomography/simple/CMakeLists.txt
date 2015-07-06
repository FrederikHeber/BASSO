# does simple test call work?
add_test (ComputerTomographySimple 
	${PROJECT_BINARY_DIR}/src/ComputerTomographyBase/ComputerTomography 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01 
		--matrix ${PROJECT_SOURCE_DIR}/tests/ComputerTomography/simple/pre/matrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/ComputerTomography/simple/pre/vector.m)
set_tests_properties (ComputerTomographySimple
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of X.*2")
set_tests_properties (ComputerTomographySimple
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of Y.*2")
set_tests_properties (ComputerTomographySimple
	PROPERTIES PASS_REGULAR_EXPRESSION ".*noise.*0.01")
set_tests_properties (ComputerTomographySimple
	PROPERTIES PASS_REGULAR_EXPRESSION ".*Ax = y.*1 0")
