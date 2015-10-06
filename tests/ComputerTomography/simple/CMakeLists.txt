# does simple test call work?
add_test (ComputerTomographySimple 
	${PROJECT_BINARY_DIR}/src/ComputerTomographyBase/ComputerTomography 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01
		--maxiter 10
		--num-angles 2
		--num-offsets 3
		--num-pixels-x 11
		--num-pixels-y 11		 
		--rhs ${PROJECT_SOURCE_DIR}/tests/ComputerTomography/simple/pre/vector.m)
set_tests_properties (ComputerTomographySimple
	PROPERTIES PASS_REGULAR_EXPRESSION "olution after 10 with .* residual of 0.23")
