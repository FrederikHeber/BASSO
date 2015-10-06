# fails because no measurement dimensions
add_test (ComputerTomographySimple_failmeasurementdimensions 
	${PROJECT_BINARY_DIR}/src/ComputerTomographyBase/ComputerTomography 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01 
		--rhs ${PROJECT_SOURCE_DIR}/tests/ComputerTomography/simple/pre/sinogram.mat)
set_tests_properties (ComputerTomographySimple_failmeasurementdimensions
	PROPERTIES WILL_FAIL true)

# fails because of missing rhs file 	
add_test (ComputerTomographySimple_failmissingfile 
	${PROJECT_BINARY_DIR}/src/ComputerTomographyBase/ComputerTomography 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01 
		--num-angles 60
		--num-offsets 61
		--num-pixels-x 41
		--num-pixels-y 41)
set_tests_properties (ComputerTomographySimple_failmissingfile
	PROPERTIES WILL_FAIL true)

# succeeds 	
add_test (ComputerTomographySimple_succeeds 
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
set_tests_properties (ComputerTomographySimple_succeeds
	PROPERTIES PASS_REGULAR_EXPRESSION "olution after 10 with .* residual of 0.23")
