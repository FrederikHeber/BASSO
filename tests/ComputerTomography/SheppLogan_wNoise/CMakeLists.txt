# does simple test call work?
add_test (ComputerTomographySheppLogan_wNoise 
	${PROJECT_BINARY_DIR}/Examples/ComputerTomography/ComputerTomography
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01 
		--maxiter 50
		--noise-level 0.01
		--num-angles 40
		--num-offsets 41
		--num-pixels-x 41
		--num-pixels-y 41
		--solution solution.mat
		--phantom ${PROJECT_SOURCE_DIR}/tests/ComputerTomography/SheppLogan_wNoise/pre/phantom.mat)
set_tests_properties (ComputerTomographySheppLogan_wNoise
	PROPERTIES PASS_REGULAR_EXPRESSION "olution after 50.* residual of .*0.03")

# test 1b: checking results
add_test (ComputerTomographySheppLogan_wNoise_check
	diff ${PROJECT_SOURCE_DIR}/tests/ComputerTomography/SheppLogan_wNoise/post/solution.mat solution.mat)
set_tests_properties (ComputerTomographySheppLogan_wNoise_check
	PROPERTIES DEPENDS ComputerTomographySheppLogan_wNoise)
