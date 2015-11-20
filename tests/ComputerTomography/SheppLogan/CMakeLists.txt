# does simple test call work?
add_test (ComputerTomographySheppLogan 
	${PROJECT_BINARY_DIR}/Examples/ComputerTomography/ComputerTomography
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01 
		--maxiter 50
		--num-angles 60
		--num-offsets 61
		--num-pixels-x 21
		--num-pixels-y 21
		--solution solution.mat
		--sinogram ${CMAKE_CURRENT_SOURCE_DIR}/pre/sinogram.mat)
set_tests_properties (ComputerTomographySheppLogan
	PROPERTIES PASS_REGULAR_EXPRESSION "olution after 50.* residual of .*0.29")

# test 1b: checking results
add_test (ComputerTomographySheppLogan_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/solution.mat solution.mat)
set_tests_properties (ComputerTomographySheppLogan_check
	PROPERTIES DEPENDS ComputerTomographySheppLogan)
