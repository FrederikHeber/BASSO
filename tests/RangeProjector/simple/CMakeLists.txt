# does simple call to configurator work?
add_test (RangeProjectorConfiguratorSimple_Vector 
	${PROJECT_BINARY_DIR}/Examples/RangeProjector/RangeProjectorConfigurator
		--config RangeProjector.cfg
		--verbose 2
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 1e-8 
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/matrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/vector.m
		--solution solution.mat)
set_tests_properties (RangeProjectorConfiguratorSimple_Vector
	PROPERTIES PASS_REGULAR_EXPRESSION "Magnitude of noise.*.1")

# test 2b: checking results
add_test (RangeProjectorConfiguratorSimple_Vector_check
	diff
		--ignore-matching-lines=matrix
		--ignore-matching-lines=rhs
		RangeProjector.cfg 
		${CMAKE_CURRENT_SOURCE_DIR}/post/RangeProjector.cfg)
set_tests_properties (RangeProjectorConfiguratorSimple_Vector_check
	PROPERTIES DEPENDS RangeProjectorConfiguratorSimple_Vector)

# does simple test call work?
add_test (RangeProjectorSimple 
	${PROJECT_BINARY_DIR}/Examples/RangeProjector/RangeProjector
		--config RangeProjector.cfg)
set_tests_properties (RangeProjectorSimple
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of X.*2")
set_tests_properties (RangeProjectorSimple
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of Y.*2")
set_tests_properties (RangeProjectorSimple
	PROPERTIES PASS_REGULAR_EXPRESSION ".*noise.*0.01")
set_tests_properties (RangeProjectorSimple
	PROPERTIES PASS_REGULAR_EXPRESSION ".*Ax = y.*1 0")
set_tests_properties (RangeProjectorSimple
	PROPERTIES DEPENDS RangeProjectorConfiguratorSimple_Vector)

# test 2: project twice does not change it anymore
add_test (RangeProjectorSecondProjection
	${PROJECT_BINARY_DIR}/Examples/RangeProjector/RangeProjector 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 1e-8 
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/matrix.m 
		--rhs solution.mat
		--solution solution2.mat)
set_tests_properties (RangeProjectorSecondProjection
	PROPERTIES DEPENDS RangeProjectorSimple)

# test 2b: checking results
add_test (RangeProjectorSecondProjection_check
	${DIFFCMD} ${NDIFF_OPTIONS} solution.mat solution2.mat)
set_tests_properties (RangeProjectorSecondProjection_check
	PROPERTIES DEPENDS RangeProjectorSecondProjection)

# does simple test call fail when matrix does not exist?
add_test (RangeProjectorSimple_NonExistantMatrix 
	${PROJECT_BINARY_DIR}/Examples/RangeProjector/RangeProjector 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01 
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/nomatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/vector.m)
set_tests_properties (RangeProjectorSimple_NonExistantMatrix
	PROPERTIES WILL_FAIL true)
