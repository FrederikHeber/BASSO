# does simple call to configurator work?
add_test (NoiseAdderConfiguratorSimple_Vector 
	${PROJECT_BINARY_DIR}/Examples/NoiseAdder/NoiseAdderConfigurator
		--config NoiseAdder.cfg
		--input ${PROJECT_SOURCE_DIR}/tests/NoiseAdder/simple/pre/vector.m
		--output solution_vector.mat
		--noise-level 0.1
		--relative-level 1)
set_tests_properties (NoiseAdderConfiguratorSimple_Vector
	PROPERTIES PASS_REGULAR_EXPRESSION "Noise level.*.1")

# test 2b: checking results
add_test (NoiseAdderConfiguratorSimple_Vector_check
	diff
		--ignore-matching-lines=input 
		NoiseAdder.cfg 
		${PROJECT_SOURCE_DIR}/tests/NoiseAdder/simple/post/NoiseAdder.cfg)
set_tests_properties (NoiseAdderConfiguratorSimple_Vector_check
	PROPERTIES DEPENDS NoiseAdderConfiguratorSimple_Vector)

# does simple test call work?
add_test (NoiseAdderSimple_Vector 
	${PROJECT_BINARY_DIR}/Examples/NoiseAdder/NoiseAdder
		--config NoiseAdder.cfg)
set_tests_properties (NoiseAdderSimple_Vector
	PROPERTIES PASS_REGULAR_EXPRESSION "Noise level.*.1")
set_tests_properties (NoiseAdderSimple_Vector
	PROPERTIES DEPENDS NoiseAdderConfiguratorSimple_Vector)

# does simple test call work?
add_test (NoiseAdderSimple_Matrix
	${PROJECT_BINARY_DIR}/Examples/NoiseAdder/NoiseAdder 
		--input ${PROJECT_SOURCE_DIR}/tests/NoiseAdder/simple/pre/matrix.m
		--output solution_matrix.mat
		--noise-level 0.2)
set_tests_properties (NoiseAdderSimple_Matrix
	PROPERTIES PASS_REGULAR_EXPRESSION "Noise level.*.1")

# does simple test call fail when vector does not exist?
add_test (NoiseAdderSimple_NonExistantVector 
	${PROJECT_BINARY_DIR}/Examples/NoiseAdder/NoiseAdder 
		--input ${PROJECT_SOURCE_DIR}/tests/NoiseAdder/simple/pre/noinput.m 
		--output solution.mat
		--noise-level 0.1)
set_tests_properties (NoiseAdderSimple_NonExistantVector
	PROPERTIES WILL_FAIL true)

# does simple test call fail when solution does exist?
add_test (NoiseAdderSimple_ExistantOutput 
	${PROJECT_BINARY_DIR}/Examples/NoiseAdder/NoiseAdder 
		--input ${PROJECT_SOURCE_DIR}/tests/NoiseAdder/simple/pre/vector.m
		--output solution_vector.mat
		--noise-level 0.1)
set_tests_properties (NoiseAdderSimple_ExistantOutput
	PROPERTIES WILL_FAIL true)
set_tests_properties (NoiseAdderSimple_ExistantOutput
	PROPERTIES DEPENDS NoiseAdderSimple_Vector)
	
