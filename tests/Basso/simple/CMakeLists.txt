# does simple test call work?
add_test (BassoSimple 
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01 
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/matrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/vector.m)
set_tests_properties (BassoSimple
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of X.*2")
set_tests_properties (BassoSimple
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of Y.*2")
set_tests_properties (BassoSimple
	PROPERTIES PASS_REGULAR_EXPRESSION ".*noise.*0.01")
set_tests_properties (BassoSimple
	PROPERTIES PASS_REGULAR_EXPRESSION ".*Ax = y.*1 0")

# does simple test call fail when matrix does not exist?
add_test (BassoSimple_NonExistantMatrix 
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01 
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/nomatrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/vector.m)
set_tests_properties (BassoSimple_NonExistantMatrix
	PROPERTIES WILL_FAIL true)

# tuple params working
add_test (BassoSimple_TupleParameters 
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01 
		--iteration-file simple.db
		--tuple-parameters "seed" "426"
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/matrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/vector.m)
set_tests_properties (BassoSimple_TupleParameters
	PROPERTIES PASS_REGULAR_EXPRESSION "tuple parameters.*seed,426,")
	
