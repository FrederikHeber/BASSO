# does simple test call work?
add_test (NonnegativeMatrixFactorsSimple 
	${PROJECT_BINARY_DIR}/Examples/NonnegativeMatrixFactors/NonnegativeMatrixFactors 
		--source-first-factor ${CMAKE_CURRENT_SOURCE_DIR}/pre/spectralmatrix.m
		--source-second-factor ${CMAKE_CURRENT_SOURCE_DIR}/pre/pixelmatrix.m
		--destination-first-factor spectralmatrix.m
		--destination-second-factor pixelmatrix.m
		)
set_tests_properties (NonnegativeMatrixFactorsSimple
	PROPERTIES WILL_FAIL false)
