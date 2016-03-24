# test 1: does simple test call work?
add_test (NonnegativeMatrixFactorsSimple 
	${PROJECT_BINARY_DIR}/Examples/NonnegativeMatrixFactors/NonnegativeMatrixFactors 
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/matrix.mat
		--truncation-dimension 3
		--destination-first-factor spectralmatrix.mat
		--destination-second-factor pixelmatrix.mat
		)
set_tests_properties (NonnegativeMatrixFactorsSimple
	PROPERTIES WILL_FAIL false)

# test 1b: checking factor was copied
add_test (NonnegativeMatrixFactorsSimple_check_first
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/spectralmatrix.mat spectralmatrix.mat)
set_tests_properties (NonnegativeMatrixFactorsSimple_check_first
	PROPERTIES WILL_FAIL false)
set_tests_properties (NonnegativeMatrixFactorsSimple_check_first
	PROPERTIES DEPENDS NonnegativeMatrixFactorsSimple)

# test 1c: checking factor was copied
add_test (NonnegativeMatrixFactorsSimple_check_second
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/pixelmatrix.mat pixelmatrix.mat)
set_tests_properties (NonnegativeMatrixFactorsSimple_check_second
	PROPERTIES WILL_FAIL false)
set_tests_properties (NonnegativeMatrixFactorsSimple_check_second
	PROPERTIES DEPENDS NonnegativeMatrixFactorsSimple)
	
