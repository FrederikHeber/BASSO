# fails because no measurement dimensions
add_test (RadonMatrixWriterSimple_failmeasurementdimensions 
	${PROJECT_BINARY_DIR}/Examples/RadonMatrixWriter/RadonMatrixWriter 
		--radon-matrix radon.mat)
set_tests_properties (RadonMatrixWriterSimple_failmeasurementdimensions
	PROPERTIES WILL_FAIL true)

# fails because of missing rhs file 	
add_test (RadonMatrixWriterSimple_failmissingfile 
	${PROJECT_BINARY_DIR}/Examples/RadonMatrixWriter/RadonMatrixWriter 
		--num-pixels-x 2 
		--num-pixels-y 2
		--num-angles 2 
		--num-offsets 3)
set_tests_properties (RadonMatrixWriterSimple_failmissingfile
	PROPERTIES WILL_FAIL true)

# succeeds 	
add_test (RadonMatrixWriterSimple_succeeds 
	${PROJECT_BINARY_DIR}/Examples/RadonMatrixWriter/RadonMatrixWriter 
		--num-pixels-x 2 
		--num-pixels-y 2
		--num-angles 2 
		--num-offsets 3
		--radon-matrix radon.mat)
set_tests_properties (RadonMatrixWriterSimple_succeeds
	PROPERTIES WILL_FAIL false)

# test 1b: checking results
add_test (RadonMatrixWriterSimple_succeeds_check
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/radon.mat radon.mat)
set_tests_properties (RadonMatrixWriterSimple_succeeds_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (RadonMatrixWriterSimple_succeeds_check
	PROPERTIES DEPENDS RadonMatrixWriterSimple_succeeds)
