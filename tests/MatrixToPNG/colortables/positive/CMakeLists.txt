# test 1a: blackwhite positive
add_test (MatrixToPNG_blackwhite_positive
	${PROJECT_BINARY_DIR}/Examples/MatrixToPNG/MatrixToPNG
		--image blackwhite_positive.png
		--num-pixels-x 11 
		--num-pixels-y 1
		--colorize blackwhite
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/positive.mat)
set_tests_properties (MatrixToPNG_blackwhite_positive
	PROPERTIES WILL_FAIL false)

# test 1b: checking results
add_test (MatrixToPNG_blackwhite_positive_check
	diff ${CMAKE_CURRENT_SOURCE_DIR}/post/blackwhite_positive.png blackwhite_positive.png)
set_tests_properties (MatrixToPNG_blackwhite_positive_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixToPNG_blackwhite_positive_check
	PROPERTIES DEPENDS MatrixToPNG_blackwhite_positive)


# test 2a: redgreen positive
add_test (MatrixToPNG_redgreen_positive
	${PROJECT_BINARY_DIR}/Examples/MatrixToPNG/MatrixToPNG
		--image redgreen_positive.png
		--num-pixels-x 11 
		--num-pixels-y 1
		--colorize redgreen
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/positive.mat)
set_tests_properties (MatrixToPNG_redgreen_positive
	PROPERTIES WILL_FAIL false)

# test 2b: checking results
add_test (MatrixToPNG_redgreen_positive_check
	diff ${CMAKE_CURRENT_SOURCE_DIR}/post/redgreen_positive.png redgreen_positive.png)
set_tests_properties (MatrixToPNG_redgreen_positive_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixToPNG_redgreen_positive_check
	PROPERTIES DEPENDS MatrixToPNG_redgreen_positive)


# test 3a: redblue positive
add_test (MatrixToPNG_redblue_positive
	${PROJECT_BINARY_DIR}/Examples/MatrixToPNG/MatrixToPNG
		--image redblue_positive.png
		--num-pixels-x 11 
		--num-pixels-y 1
		--colorize redblue
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/positive.mat)
set_tests_properties (MatrixToPNG_redblue_positive
	PROPERTIES WILL_FAIL false)

# test 3b: checking results
add_test (MatrixToPNG_redblue_positive_check
	diff ${CMAKE_CURRENT_SOURCE_DIR}/post/redblue_positive.png redblue_positive.png)
set_tests_properties (MatrixToPNG_redblue_positive_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixToPNG_redblue_positive_check
	PROPERTIES DEPENDS MatrixToPNG_redblue_positive)


# test 4a: bluegreenred positive
add_test (MatrixToPNG_bluegreenred_positive
	${PROJECT_BINARY_DIR}/Examples/MatrixToPNG/MatrixToPNG
		--image bluegreenred_positive.png
		--num-pixels-x 11 
		--num-pixels-y 1
		--colorize bluegreenred
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/positive.mat)
set_tests_properties (MatrixToPNG_bluegreenred_positive
	PROPERTIES WILL_FAIL false)

# test 4b: checking results
add_test (MatrixToPNG_bluegreenred_positive_check
	diff ${CMAKE_CURRENT_SOURCE_DIR}/post/bluegreenred_positive.png bluegreenred_positive.png)
set_tests_properties (MatrixToPNG_bluegreenred_positive_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixToPNG_bluegreenred_positive_check
	PROPERTIES DEPENDS MatrixToPNG_bluegreenred_positive)
		