# test 1a: blackwhite more_negative
add_test (MatrixToPNG_blackwhite_more_negative
	${PROJECT_BINARY_DIR}/Examples/MatrixToPNG/MatrixToPNG
		--image blackwhite_more_negative.png
		--num-pixels-x 17
		--num-pixels-y 1
		--colorize blackwhite
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/more_negative.mat)
set_tests_properties (MatrixToPNG_blackwhite_more_negative
	PROPERTIES WILL_FAIL false)

# test 1b: checking results
add_test (MatrixToPNG_blackwhite_more_negative_check
	diff ${CMAKE_CURRENT_SOURCE_DIR}/post/blackwhite_more_negative.png blackwhite_more_negative.png)
set_tests_properties (MatrixToPNG_blackwhite_more_negative_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixToPNG_blackwhite_more_negative_check
	PROPERTIES DEPENDS MatrixToPNG_blackwhite_more_negative)


# test 2a: redgreen more_negative
add_test (MatrixToPNG_redgreen_more_negative
	${PROJECT_BINARY_DIR}/Examples/MatrixToPNG/MatrixToPNG
		--image redgreen_more_negative.png
		--num-pixels-x 17
		--num-pixels-y 1
		--colorize redgreen
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/more_negative.mat)
set_tests_properties (MatrixToPNG_redgreen_more_negative
	PROPERTIES WILL_FAIL false)

# test 2b: checking results
add_test (MatrixToPNG_redgreen_more_negative_check
	diff ${CMAKE_CURRENT_SOURCE_DIR}/post/redgreen_more_negative.png redgreen_more_negative.png)
set_tests_properties (MatrixToPNG_redgreen_more_negative_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixToPNG_redgreen_more_negative_check
	PROPERTIES DEPENDS MatrixToPNG_redgreen_more_negative)


# test 3a: redblue more_negative
add_test (MatrixToPNG_redblue_more_negative
	${PROJECT_BINARY_DIR}/Examples/MatrixToPNG/MatrixToPNG
		--image redblue_more_negative.png
		--num-pixels-x 17
		--num-pixels-y 1
		--colorize redblue
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/more_negative.mat)
set_tests_properties (MatrixToPNG_redblue_more_negative
	PROPERTIES WILL_FAIL false)

# test 3b: checking results
add_test (MatrixToPNG_redblue_more_negative_check
	diff ${CMAKE_CURRENT_SOURCE_DIR}/post/redblue_more_negative.png redblue_more_negative.png)
set_tests_properties (MatrixToPNG_redblue_more_negative_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixToPNG_redblue_more_negative_check
	PROPERTIES DEPENDS MatrixToPNG_redblue_more_negative)


# test 4a: bluegreenred more_negative
add_test (MatrixToPNG_bluegreenred_more_negative
	${PROJECT_BINARY_DIR}/Examples/MatrixToPNG/MatrixToPNG
		--image bluegreenred_more_negative.png
		--num-pixels-x 17
		--num-pixels-y 1
		--colorize bluegreenred
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/more_negative.mat)
set_tests_properties (MatrixToPNG_bluegreenred_more_negative
	PROPERTIES WILL_FAIL false)

# test 4b: checking results
add_test (MatrixToPNG_bluegreenred_more_negative_check
	diff ${CMAKE_CURRENT_SOURCE_DIR}/post/bluegreenred_more_negative.png bluegreenred_more_negative.png)
set_tests_properties (MatrixToPNG_bluegreenred_more_negative_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixToPNG_bluegreenred_more_negative_check
	PROPERTIES DEPENDS MatrixToPNG_bluegreenred_more_negative)
		