# test 1a: blackwhite mixed_equal
add_test (MatrixToPNG_blackwhite_mixed_equal
	${PROJECT_BINARY_DIR}/Examples/MatrixToPNG/MatrixToPNG
		--image blackwhite_mixed_equal.png
		--num-pixels-x 21
		--num-pixels-y 1
		--colorize blackwhite
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/mixed_equal.mat)
set_tests_properties (MatrixToPNG_blackwhite_mixed_equal
	PROPERTIES WILL_FAIL false)

# test 1b: checking results
add_test (MatrixToPNG_blackwhite_mixed_equal_check
	diff ${CMAKE_CURRENT_SOURCE_DIR}/post/blackwhite_mixed_equal.png blackwhite_mixed_equal.png)
set_tests_properties (MatrixToPNG_blackwhite_mixed_equal_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixToPNG_blackwhite_mixed_equal_check
	PROPERTIES DEPENDS MatrixToPNG_blackwhite_mixed_equal)


# test 2a: redgreen mixed_equal
add_test (MatrixToPNG_redgreen_mixed_equal
	${PROJECT_BINARY_DIR}/Examples/MatrixToPNG/MatrixToPNG
		--image redgreen_mixed_equal.png
		--num-pixels-x 21
		--num-pixels-y 1
		--colorize redgreen
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/mixed_equal.mat)
set_tests_properties (MatrixToPNG_redgreen_mixed_equal
	PROPERTIES WILL_FAIL false)

# test 2b: checking results
add_test (MatrixToPNG_redgreen_mixed_equal_check
	diff ${CMAKE_CURRENT_SOURCE_DIR}/post/redgreen_mixed_equal.png redgreen_mixed_equal.png)
set_tests_properties (MatrixToPNG_redgreen_mixed_equal_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixToPNG_redgreen_mixed_equal_check
	PROPERTIES DEPENDS MatrixToPNG_redgreen_mixed_equal)


# test 3a: redblue mixed_equal
add_test (MatrixToPNG_redblue_mixed_equal
	${PROJECT_BINARY_DIR}/Examples/MatrixToPNG/MatrixToPNG
		--image redblue_mixed_equal.png
		--num-pixels-x 21
		--num-pixels-y 1
		--colorize redblue
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/mixed_equal.mat)
set_tests_properties (MatrixToPNG_redblue_mixed_equal
	PROPERTIES WILL_FAIL false)

# test 3b: checking results
add_test (MatrixToPNG_redblue_mixed_equal_check
	diff ${CMAKE_CURRENT_SOURCE_DIR}/post/redblue_mixed_equal.png redblue_mixed_equal.png)
set_tests_properties (MatrixToPNG_redblue_mixed_equal_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixToPNG_redblue_mixed_equal_check
	PROPERTIES DEPENDS MatrixToPNG_redblue_mixed_equal)


# test 4a: bluegreenred mixed_equal
add_test (MatrixToPNG_bluegreenred_mixed_equal
	${PROJECT_BINARY_DIR}/Examples/MatrixToPNG/MatrixToPNG
		--image bluegreenred_mixed_equal.png
		--num-pixels-x 21
		--num-pixels-y 1
		--colorize bluegreenred
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/mixed_equal.mat)
set_tests_properties (MatrixToPNG_bluegreenred_mixed_equal
	PROPERTIES WILL_FAIL false)

# test 4b: checking results
add_test (MatrixToPNG_bluegreenred_mixed_equal_check
	diff ${CMAKE_CURRENT_SOURCE_DIR}/post/bluegreenred_mixed_equal.png bluegreenred_mixed_equal.png)
set_tests_properties (MatrixToPNG_bluegreenred_mixed_equal_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixToPNG_bluegreenred_mixed_equal_check
	PROPERTIES DEPENDS MatrixToPNG_bluegreenred_mixed_equal)
		