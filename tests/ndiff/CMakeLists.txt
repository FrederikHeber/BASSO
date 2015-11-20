# does the application run
add_test (ndiffNotXml ${NDIFF})
set_tests_properties (ndiffNotXml
	PROPERTIES PASS_REGULAR_EXPRESSION ".*-precision number-of-bits")

add_test (ndiff_simple1
	${DIFFCMD}
	-q -relerr 1e-5
	${CMAKE_CURRENT_SOURCE_DIR}/pre/infile1
	${CMAKE_CURRENT_SOURCE_DIR}/pre/infile2)
set_tests_properties (ndiff_simple1
	PROPERTIES WILL_FAIL TRUE)

add_test (ndiff_simple2 
	${DIFFCMD}
	-q -relerr 1e-4
	${CMAKE_CURRENT_SOURCE_DIR}/pre/infile1
	${CMAKE_CURRENT_SOURCE_DIR}/pre/infile2)
set_tests_properties (ndiff_simple2
	PROPERTIES WILL_FAIL FALSE)
