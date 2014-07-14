### TESTS ###

########################## SignumUnitTest #######################
add_executable (SignumUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	SignumUnitTest.cpp) 
target_link_libraries (SignumUnitTest 
	${EXTRA_LIBS}
	"${CPPUNIT_LIBRARY}")

add_test (SignumUnitTest 
	${PROJECT_BINARY_DIR}/src/Math/unittests/SignumUnitTest)


# mark as to fail if required
#set_tests_properties (SignumUnitTest
#	PROPERTIES WILL_FAIL true)
