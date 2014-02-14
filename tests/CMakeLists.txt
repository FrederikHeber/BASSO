### TESTS ###
#add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
#                  DEPENDS Basso)


# does the application run
add_test (BassoRuns ${PROJECT_BINARY_DIR}/src/Basso 25)

# does the usage message work?
add_test (BassoUsage ${PROJECT_BINARY_DIR}/src/Basso)
set_tests_properties (BassoUsage
	PROPERTIES PASS_REGULAR_EXPRESSION "Usage:.*number")

# define a macro to simplify adding tests, the use it
macro (do_test arg result)
	add_test (BassoComp${arg} ${PROJECT_BINARY_DIR}/src/Basso ${arg})
	set_tests_properties (BassoComp${arg}
		PROPERTIES PASS_REGULAR_EXPRESSION ${result})
endmacro (do_test)

# do a bunch of result based tests
do_test (4 "4 is 2")
do_test (9 "9 is 3")
do_test (5 "5 is 2.236")
do_test (7 "7 is 2.645")
do_test (25 "25 is 5")
do_test (-25 "-25 is 0")
do_test (0.0001 "0.0001 is 0.01")
