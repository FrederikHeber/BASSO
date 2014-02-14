cmake_minimum_required (VERSION 2.6)
project (Basso)

# The version number.
set (Basso_VERSION_MAJOR 1)
set (Basso_VERSION_MINOR 0)

# does the system provide the log and exp functions?
include (${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
check_function_exists (log HAVE_LOG)
check_function_exists (exp HAVE_EXP)

# configure a header file to pass some CMake settings
# to the source code
configure_file (
	"${PROJECT_SOURCE_DIR}/BassoConfig.h.in"
	"${PROJECT_BINARY_DIR}/BassoConfig.h"
)

# add the binary tree to the search path for include files
# so that we will find BassoConfig.h
include_directories ("${PROJECT_BINARY_DIR}")

# should we use our own math functions?
option (USE_MYMATH 
	"Use tutorial provided math implementation" ON) 

# add the Minimizations library?
#
if (USE_MYMATH)
  include_directories ("${PROJECT_SOURCE_DIR}/Minimizations")
  add_subdirectory (Minimizations)
  set (EXTRA_LIBS ${EXTRA_LIBS} Minimizations)
endif (USE_MYMATH)

# add the executable
add_executable (Basso Basso.cpp)
target_link_libraries (Basso ${EXTRA_LIBS})

# add the install targets
install (TARGETS Basso DESTINATION bin)
install (FILES "${PROJECT_BINARY_DIR}/BassoConfig.h" DESTINATION include)

### TESTS ###
enable_testing ()

# does the application run
add_test (BassoRuns Basso 25)

# does the usage message work?
add_test (BassoUsage Basso)
set_tests_properties (BassoUsage
	PROPERTIES PASS_REGULAR_EXPRESSION "Usage:.*number")

# define a macro to simplify adding tests, the use it
macro (do_test arg result)
	add_test (BassoComp${arg} Basso ${arg})
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

# build a CPack driven installer package
include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE
	"${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set (CPACK_PACKAGE_VERSION_MAJOR "${Basso_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${Basso_VERSION_MINOR}")
include (CPack)
