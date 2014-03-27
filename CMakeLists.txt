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

# find eigen3 package for linear algebra routines
set ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/  )
FIND_PACKAGE(Eigen3 REQUIRED)
include_directories( "${EIGEN3_INCLUDE_DIR}" )

# add the Minimizations library?
#
if (USE_MYMATH)
  include_directories ("${PROJECT_SOURCE_DIR}/Minimizations")
  add_subdirectory (src/Minimizations)
  set (EXTRA_LIBS ${EXTRA_LIBS} Minimizations)
endif (USE_MYMATH)

# install targets for config
install (FILES "${PROJECT_BINARY_DIR}/BassoConfig.h" DESTINATION include)

# add src folder
add_subdirectory (src)

# add test folder
enable_testing ()	# must reside here
add_subdirectory (tests)

# build a CPack driven installer package
include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE
	"${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set (CPACK_PACKAGE_VERSION_MAJOR "${Basso_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${Basso_VERSION_MINOR}")
include (CPack)
