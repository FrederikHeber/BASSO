##########################################################################
########################      Basso      #################################
##########################################################################


cmake_minimum_required (VERSION 2.6)
project (Basso C CXX)

# set debug and release flags
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g3 -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-g0 -O3")

# The version number.
set (Basso_VERSION_MAJOR 1)
set (Basso_VERSION_MINOR 0)

# should we count and measure the time of all operations?
option (USE_TIMINGS 
	"Use counter and timings on all vector space operations" OFF)

# use operator norm instead of matrix norm
option (FULLMATRIXNORM
	"Use operator norm instead of matrix norm, requires SVD calculation" OFF)

######################## CHECK HEADERS  ##################################

# does the system provide the log and exp functions?
set(CMAKE_REQUIRED_INCLUDES math.h)
set(CMAKE_REQUIRED_LIBRARIES m)
include (${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
check_function_exists (log HAVE_LOG)
check_function_exists (exp HAVE_EXP)

####################### CHECK LIBRARIES  #################################

# find eigen3 package for linear algebra routines
set ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/  )
FIND_PACKAGE(Eigen3 REQUIRED)
include_directories( "${EIGEN3_INCLUDE_DIR}" )

# find cppunit for unit tests
set ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/  )
FIND_PACKAGE(CppUnit REQUIRED)
if (CPPUNIT_FOUND)
	include_directories( "${CPPUNIT_INCLUDE_DIR}" )
	#set (EXTRA_LIBS ${EXTRA_LIBS} "${CPPUNIT_LIBRARY}")
endif (CPPUNIT_FOUND)

# looking for blas
set (BLA_VENDOR "Generic")
FIND_PACKAGE(BLAS REQUIRED)
if (BLAS_FOUND)
       # find gsl for function minimization
       set ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/  )
       FIND_PACKAGE(GSL REQUIRED)
       if (GSL_FOUND)
               include_directories( "${GSL_INCLUDE_DIR}" )
               set (GSL_LIBRARY ${GSL_LIBRARIES} ${BLAS_blas_LIBRARY})
       endif (GSL_FOUND)
endif (BLAS_FOUND)

# find NLopt
set ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/  )
find_package (NLopt)
if (NLOPT_FOUND)
               include_directories( "${NLOPT_INCLUDE_DIR}" )
endif (NLOPT_FOUND)

# find (p)thread
find_package (Threads)

# find librt (required by boost_thread)
set ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/  )
find_package (Librt)

# find boost
set ( Boost_ADDITIONAL_VERSIONS "1.54" "1.54.0" "1.55" "1.55.0" )
set ( Boost_USE_STATIC_LIBS ON )
find_package ( Boost 1.54 COMPONENTS
	chrono
	log_setup 
	log 
	filesystem
	thread
	program_options 
	system 
	REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )

# find pocosqlite
set ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/  )
find_package (Poco)
include_directories(${Poco_INCLUDE_DIRS})
link_directories(${Poco_LIBRARY_DIR})
set (POCO_SQL_LIBRARIES
	${Poco_Data_LIBRARY_RELEASE}
	${Poco_SQLite_LIBRARY_RELEASE}
	${Poco_Foundation_LIBRARY_RELEASE}
	)

######################## BassoConfig.hpp #################################

# configure a header file to pass some CMake settings
# to the source code
configure_file (
	"${PROJECT_SOURCE_DIR}/BassoConfig.h.in"
	"${PROJECT_BINARY_DIR}/BassoConfig.h"
)

# add the binary tree to the search path for include files
# so that we will find BassoConfig.h
include_directories ("${PROJECT_BINARY_DIR}")

# install targets for config
install (FILES "${PROJECT_BINARY_DIR}/BassoConfig.h" DESTINATION include)

################## BRANCH INTO SUBDIRECTORIES ############################

#### LIBRARIES ####

# include source folder
include_directories ("${PROJECT_SOURCE_DIR}/src")

# add the CommandLineOptions library
include_directories ("${PROJECT_SOURCE_DIR}/CommandLineOptions")
add_subdirectory (src/CommandLineOptions)
set (EXTRA_LIBS ${EXTRA_LIBS} CommandLineOptions)

# add the MatrixIO library
include_directories ("${PROJECT_SOURCE_DIR}/MatrixIO")
add_subdirectory (src/MatrixIO)
# is header only

# add the SolutionFactory folder
add_subdirectory (src/SolutionFactory)
set (EXTRA_LIBS ${EXTRA_LIBS} SolutionFactory)

# add the Minimizations library
include_directories ("${PROJECT_SOURCE_DIR}/Minimizations")
add_subdirectory (src/Minimizations)
set (EXTRA_LIBS ${EXTRA_LIBS} Minimizations)

# add the Math library
add_subdirectory (src/Math)
set (EXTRA_LIBS ${EXTRA_LIBS} Math)

# add the Log folder
add_subdirectory (src/Log)
set (EXTRA_LIBS ${EXTRA_LIBS} Log)

# add the Database folder
add_subdirectory (src/Database)
set (EXTRA_LIBS ${EXTRA_LIBS} Database)

#### SOURCE ####

# add src folder
add_subdirectory (src)

#### TESTS AND UNITTESTS ####

# add test folder
enable_testing ()	# must reside here
add_subdirectory (tests)
add_subdirectory (src/Database/unittests)
add_subdirectory (src/Math/unittests)
add_subdirectory (src/Minimizations/Elements/unittests)
add_subdirectory (src/Minimizations/Functions/unittests)
add_subdirectory (src/Minimizations/Mappings/unittests)
add_subdirectory (src/Minimizations/Norms/unittests)

##################### PACKING INSTRUCTIONS ###############################

# build a CPack driven installer package
include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE
	"${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set (CPACK_PACKAGE_VERSION_MAJOR "${Basso_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${Basso_VERSION_MINOR}")
include (CPack)

##########################################################################
