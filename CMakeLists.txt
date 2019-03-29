##########################################################################
########################      Basso      #################################
##########################################################################


cmake_minimum_required (VERSION 2.6)
project (Basso C CXX)

# set debug and release flags
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g3 -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-g0 -O3 -DNDEBUG")

# The version number.
set (Basso_VERSION_MAJOR 1)
set (Basso_VERSION_MINOR 11)
set (Basso_VERSION_PATCH 1)

set (ISSUE_REPORT_MSG "Please consider opening an issue at https://github.com/FrederikHeber/BASSO with anything that is helpful.")

# should we count and measure the time of all operations?
option (USE_TIMINGS 
	"Use counter and timings on all vector space operations" OFF)

# use operator norm instead of matrix norm
option (FULLMATRIXNORM
	"Use operator norm instead of matrix norm, requires SVD calculation" OFF)

# use MPI or OpenMP
option (USE_MPI
	"Use MPI parallelization for MatrixFactorizer" OFF)
option (USE_OPENMP
	"Use OPENMP parallelization for MatrixFactorizer" OFF)
if (USE_MPI)
	if (USE_OPENMP)
		error("Cannot activate both MPI and OpenMP parallelization.")
	endif (USE_OPENMP)
	set (USE_MPI $USE_MPI)
endif (USE_MPI)
if (USE_OPENMP)
	set (USE_OPENMP $USE_OPENMP)
endif (USE_OPENMP)

###################### GIT VERSION HEADER  ################################

# create a version target which must be added to every executable
FIND_PACKAGE(Git)
FILE(WRITE ${CMAKE_BINARY_DIR}/version.h.in
"\#define VERSION \"@VERSION@\""
)
FILE(WRITE ${CMAKE_BINARY_DIR}/version.cmake
"EXECUTE_PROCESS(
     COMMAND ${GIT_EXECUTABLE} describe --dirty --always
     OUTPUT_VARIABLE VERSION
     OUTPUT_STRIP_TRAILING_WHITESPACE
 )
 CONFIGURE_FILE(\${SRC} \${DST} @ONLY)
")
#INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})
ADD_CUSTOM_TARGET(
    version
    ${CMAKE_COMMAND} -D SRC=${CMAKE_BINARY_DIR}/version.h.in
                     -D DST=${CMAKE_BINARY_DIR}/version.h
                     -P ${CMAKE_BINARY_DIR}/version.cmake
)

######################## CHECK HEADERS  ##################################

# does the system provide the log and exp functions?
set(CMAKE_REQUIRED_INCLUDES math.h)
set(CMAKE_REQUIRED_LIBRARIES m)
include (${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
check_function_exists (log HAVE_LOG)
check_function_exists (exp HAVE_EXP)

# find MPI
if (USE_MPI)
	set ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/  )
	FIND_PACKAGE(MPI REQUIRED)
	if (MPI_FOUND)
		include_directories( ${MPI_CXX_INCLUDE_PATH})
		FIND_PROGRAM(MPIRUN NAME mpirun PATHS usr/bin/ /bin ENV PATH NO_DEFAULT_PATH)
		set (MPIRUN_FOUND ${MPIRUN_FOUND})
		set (MPIRUN ${MPIRUN} CACHE STRING "command to launch multiple mpi jobs")
	endif (MPI_FOUND)
endif (USE_MPI)

# find OpenMP
if (USE_OPENMP)
	set ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/  )
	FIND_PACKAGE(OpenMP REQUIRED)
endif (USE_OPENMP)

####################### CHECK LIBRARIES  #################################

# find eigen3 package for linear algebra routines
set ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/  )
FIND_PACKAGE(Eigen3 3.2.0)
include_directories( ${EIGEN3_INCLUDE_DIR} )

# find cppunit for unit tests
set ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/  )
FIND_PACKAGE(CppUnit REQUIRED)
if (CPPUNIT_FOUND)
	include_directories( ${CPPUNIT_INCLUDE_DIR} )
	#set (EXTRA_LIBS ${EXTRA_LIBS} ${CPPUNIT_LIBRARY})
endif (CPPUNIT_FOUND)

# looking for blas
set (BLA_VENDOR "Generic")
FIND_PACKAGE(BLAS REQUIRED)
if (BLAS_FOUND)
       # find gsl for function minimization
       set ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/  )
       FIND_PACKAGE(GSL REQUIRED)
       if (GSL_FOUND)
               include_directories( ${GSL_INCLUDE_DIR} )
               set (GSL_LIBRARY ${GSL_LIBRARIES} ${BLAS_blas_LIBRARY})
       endif (GSL_FOUND)
endif (BLAS_FOUND)

# find NLopt
set ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/  )
find_package (NLopt)
if (NLOPT_FOUND)
               include_directories( ${NLOPT_INCLUDE_DIR} )
endif (NLOPT_FOUND)

# find (p)thread
find_package (Threads REQUIRED)

# find librt (required by boost_thread)
set ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/  )
find_package (Librt REQUIRED)

# find boost
ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)
set ( Boost_ADDITIONAL_VERSIONS "1.61" "1.61.0" )
set ( Boost_USE_STATIC_LIBS OFF )
if (MPI_FOUND)
	find_package ( Boost 1.61 COMPONENTS
		chrono
		log_setup 
		log 
		filesystem
		mpi
		program_options
		python3
		serialization 
		system 
		thread
		REQUIRED )
else (MPI_FOUND)
	find_package ( Boost 1.61 COMPONENTS
		chrono
		log_setup 
		log 
		filesystem
		program_options
		python3
		serialization 
		system 
		thread
		REQUIRED )
endif (MPI_FOUND)
set (BOOST_CORE_LIBRARIES
	${Boost_CHRONO_LIBRARY}
	${Boost_LOG_LIBRARY}
	${Boost_FILESYSTEM_LIBRARY}
	${Boost_PROGRAM_OPTIONS_LIBRARY}
	${Boost_SYSTEM_LIBRARY})
include_directories( ${Boost_INCLUDE_DIR} )

# find pocosqlite
set ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/  )
find_package (Poco REQUIRED)
include_directories(${Poco_INCLUDE_DIRS})
link_directories(${Poco_LIBRARY_DIR})
if (Poco_SQLite_LIBRARY_RELEASE)
	set (POCO_SQL_LIBRARIES
		${Poco_Data_LIBRARY_RELEASE}
		${Poco_SQLite_LIBRARY_RELEASE}
		${Poco_Foundation_LIBRARY_RELEASE}
		)
else(Poco_SQLite_LIBRARY_RELEASE)
	set (POCO_SQL_LIBRARIES
		${Poco_Data_LIBRARY_RELEASE}
		${Poco_DataSQLite_LIBRARY_RELEASE}
		${Poco_Foundation_LIBRARY_RELEASE}
		)
endif(Poco_SQLite_LIBRARY_RELEASE)

# find libpng
find_package(PNG)
if (PNG_FOUND)
	find_package(png++)
	if (png++_FOUND)
		include_directories(${PNG_INCLUDE_DIRS} ${png++_INCLUDE_DIRS})
		set (PNG_LIBRARIES ${PNG_LIBRARIES})
		set( PNG_FOUND ${PNG_FOUND})
		set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${PNG_DEFINITIONS})
	endif (png++_FOUND)
endif (PNG_FOUND)

# find libarpack
find_package(ARPACK)
if (ARPACK_FOUND)
	set (ARPACK_LIBRARIES ${ARPACK_LIBRARIES})
	set( ARPACK_FOUND ${ARPACK_FOUND})
endif (ARPACK_FOUND)

####################### CHECK given DIFF #################################

# should ndiff be used instead of diff
option (USE_NDIFF
	"Use ndiff instead of diff" OFF)

if (USE_NDIFF)
	set (NDIFF_PATH ${NDIFF_PATH} CACHE STRING "path to numerical diff")
	FIND_PROGRAM(NDIFF NAME ndiff PATHS ${NDIFF_PATH} /usr/bin/ /bin ENV PATH NO_DEFAULT_PATH)
	if (NDIFF EQUAL NDIFF-NOTFOUND)
		set (DIFFCMD diff CACHE STRING "command to compare result files")
	else (NDIFF EQUAL NDIFF-NOTFOUND)
		set (NDIFF_FOUND ${NDIFF_FOUND})
		set (DIFFCMD ${NDIFF} CACHE STRING "command to compare result files")
	endif (NDIFF EQUAL NDIFF-NOTFOUND)
	if (NDIFF_OPTIONS)
		set (NDIFF_OPTIONS ${NDIFF_OPTIONS} CACHE STRING "Options for (numerical) diff")
	endif (NDIFF_OPTIONS)
else (USE_NDIFF)
	set (DIFFCMD diff CACHE STRING "command to compare result files")
endif (USE_NDIFF)

###################### CHECK given SQLite ################################

# sqlite is only used in tests which are skipped of not found
FIND_PROGRAM(SQLITE NAME sqlite sqlite3 PATHS /usr/bin/ /bin ENV PATH NO_DEFAULT_PATH)
if (NOT SQLITE EQUAL SQLITE-NOTFOUND)
	set (SQLITECMD ${SQLITE} CACHE STRING "command to extract values from sqlite database")
	set (SQLITE_FOUND TRUE CACHE BOOL "whether sqlite was found or not")
else (NOT SQLITE EQUAL SQLITE-NOTFOUND)
	set (SQLITE_FOUND FALSE CACHE BOOL "whether sqlite was found or not")
endif (NOT SQLITE EQUAL SQLITE-NOTFOUND)

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

######################## RPATH settings ##################################
#
# see https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
#
##########################################################################

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

################## BRANCH INTO SUBDIRECTORIES ############################

#### LIBRARIES ####

# include source folder
include_directories ("${PROJECT_SOURCE_DIR}/src")

# include Examples folder
include_directories ("${PROJECT_SOURCE_DIR}/Examples")

# add the Options library
add_subdirectory (src/Options)
set (EXTRA_LIBS ${EXTRA_LIBS} Options)

# add the MatrixIO library
add_subdirectory (src/MatrixIO)
# is header only

# add the Solvers folder
add_subdirectory (src/Solvers)
set (EXTRA_LIBS ${EXTRA_LIBS} Solvers)

# add the Minimizations library
add_subdirectory (src/Minimizations)
set (EXTRA_LIBS ${EXTRA_LIBS} Minimizations)

# add the Math library
add_subdirectory (src/Math)
set (EXTRA_LIBS ${EXTRA_LIBS} Math)

# add the Database folder
add_subdirectory (src/Database)
set (EXTRA_LIBS ${EXTRA_LIBS} Database)

# add the Log folder
add_subdirectory (src/Log)
set (EXTRA_LIBS ${EXTRA_LIBS} Log)

#### PYTHON BINDINGS ####

add_subdirectory(src/Python)

#### EXECUTABLES ####

# add Examples folders of executables
add_subdirectory (Examples)

#### TESTS AND UNITTESTS ####

# add test folder
enable_testing ()	# must reside here
add_subdirectory (tests)
add_subdirectory (src/Database/unittests)
add_subdirectory (src/Math/unittests)
add_subdirectory (src/Minimizations/Elements/unittests)
add_subdirectory (src/Minimizations/Functions/unittests)
add_subdirectory (src/Minimizations/Mappings/unittests)
add_subdirectory (src/Minimizations/Minimizers/StoppingCriteria/unittests)
add_subdirectory (src/Minimizations/Norms/unittests)
add_subdirectory (src/Solvers/AuxiliaryConstraints/unittests)
add_subdirectory (Examples/ComputerTomography/DiscretizedRadon/unittests)

### DOXYGEN DOCUMENTATION

add_subdirectory (src/Documentation)

### Userguide

add_subdirectory (doc)

######################## RPATH settings ##################################
#
# see https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
#
##########################################################################

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

##################### PACKING INSTRUCTIONS ###############################

# build a CPack driven installer package
include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE
	"${CMAKE_CURRENT_SOURCE_DIR}/License.txt")

SET(CPACK_SET_DESTDIR "on")
SET(CPACK_PACKAGING_INSTALL_PREFIX "/tmp")
set(CPACK_GENERATOR "DEB")

SET(CPACK_PACKAGE_DESCRIPTION "Library for non-linear optimization in Banach spaces")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Solve (non-)linear system of equations in arbitrary lp-spaces")
SET(CPACK_PACKAGE_VENDOR "Frederik Heber")
SET(CPACK_PACKAGE_CONTACT "frederik.heber@gmail.com")
set (CPACK_PACKAGE_VERSION_MAJOR "${Basso_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${Basso_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${Basso_VERSION_PATCH}")
SET(CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${Basso_VERSION_MAJOR}.${Basso_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${Basso_VERSION_MAJOR}.${Basso_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Frederik Heber")
SET(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
SET(CPACK_DEBIAN_PACKAGE_SECTION "science")
SET(CPACK_DEBIAN_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSOR})

#set (CPACK_DEBIAN_PACKAGE_DEPENDS 
#	"libboost-chrono (>=1.61), 
#	libboost-log (>=1.61),
#	libboost-filesystem (>=1.61),
#	libboost-program_options (>=1.61),
#	libboost-system (>=1.61)"
#	)

include (CPack)

##########################################################################
