# look for python interpreter
set(Python_ADDITIONAL_VERSIONS 3)
include(FindPythonInterp)

# python libs
find_package(PythonLibs 3 REQUIRED)
include_directories( ${PYTHON_INCLUDE_DIRS} )

# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (pyBasso SHARED 
	pybasso.cpp
	utility.cpp
	Wrappers/CommandLineOptions.cpp
	Wrappers/InverseProblem.cpp
	Wrappers/Mapping.cpp
	Wrappers/Norm.cpp
	Wrappers/NormedSpace.cpp
	Wrappers/SpaceElement.cpp
	)

target_link_libraries(pyBasso 
	Solvers
	Database
	Options
	Minimizations
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${BOOST_CORE_LIBRARIES}
	${Boost_PYTHON3_LIBRARY}
	${PYTHON_LIBRARY} 
)

# Copy the __init__.py file
configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py COPYONLY)

# Suppress prefix "lib" because Python does not allow this prefix
set_target_properties(pyBasso PROPERTIES PREFIX "")

#install(TARGETS pyBasso DESTINATION lib)
#install(CODE "execute_process(COMMAND python setup.py install of --prefix=$${CMAKE_INSTALL_PREFIX} WORKING_DIRECTORY ./)")

# the following is taken from https://cmake.org/pipermail/cmake/2013-March/054138.html
execute_process (
   COMMAND ${PYTHON_EXECUTABLE} -c
   	"import site, sys; sys.stdout.write(site.PREFIXES[-1])"
   OUTPUT_VARIABLE PYTHON_PREFIX
)
file ( TO_CMAKE_PATH "${PYTHON_PREFIX}" PYTHON_PREFIX )
execute_process (
   COMMAND ${PYTHON_EXECUTABLE} -c
   	"import site, sys; sys.stdout.write(site.getsitepackages()[-1])"
   OUTPUT_VARIABLE PYTHON_SITE_DIR
)
file ( TO_CMAKE_PATH "${PYTHON_SITE_DIR}" PYTHON_SITE_DIR )
string ( REGEX REPLACE "^${PYTHON_PREFIX}/" ""
   PYTHON_SITE_DIR "${PYTHON_SITE_DIR}"
)

install ( TARGETS pyBasso
   LIBRARY
     DESTINATION ${PYTHON_SITE_DIR}
     COMPONENT library
)

