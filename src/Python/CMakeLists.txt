# python libs
find_package(PythonLibs 3 REQUIRED)
include_directories( ${PYTHON_INCLUDE_DIRS} )

# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (pyBasso SHARED 
	pybasso.cpp
	utility.cpp)

target_link_libraries(pyBasso 
	Minimizations
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${Boost_LIBRARIES}
	${PYTHON_LIBRARY} 
)

# Copy the __init__.py file
configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py COPYONLY)

# Suppress prefix "lib" because Python does not allow this prefix
set_target_properties(pyBasso PROPERTIES PREFIX "")

#install(TARGETS pyBasso DESTINATION "${PYTHON_INSTALL_PATH}")
