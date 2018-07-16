# python libs
find_package(PythonLibs 2.7 REQUIRED)
include_directories( ${PYTHON_INCLUDE_DIRS} )

# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (pyBasso SHARED 
	pybasso.cpp)

target_link_libraries(pyBasso ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

# Copy the __init__.py file
configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py COPYONLY)

# Suppress prefix "lib" because Python does not allow this prefix
set_target_properties(pyBasso PROPERTIES PREFIX "")

#install(TARGETS pyBasso DESTINATION "${PYTHON_INSTALL_PATH}")
