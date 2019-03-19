# look for python interpreter
set(Python_ADDITIONAL_VERSIONS 3)
include(FindPythonInterp)

# look for python modules
include (FindPythonModules)
find_python_module(numpy REQUIRED)
find_python_module(minieigen REQUIRED)

# can we instantiate NormedSpace and play around
add_test(
	NAME pyBassoSpace 
	COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pre/space_test.py 
	WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH} )
set_tests_properties (pyBassoSpace
	PROPERTIES WILL_FAIL false)	# return code 0
set_tests_properties (pyBassoSpace PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/src/Python)
    
# can we instantiate SpaceElement and play around
add_test(
	NAME pyBassoVector 
	COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pre/vector_test.py 
	WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH} )
set_tests_properties (pyBassoVector
	PROPERTIES WILL_FAIL false)	# return code 0
set_tests_properties (pyBassoVector PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/src/Python)

# can we instantiate LinearMapping and play around
add_test(
	NAME pyBassoLinearMapping 
	COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pre/linearmapping_test.py 
	WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH} )
set_tests_properties (pyBassoLinearMapping
	PROPERTIES WILL_FAIL false)	# return code 0
set_tests_properties (pyBassoLinearMapping PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/src/Python)

# can we instantiate CommandLineOptions and modify the options
add_test(
	NAME pyBassoCommandLineOptions 
	COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pre/commandlineoptions_test.py 
	WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH} )
set_tests_properties (pyBassoCommandLineOptions
	PROPERTIES WILL_FAIL false)	# return code 0
set_tests_properties (pyBassoCommandLineOptions PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/src/Python)

# can we instantiate InverseProblem and set up a random matrix toy problem
add_test(
	NAME pyBassoInverseProblem 
	COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pre/inverseproblem_test.py 
	WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH} )
set_tests_properties (pyBassoInverseProblem
	PROPERTIES WILL_FAIL false)	# return code 0
set_tests_properties (pyBassoInverseProblem PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/src/Python)

# can we instantiate InverseProblem and set up a (non-)linear operator toy problem
add_test(
	NAME pyBassoNonLinearInverseProblem 
	COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pre/nonlinear_inverseproblem_test.py 
	WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH} )
set_tests_properties (pyBassoNonLinearInverseProblem
	PROPERTIES WILL_FAIL false)	# return code 0
set_tests_properties (pyBassoNonLinearInverseProblem PROPERTIES
    ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/src/Python)

