# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (SolutionFactory 
	SolutionFactory.cpp
	)

install (FILES SolutionFactory.hpp DESTINATION include)
