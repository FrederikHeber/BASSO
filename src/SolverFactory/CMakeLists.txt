# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (SolverFactory 
	SolverFactory.cpp
	)

install (FILES SolverFactory.hpp DESTINATION include)
