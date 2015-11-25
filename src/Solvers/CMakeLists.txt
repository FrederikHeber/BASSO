# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Solvers 
	InverseProblemSolver.cpp
	SolverFactory/SolverFactory.cpp
	)

install (FILES FeasibilityProblem.hpp DESTINATION include)
install (FILES InverseProblemSolver.hpp DESTINATION include)
install (FILES SolverFactory/SolverFactory.hpp DESTINATION include)
