# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Solvers 
	AuxiliaryConstraints/AuxiliaryConstraints.cpp
	AuxiliaryConstraints/AuxiliaryConstraintsFactory.cpp
	InverseProblemSolver.cpp
	SolverFactory/SolverFactory.cpp
	)

install (FILES AuxiliaryConstraints/AuxiliaryConstraints.hpp DESTINATION include)
install (FILES AuxiliaryConstraints/AuxiliaryConstraint_impl.hpp DESTINATION include)
install (FILES AuxiliaryConstraints/AuxiliaryConstraints_AND.hpp DESTINATION include)
install (FILES AuxiliaryConstraints/AuxiliaryConstraintsFactory.hpp DESTINATION include)
install (FILES AuxiliaryConstraints/NonnegativeConstraint.hpp DESTINATION include)
install (FILES AuxiliaryConstraints/NonpositiveConstraint.hpp DESTINATION include)
install (FILES AuxiliaryConstraints/UnityConstraint.hpp DESTINATION include)
install (FILES FeasibilityProblem.hpp DESTINATION include)
install (FILES InverseProblemSolver.hpp DESTINATION include)
install (FILES SolverFactory/SolverFactory.hpp DESTINATION include)
