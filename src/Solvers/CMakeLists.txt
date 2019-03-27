# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Solvers SHARED
	AuxiliaryConstraints/AuxiliaryConstraints.cpp
	AuxiliaryConstraints/AuxiliaryConstraintsFactory.cpp
	InverseProblemSolver.cpp
	RangeProjectionSolver.cpp
	SolverFactory/SolverFactory.cpp
	SplitFeasibilitySolver.cpp
	)

target_link_libraries(Solvers
	Minimizations
	Log
)

install (TARGETS Solvers DESTINATION lib)

install (FILES AuxiliaryConstraints/AuxiliaryConstraints.hpp DESTINATION include)
install (FILES AuxiliaryConstraints/AuxiliaryConstraint_impl.hpp DESTINATION include)
install (FILES AuxiliaryConstraints/AuxiliaryConstraints_AND.hpp DESTINATION include)
install (FILES AuxiliaryConstraints/AuxiliaryConstraintsFactory.hpp DESTINATION include)
install (FILES AuxiliaryConstraints/NonnegativeConstraint.hpp DESTINATION include)
install (FILES AuxiliaryConstraints/NonpositiveConstraint.hpp DESTINATION include)
install (FILES AuxiliaryConstraints/UnityConstraint.hpp DESTINATION include)
install (FILES AuxiliaryConstraintsProblem.hpp DESTINATION include)
install (FILES FeasibilityProblem.hpp DESTINATION include)
install (FILES InverseProblemSolver.hpp DESTINATION include)
install (FILES RangeProjectionSolver.hpp DESTINATION include)
install (FILES SolverFactory/SolverFactory.hpp DESTINATION include)
install (FILES SplitFeasibilitySolver.hpp DESTINATION include)
