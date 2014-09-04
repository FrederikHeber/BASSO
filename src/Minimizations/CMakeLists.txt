# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Minimizations 
	BregmanDistance.cpp
	BregmanProjectionFunctional.cpp 
	DualityMapping.cpp
	GeneralMinimizer.cpp
	LandweberMinimizer.cpp
	LandweberMinimizer_calculateMatchingTau.cpp
	LandweberMinimizer_calculateOptimalStepwidth.cpp
	MinimizerFactory.cpp
	SequentialSubspaceMinimizer.cpp
	SequentialSubspaceMinimizerNoise.cpp
	SmoothnessModulus.cpp)

install (TARGETS Minimizations DESTINATION bin)
install (FILES BregmanDistance.hpp DESTINATION include)
install (FILES BregmanProjectionFunctional.hpp DESTINATION include)
install (FILES DualityMapping.hpp DESTINATION include)
install (FILES GeneralMinimizer.hpp DESTINATION include)
install (FILES LandweberMinimizer.hpp DESTINATION include)
install (FILES MinimizationExceptions.hpp DESTINATION include)
install (FILES MinimizerFactory.hpp DESTINATION include)
install (FILES SequentialSubspaceMinimizer.hpp DESTINATION include)
install (FILES SequentialSubspaceMinimizerNoise.hpp DESTINATION include)
install (FILES SmoothnessModulus.hpp DESTINATION include)
