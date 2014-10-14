# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Minimizations 
	BregmanDistance.cpp
	BregmanProjectionFunctional.cpp 
	DualityMappings/DualityMapping.cpp
	DualityMappings/IllegalDualityMapping.cpp
	DualityMappings/SoftThresholdingOperator.cpp
	GeneralMinimizer.cpp
	LandweberMinimizer.cpp
	LandweberMinimizer_calculateMatchingTau.cpp
	LandweberMinimizer_calculateOptimalStepwidth.cpp
	MinimizerFactory.cpp
	SequentialSubspaceMinimizer.cpp
	SequentialSubspaceMinimizer_HyperplaneProjection.cpp
	SequentialSubspaceMinimizerNoise.cpp
	SmoothnessModulus.cpp
	)

install (TARGETS Minimizations DESTINATION bin)
install (FILES BregmanDistance.hpp DESTINATION include)
install (FILES BregmanProjectionFunctional.hpp DESTINATION include)
install (FILES DualityMappings/DualityMapping.hpp DESTINATION include)
install (FILES DualityMappings/IllegalDualityMapping.hpp DESTINATION include)
install (FILES DualityMappings/SoftThresholdingOperator.hpp DESTINATION include)
install (FILES Functions/FunctionMinimizer.hpp DESTINATION include)
install (FILES Functions/FunctionMinimizer_impl.hpp DESTINATION include)
install (FILES Functions/MinimizationFunctional.hpp DESTINATION include)
install (FILES GeneralMinimizer.hpp DESTINATION include)
install (FILES LandweberMinimizer.hpp DESTINATION include)
install (FILES MinimizationExceptions.hpp DESTINATION include)
install (FILES MinimizerFactory.hpp DESTINATION include)
install (FILES SequentialSubspaceMinimizer.hpp DESTINATION include)
install (FILES SequentialSubspaceMinimizer_HyperplaneProjection.hpp DESTINATION include)
install (FILES SequentialSubspaceMinimizerNoise.hpp DESTINATION include)
install (FILES SmoothnessModulus.hpp DESTINATION include)
