# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Minimizations 
	BregmanDistance.cpp
	BregmanProjectionFunctional.cpp 
	DualityMappings/DualityMapping.cpp
	DualityMappings/IllegalDualityMapping.cpp
	DualityMappings/SoftThresholdingOperator.cpp
	Elements/SpaceElement.cpp
	GeneralMinimizer.cpp
	LandweberMinimizer.cpp
	LandweberMinimizer_calculateMatchingTau.cpp
	LandweberMinimizer_calculateOptimalStepwidth.cpp
	Mappings/LinearMapping.cpp
	MinimizerFactory.cpp
	SequentialSubspaceMinimizer.cpp
	SequentialSubspaceMinimizer_HyperplaneProjection.cpp
	SequentialSubspaceMinimizerNoise.cpp
	SmoothnessModulus.cpp
	Spaces/NormedSpace.cpp
	)

install (TARGETS Minimizations DESTINATION bin)
install (FILES BregmanDistance.hpp DESTINATION include)
install (FILES BregmanProjectionFunctional.hpp DESTINATION include)
install (FILES DualityMappings/DualityMapping.hpp DESTINATION include)
install (FILES DualityMappings/IllegalDualityMapping.hpp DESTINATION include)
install (FILES DualityMappings/SoftThresholdingOperator.hpp DESTINATION include)

install (FILES Elements/SpaceElement.hpp DESTINATION include)

install (FILES Functions/FunctionMinimizer.hpp DESTINATION include)
install (FILES Functions/FunctionMinimizer_impl.hpp DESTINATION include)
install (FILES Functions/MinimizationFunctional.hpp DESTINATION include)
install (FILES Functions/SmoothnessFunctional.hpp DESTINATION include)
install (FILES GeneralMinimizer.hpp DESTINATION include)
install (FILES GeneralMinimizer_DefaultDualityMappings.hpp DESTINATION include)
install (FILES GeneralMinimizer_RegularizedL1Norm.hpp DESTINATION include)
install (FILES LandweberMinimizer.hpp DESTINATION include)
install (FILES MinimizationExceptions.hpp DESTINATION include)
install (FILES MinimizerFactory.hpp DESTINATION include)

install (FILES Mappings/LinearMapping.hpp DESTINATION include)
install (FILES Mappings/Mapping.hpp DESTINATION include)

install (FILES Norms/L1Norm.hpp DESTINATION include)
install (FILES Norms/LInfinityNorm.hpp DESTINATION include)
install (FILES Norms/LpNorm.hpp DESTINATION include)
install (FILES Norms/Norm.hpp DESTINATION include)
install (FILES Norms/NormExceptions.hpp DESTINATION include)

install (FILES SequentialSubspaceMinimizer.hpp DESTINATION include)
install (FILES SequentialSubspaceMinimizer_HyperplaneProjection.hpp DESTINATION include)
install (FILES SequentialSubspaceMinimizerNoise.hpp DESTINATION include)
install (FILES SmoothnessModulus.hpp DESTINATION include)

install (FILES Spaces/NormedSpace.hpp DESTINATION include)

install (FILES types.hpp DESTINATION include)
