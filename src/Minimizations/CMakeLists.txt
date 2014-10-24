# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Minimizations 
	DualityMappings/IllegalDualityMapping.cpp
	DualityMappings/SoftThresholdingOperator.cpp
	Elements/SpaceElement.cpp
	Functions/BregmanDistance.cpp
	Functions/BregmanProjectionFunctional.cpp 
	Functions/HyperplaneProjection.cpp
	Functions/SmoothnessModulus.cpp
	Mappings/LinearMapping.cpp
	Mappings/L1DualityMapping.cpp
	Mappings/LpDualityMapping.cpp
	Mappings/LInfinityDualityMapping.cpp
	Minimizers/GeneralMinimizer.cpp
	Minimizers/LandweberMinimizer.cpp
	Minimizers/LandweberMinimizer_calculateMatchingTau.cpp
	Minimizers/LandweberMinimizer_calculateOptimalStepwidth.cpp
	MinimizerFactory.cpp
	Minimizers/SequentialSubspaceMinimizer.cpp
	Minimizers/SequentialSubspaceMinimizerNoise.cpp
	Spaces/NormedSpace.cpp
	)

install (TARGETS Minimizations DESTINATION bin)
install (FILES DualityMappings/IllegalDualityMapping.hpp DESTINATION include)
install (FILES DualityMappings/SoftThresholdingOperator.hpp DESTINATION include)

install (FILES Elements/SpaceElement.hpp DESTINATION include)

install (FILES Functions/BregmanDistance.hpp DESTINATION include)
install (FILES Functions/BregmanProjectionFunctional.hpp DESTINATION include)
install (FILES Functions/FunctionMinimizer.hpp DESTINATION include)
install (FILES Functions/FunctionMinimizer_impl.hpp DESTINATION include)
install (FILES Functions/HyperplaneProjection.hpp DESTINATION include)
install (FILES Functions/MinimizationFunctional.hpp DESTINATION include)
install (FILES Functions/ResidualFunctional.hpp DESTINATION include)
install (FILES Functions/SmoothnessFunctional.hpp DESTINATION include)
install (FILES Functions/SmoothnessModulus.hpp DESTINATION include)

install (FILES Mappings/LinearMapping.hpp DESTINATION include)
install (FILES Mappings/L1DualityMapping.hpp DESTINATION include)
install (FILES Mappings/LpDualityMapping.hpp DESTINATION include)
install (FILES Mappings/LInfinityDualityMapping.hpp DESTINATION include)
install (FILES Mappings/Mapping.hpp DESTINATION include)
install (FILES Mappings/PowerTypeDualityMapping.hpp DESTINATION include)

install (FILES Minimizers/GeneralMinimizer.hpp DESTINATION include)
install (FILES Minimizers/GeneralMinimizer_DefaultDualityMappings.hpp DESTINATION include)
install (FILES Minimizers/GeneralMinimizer_RegularizedL1Norm.hpp DESTINATION include)
install (FILES Minimizers/LandweberMinimizer.hpp DESTINATION include)
install (FILES Minimizers/SequentialSubspaceMinimizer.hpp DESTINATION include)
install (FILES Minimizers/SequentialSubspaceMinimizerNoise.hpp DESTINATION include)

install (FILES MinimizationExceptions.hpp DESTINATION include)
install (FILES MinimizerFactory.hpp DESTINATION include)

install (FILES Norms/IllegalNorm.hpp DESTINATION include)
install (FILES Norms/L1Norm.hpp DESTINATION include)
install (FILES Norms/LInfinityNorm.hpp DESTINATION include)
install (FILES Norms/LpNorm.hpp DESTINATION include)
install (FILES Norms/Norm.hpp DESTINATION include)
install (FILES Norms/NormExceptions.hpp DESTINATION include)
install (FILES Norms/RegularizedL1Norm.hpp DESTINATION include)

install (FILES Spaces/NormedSpace.hpp DESTINATION include)

install (FILES types.hpp DESTINATION include)
