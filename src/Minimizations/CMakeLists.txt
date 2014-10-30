# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Minimizations 
	Elements/SpaceElement.cpp
	Functions/BregmanDistance.cpp
	Functions/BregmanProjectionFunctional.cpp 
	Functions/HyperplaneProjection.cpp
	Functions/Minimizers/Minimizer_gsl.cpp
	Functions/Minimizers/Minimizer_nlopt.cpp
	Functions/ResidualFunctional.cpp
	Functions/SmoothnessModulus.cpp
	Functions/VectorProjection.cpp
	Functions/VectorProjection_BregmanDistanceToLine.cpp
	InverseProblems/InverseProblem.cpp
	InverseProblems/InverseProblemFactory.cpp
	Mappings/IllegalDualityMapping.cpp
	Mappings/LinearMapping.cpp
	Mappings/L1DualityMapping.cpp
	Mappings/LpDualityMapping.cpp
	Mappings/LInfinityDualityMapping.cpp
	Mappings/PowerTypeDualityMapping.cpp
	Mappings/PowerTypeDualityMappingFactory.cpp
	Mappings/SoftThresholdingMapping.cpp
	Minimizers/GeneralMinimizer.cpp
	Minimizers/LandweberMinimizer.cpp
	Minimizers/MinimizerFactory.cpp
	Minimizers/SequentialSubspaceMinimizer.cpp
	Minimizers/SequentialSubspaceMinimizer_IterationState.cpp
	Minimizers/SequentialSubspaceMinimizerNoise.cpp
	Minimizers/StepWidths/ConstantRegularizedL1NormStepWidth.cpp
	Minimizers/StepWidths/DynamicRegularizedL1NormStepWidth.cpp
	Minimizers/StepWidths/DetermineStepWidthFactory.cpp
	Minimizers/StepWidths/LandweberFixedStepWidth.cpp
	Minimizers/StepWidths/ResidualMinimizingStepwidth.cpp
	Norms/NormFactory.cpp
	Spaces/NormedSpaceFactory.cpp
	Spaces/NormedSpace.cpp
	)

install (TARGETS Minimizations DESTINATION bin)

install (FILES Elements/SpaceElement.hpp DESTINATION include)

install (FILES Functions/BregmanDistance.hpp DESTINATION include)
install (FILES Functions/BregmanProjectionFunctional.hpp DESTINATION include)
install (FILES Functions/HyperplaneProjection.hpp DESTINATION include)
install (FILES Functions/Minimizers/FunctionalMinimizer.hpp DESTINATION include)
install (FILES Functions/Minimizers/FunctionalMinimizer_impl.hpp DESTINATION include)
install (FILES Functions/Minimizers/MinimizationFunctional.hpp DESTINATION include)
install (FILES Functions/Minimizers/Minimizer.hpp DESTINATION include)
install (FILES Functions/Minimizers/Minimizer_gsl.hpp DESTINATION include)
install (FILES Functions/Minimizers/Minimizer_nlopt.hpp DESTINATION include)
install (FILES Functions/ResidualFunctional.hpp DESTINATION include)
install (FILES Functions/SmoothnessFunctional.hpp DESTINATION include)
install (FILES Functions/SmoothnessModulus.hpp DESTINATION include)
install (FILES Functions/VectorProjection.hpp DESTINATION include)
install (FILES Functions/VectorProjection_BregmanDistanceToLine.hpp DESTINATION include)

install (FILES InverseProblems/InverseProblem.hpp DESTINATION include)
install (FILES InverseProblems/InverseProblemFactory.hpp DESTINATION include)

install (FILES Mappings/IllegalDualityMapping.hpp DESTINATION include)
install (FILES Mappings/LinearMapping.hpp DESTINATION include)
install (FILES Mappings/L1DualityMapping.hpp DESTINATION include)
install (FILES Mappings/LpDualityMapping.hpp DESTINATION include)
install (FILES Mappings/LInfinityDualityMapping.hpp DESTINATION include)
install (FILES Mappings/Mapping.hpp DESTINATION include)
install (FILES Mappings/PowerTypeDualityMapping.hpp DESTINATION include)
install (FILES Mappings/PowerTypeDualityMappingFactory.hpp DESTINATION include)
install (FILES Mappings/SoftThresholdingMapping.hpp DESTINATION include)

install (FILES Minimizers/GeneralMinimizer.hpp DESTINATION include)
install (FILES Minimizers/GeneralMinimizer_DefaultDualityMappings.hpp DESTINATION include)
install (FILES Minimizers/GeneralMinimizer_RegularizedL1Norm.hpp DESTINATION include)
install (FILES Minimizers/LandweberMinimizer.hpp DESTINATION include)
install (FILES Minimizers/MinimizationExceptions.hpp DESTINATION include)
install (FILES Minimizers/MinimizerFactory.hpp DESTINATION include)
install (FILES Minimizers/SequentialSubspaceMinimizer.hpp DESTINATION include)
install (FILES Minimizers/SequentialSubspaceMinimizerNoise.hpp DESTINATION include)
install (FILES Minimizers/StepWidths/ConstantRegularizedL1NormStepWidth.hpp DESTINATION include)
install (FILES Minimizers/StepWidths/DetermineStepWidth.hpp DESTINATION include)
install (FILES Minimizers/StepWidths/DetermineStepWidthFactory.hpp DESTINATION include)
install (FILES Minimizers/StepWidths/DynamicRegularizedL1NormStepWidth.hpp DESTINATION include)
install (FILES Minimizers/StepWidths/LandweberFixedStepWidth.hpp DESTINATION include)
install (FILES Minimizers/StepWidths/ResidualMinimizingStepwidth.hpp DESTINATION include)

install (FILES Norms/IllegalNorm.hpp DESTINATION include)
install (FILES Norms/L1Norm.hpp DESTINATION include)
install (FILES Norms/LInfinityNorm.hpp DESTINATION include)
install (FILES Norms/LpNorm.hpp DESTINATION include)
install (FILES Norms/Norm.hpp DESTINATION include)
install (FILES Norms/NormExceptions.hpp DESTINATION include)
install (FILES Norms/NormFactory.hpp DESTINATION include)
install (FILES Norms/RegularizedL1Norm.hpp DESTINATION include)

install (FILES Spaces/NormedSpaceFactory.hpp DESTINATION include)
install (FILES Spaces/NormedSpace.hpp DESTINATION include)

install (FILES types.hpp DESTINATION include)
