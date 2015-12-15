# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Minimizations 
	Elements/ElementCreator.cpp
	Elements/LinearDependencyChecker.cpp
	Elements/SpaceElement.cpp
	Elements/SpaceElementFactory.cpp
	Elements/SpaceElementIO.cpp
	Functions/BregmanDistance.cpp
	Functions/BregmanProjectionFunctional.cpp 
	Functions/MetricProjectionFunctional.cpp 
	Functions/Minimizers/FunctionalMinimizerFactory.cpp
	Functions/Minimizers/Minimizer_gsl.cpp
	Functions/Minimizers/Minimizer_nlopt.cpp
	Functions/ResidualFunctional.cpp
	Functions/SmoothnessModulus.cpp
	Functions/VectorProjection.cpp
	Functions/VectorProjection_BregmanDistanceToLine.cpp
	InverseProblems/InverseProblem.cpp
	InverseProblems/InverseProblemFactory.cpp
	Mappings/DualityMappingFactory.cpp
	Mappings/IllegalDualityMapping.cpp
	Mappings/LinearMapping.cpp
	Mappings/LinearMappingFactory.cpp
	Mappings/SingularValueDecomposition.cpp
	Mappings/SingularValueDecomposition_impl.cpp
	Mappings/Specifics/L1DualityMapping.cpp
	Mappings/Specifics/LpDualityMapping.cpp
	Mappings/Specifics/LInfinityDualityMapping.cpp
	Mappings/Specifics/RelativeShrinkageMapping.cpp
	Minimizers/GeneralMinimizer.cpp
	Minimizers/LandweberMinimizer.cpp
	Minimizers/MinimizerFactory.cpp
	Minimizers/SequentialSubspaceMinimizer.cpp
	Minimizers/SequentialSubspaceMinimizer_IterationState.cpp
	Minimizers/SequentialSubspaceMinimizerNoise.cpp
	Minimizers/Searchspace/LastNSearchDirections.cpp
	Minimizers/Searchspace/NemirovskyDirection.cpp
	Minimizers/Searchspace/Searchspace.cpp
	Minimizers/Searchspace/SearchspaceFactory.cpp
	Minimizers/StepWidths/ConstantRegularizedL1NormStepWidth.cpp
	Minimizers/StepWidths/DynamicRegularizedL1NormStepWidth.cpp
	Minimizers/StepWidths/DetermineStepWidthFactory.cpp
	Minimizers/StepWidths/LandweberFixedStepWidth.cpp
	Minimizers/StepWidths/ResidualMinimizingStepwidth.cpp
	Minimizers/StoppingCriteria/StoppingCriterion.cpp
	Minimizers/StoppingCriteria/StoppingCriteriaFactory.cpp
	Norms/NormFactory.cpp
	Spaces/NormedSpaceFactory.cpp
	Spaces/NormedSpace.cpp
	)
	
target_link_libraries(Minimizations ${NLOPT_LIBRARIES})

install (TARGETS Minimizations DESTINATION bin)

install (FILES Elements/Eigen_matrix_serialization.hpp DESTINATION include)
install (FILES Elements/ElementCreator.hpp DESTINATION include)
install (FILES Elements/LinearDependencyChecker.hpp DESTINATION include)
install (FILES Elements/RepresentationAdvocate.hpp DESTINATION include)
install (FILES Elements/SpaceElement.hpp DESTINATION include)
install (FILES Elements/SpaceElementFactory.hpp DESTINATION include)
install (FILES Elements/SpaceElementIO.hpp DESTINATION include)
install (FILES Elements/VectorSetter.hpp DESTINATION include)

install (FILES Functions/BregmanDistance.hpp DESTINATION include)
install (FILES Functions/BregmanProjectionFunctional.hpp DESTINATION include)
install (FILES Functions/HyperplaneProjection.hpp DESTINATION include)
install (FILES Functions/MetricProjectionFunctional.hpp DESTINATION include)
install (FILES Functions/Minimizers/FunctionMinimizer.hpp DESTINATION include)
install (FILES Functions/Minimizers/FunctionalMinimizerFactory.hpp DESTINATION include)
install (FILES Functions/Minimizers/FunctionalMinimizerFactory_impl.hpp DESTINATION include)
install (FILES Functions/Minimizers/FunctionalMinimizer.hpp DESTINATION include)
install (FILES Functions/Minimizers/FunctionalMinimizer_impl.hpp DESTINATION include)
install (FILES Functions/Minimizers/FunctionalMinimizer_exactLinesearch.hpp DESTINATION include)
install (FILES Functions/Minimizers/FunctionalMinimizer_exactLinesearch_impl.hpp DESTINATION include)
install (FILES Functions/Minimizers/FunctionalMinimizer_inexactLinesearch.hpp DESTINATION include)
install (FILES Functions/Minimizers/FunctionalMinimizer_inexactLinesearch_impl.hpp DESTINATION include)
install (FILES Functions/Minimizers/FunctionalMinimizer_deprecated.hpp DESTINATION include)
install (FILES Functions/Minimizers/FunctionalMinimizer_deprecated_impl.hpp DESTINATION include)
install (FILES Functions/Minimizers/MinimizationFunctional.hpp DESTINATION include)
install (FILES Functions/Minimizers/Minimizer.hpp DESTINATION include)
install (FILES Functions/Minimizers/Minimizer_gsl.hpp DESTINATION include)
install (FILES Functions/Minimizers/Minimizer_nlopt.hpp DESTINATION include)
install (FILES Functions/Minimizers/MinimizerExceptions.hpp DESTINATION include)
install (FILES Functions/ResidualFunctional.hpp DESTINATION include)
install (FILES Functions/SmoothnessFunctional.hpp DESTINATION include)
install (FILES Functions/SmoothnessModulus.hpp DESTINATION include)
install (FILES Functions/VectorProjection.hpp DESTINATION include)
install (FILES Functions/VectorProjection_BregmanDistanceToLine.hpp DESTINATION include)

install (FILES InverseProblems/InverseProblem.hpp DESTINATION include)
install (FILES InverseProblems/InverseProblemFactory.hpp DESTINATION include)
install (FILES InverseProblems/QuickAccessReferences.hpp DESTINATION include)

install (FILES Mappings/DualityMappingFactory.hpp DESTINATION include)
install (FILES Mappings/IllegalDualityMapping.hpp DESTINATION include)
install (FILES Mappings/LinearMapping.hpp DESTINATION include)
install (FILES Mappings/LinearMappingFactory.hpp DESTINATION include)
install (FILES Mappings/Mapping.hpp DESTINATION include)
install (FILES Mappings/MappingExceptions.hpp DESTINATION include)
install (FILES Mappings/PowerTypeDualityMapping.hpp DESTINATION include)
install (FILES Mappings/SingularValueDecomposition.hpp DESTINATION include)
install (FILES Mappings/SingularValueDecomposition_impl.hpp DESTINATION include)
install (FILES Mappings/Specifics/L1DualityMapping.hpp DESTINATION include)
install (FILES Mappings/Specifics/LpDualityMapping.hpp DESTINATION include)
install (FILES Mappings/Specifics/LInfinityDualityMapping.hpp DESTINATION include)
install (FILES Mappings/Specifics/RelativeShrinkageMapping.hpp DESTINATION include)

install (FILES Minimizers/GeneralMinimizer.hpp DESTINATION include)
install (FILES Minimizers/LandweberMinimizer.hpp DESTINATION include)
install (FILES Minimizers/MinimizationExceptions.hpp DESTINATION include)
install (FILES Minimizers/MinimizerFactory.hpp DESTINATION include)
install (FILES Minimizers/SequentialSubspaceMinimizer.hpp DESTINATION include)
install (FILES Minimizers/SequentialSubspaceMinimizerNoise.hpp DESTINATION include)
install (FILES Minimizers/Searchspace/LastNSearchDirections.hpp DESTINATION include)
install (FILES Minimizers/Searchspace/NemirovskyDirection.hpp DESTINATION include)
install (FILES Minimizers/Searchspace/Searchspace.hpp DESTINATION include)
install (FILES Minimizers/Searchspace/SearchspaceFactory.hpp DESTINATION include)
install (FILES Minimizers/StepWidths/ConstantRegularizedL1NormStepWidth.hpp DESTINATION include)
install (FILES Minimizers/StepWidths/DetermineStepWidth.hpp DESTINATION include)
install (FILES Minimizers/StepWidths/DetermineStepWidthFactory.hpp DESTINATION include)
install (FILES Minimizers/StepWidths/DynamicRegularizedL1NormStepWidth.hpp DESTINATION include)
install (FILES Minimizers/StepWidths/LandweberFixedStepWidth.hpp DESTINATION include)
install (FILES Minimizers/StepWidths/ResidualMinimizingStepwidth.hpp DESTINATION include)
install (FILES Minimizers/StoppingCriteria/CheckIterationCount.hpp DESTINATION include)
install (FILES Minimizers/StoppingCriteria/CheckRelativeChangeResiduum.hpp DESTINATION include)
install (FILES Minimizers/StoppingCriteria/CheckRelativeResiduum.hpp DESTINATION include)
install (FILES Minimizers/StoppingCriteria/CheckResiduum.hpp DESTINATION include)
install (FILES Minimizers/StoppingCriteria/CheckWalltime.hpp DESTINATION include)
install (FILES Minimizers/StoppingCriteria/StoppingArguments.hpp DESTINATION include)
install (FILES Minimizers/StoppingCriteria/StoppingCriterion.hpp DESTINATION include)
install (FILES Minimizers/StoppingCriteria/StoppingCriterion_impl.hpp DESTINATION include)
install (FILES Minimizers/StoppingCriteria/StoppingCriteriaFactory.hpp DESTINATION include)

install (FILES Norms/Specifics/DualRegularizedL1Norm.hpp DESTINATION include)
install (FILES Norms/Specifics/RegularizedL1Norm.hpp DESTINATION include)
install (FILES Norms/IllegalNorm.hpp DESTINATION include)
install (FILES Norms/L1Norm.hpp DESTINATION include)
install (FILES Norms/LInfinityNorm.hpp DESTINATION include)
install (FILES Norms/LpNorm.hpp DESTINATION include)
install (FILES Norms/Norm.hpp DESTINATION include)
install (FILES Norms/NormExceptions.hpp DESTINATION include)
install (FILES Norms/NormFactory.hpp DESTINATION include)

install (FILES Spaces/NormedDualSpace.hpp DESTINATION include)
install (FILES Spaces/NormedSpaceFactory.hpp DESTINATION include)
install (FILES Spaces/NormedSpace.hpp DESTINATION include)
install (FILES Spaces/OperationCountMap.hpp DESTINATION include)
install (FILES Spaces/VectorSpaceOperationCounts.hpp DESTINATION include)

install (FILES types.hpp DESTINATION include)
