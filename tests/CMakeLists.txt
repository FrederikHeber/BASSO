### TESTS ###
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
                  DEPENDS 
                  Basso 
                  MatrixFactorizer
                  ComputerTomography
                  RangeProjector
                  NoiseAdder
)

# tests for BASSO program
add_subdirectory (Basso)

# tests for MatrixFactorizer program
add_subdirectory (MatrixFactorizer)

# tests for ComputerTomography program
add_subdirectory (ComputerTomography)

# tests for RangeProjector program
add_subdirectory (RangeProjector)

# tests for NoiseAdder program
add_subdirectory (NoiseAdder)

# tests for Gravity program
add_subdirectory (Gravity)

# tests for NonnegativeMatrixFactors program
add_subdirectory (NonnegativeMatrixFactors)

# tests for RadonMatrixWriter program
add_subdirectory (RadonMatrixWriter)

# test ndiff only when desired
if (USE_NDIFF)
	add_subdirectory (ndiff)
endif (USE_NDIFF)

# tests for MatrixToPNG program
add_subdirectory (MatrixToPNG)

# tests for Python interface
add_subdirectory (Python)

# tests for Userguide examples
add_subdirectory (Userguide)
