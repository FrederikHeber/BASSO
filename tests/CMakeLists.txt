### TESTS ###
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
                  DEPENDS 
                  Basso 
                  MatrixFactorizer)

# tests for BASSO program
add_subdirectory (Basso)

# tests for MatrixFactorizer program
add_subdirectory (MatrixFactorizer)
