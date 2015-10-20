#### EXECUTABLES ####

# add src folders of executables
add_subdirectory (Basso)
add_subdirectory (ComputerTomography)
add_subdirectory (Gravity)
add_subdirectory (MatrixFactorizer)
if (PNG_FOUND)
	add_subdirectory (MatrixToPNG)
endif (PNG_FOUND)
add_subdirectory (NoiseAdder)
add_subdirectory (RangeProjector)
