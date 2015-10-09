#### EXECUTABLES ####

# add src folders of executables
add_subdirectory (BassoBase)
add_subdirectory (ComputerTomographyBase)
add_subdirectory (Gravity)
add_subdirectory (MatrixFactorizerBase)
if (PNG_FOUND)
	add_subdirectory (MatrixToPNG)
endif (PNG_FOUND)
add_subdirectory (NoiseAdderBase)
add_subdirectory (RangeProjectorBase)
