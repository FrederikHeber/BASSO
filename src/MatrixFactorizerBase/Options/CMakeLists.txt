# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (MatrixFactorizerOptions 
	MatrixFactorizerOptions.cpp
	)

install (FILES MatrixFactorizerOptions.hpp DESTINATION include)
install (FILES filesystem_path_serialization.hpp DESTINATION include)
