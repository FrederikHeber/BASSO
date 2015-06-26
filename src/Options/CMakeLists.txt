# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Options 
	CommandLineOptions.cpp
	)

install (FILES CommandLineOptions.hpp DESTINATION include)

