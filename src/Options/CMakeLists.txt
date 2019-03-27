# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Options SHARED
	CommandLineOptions.cpp
	Options.cpp
	)

target_link_libraries(Options
	Minimizations
)

install (TARGETS Options DESTINATION lib)

install (FILES CommandLineOptions.hpp DESTINATION include)
install (FILES filesystem_path_serialization.hpp DESTINATION include)
install (FILES Options.hpp DESTINATION include)

