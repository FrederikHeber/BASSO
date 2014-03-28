# first we add the executable that generates the table
add_executable (MakeTable MakeTable.cpp)

# add command to generate the source code
add_custom_command (
	OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Table.h
	COMMAND MakeTable ${CMAKE_CURRENT_BINARY_DIR}/Table.h
	DEPENDS MakeTable
)

# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Minimizations mysqrt.cpp ${CMAKE_CURRENT_BINARY_DIR}/Table.h)

install (TARGETS Minimizations DESTINATION bin)
install (FILES BregmanFunctional.hpp DESTINATION include)
install (FILES DualityMapping.hpp DESTINATION include)
install (FILES Minimizations.h DESTINATION include)
