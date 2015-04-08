# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Database Database.cpp Table.cpp Table_Tuple.cpp)

install (FILES Database.hpp DESTINATION include)
install (FILES Database_impl.hpp DESTINATION include)
install (FILES Database_undef.hpp DESTINATION include)
install (FILES DefaultValue.hpp DESTINATION include)
install (FILES Table.hpp DESTINATION include)
install (FILES types.hpp DESTINATION include)
