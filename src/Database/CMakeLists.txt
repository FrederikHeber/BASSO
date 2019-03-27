# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Database SHARED
	SQLDatabase.cpp 
	Table.cpp 
	TableDirectory.cpp 
	Table_Tuple.cpp)

target_link_libraries(Database
	Minimizations
	Log
	${POCO_SQL_LIBRARIES}
)

install (TARGETS Database DESTINATION lib)

install (FILES Database.hpp DESTINATION include)
install (FILES DatabaseExceptions.hpp DESTINATION include)
install (FILES Database_impl.hpp DESTINATION include)
install (FILES Database_mock.hpp DESTINATION include)
install (FILES Database_undef.hpp DESTINATION include)
install (FILES DefaultValue.hpp DESTINATION include)
install (FILES SQLDatabase.hpp DESTINATION include)
install (FILES Table.hpp DESTINATION include)
install (FILES TableDirectory.hpp DESTINATION include)
install (FILES TableDirectoryDatabase.hpp DESTINATION include)
install (FILES TableDirectoryDatabase_mock.hpp DESTINATION include)
install (FILES Table_mock.hpp DESTINATION include)
install (FILES Tuple_mock.hpp DESTINATION include)
install (FILES types.hpp DESTINATION include)
