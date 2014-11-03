# add the executable
add_executable (Basso Basso.cpp)
target_link_libraries ( Basso 
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${Boost_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 


# add the install targets
install (TARGETS Basso DESTINATION bin)

# add test program to check that eigen-headers are found and compile
add_executable (EigenTest EigenTest.cpp)
target_link_libraries (EigenTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	${LIBRT_LIBRARIES})
	
