### TESTS ###

########################## BackprojectionMatrixUnitTest #######################
add_executable (BackprojectionMatrixUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	BackprojectionMatrixUnitTest.cpp) 
target_link_libraries (BackprojectionMatrixUnitTest 
	DiscretizedRadonBase
	${EXTRA_LIBS}
	${BOOST_CORE_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES}) 

add_test (BackprojectionMatrixUnitTest 
	${PROJECT_BINARY_DIR}/Examples/ComputerTomography/DiscretizedRadon/unittests/BackprojectionMatrixUnitTest)

########################## DiscretizedRadonMatrixUnitTest #######################
add_executable (DiscretizedRadonMatrixUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	DiscretizedRadonMatrixUnitTest.cpp) 
target_link_libraries (DiscretizedRadonMatrixUnitTest 
	DiscretizedRadonBase
	${EXTRA_LIBS}
	${BOOST_CORE_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES}) 

add_test (DiscretizedRadonMatrixUnitTest 
	${PROJECT_BINARY_DIR}/Examples/ComputerTomography/DiscretizedRadon/unittests/DiscretizedRadonMatrixUnitTest)

########################## point_tUnitTest #######################
add_executable (point_tUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	point_tUnitTest.cpp) 
target_link_libraries (point_tUnitTest 
	DiscretizedRadonBase
	${EXTRA_LIBS}
	${BOOST_CORE_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES}) 

add_test (point_tUnitTest 
	${PROJECT_BINARY_DIR}/Examples/ComputerTomography/DiscretizedRadon/unittests/point_tUnitTest)

