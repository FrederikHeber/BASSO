### TESTS ###

########################## BackprojectionMatrixUnitTest #######################
add_executable (BackprojectionMatrixUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	BackprojectionMatrixUnitTest.cpp) 
target_link_libraries (BackprojectionMatrixUnitTest 
	ComputerTomographyBase
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES}) 

add_test (BackprojectionMatrixUnitTest 
	${PROJECT_BINARY_DIR}/src/ComputerTomographyBase/DiscretizedRadon/unittests/BackprojectionMatrixUnitTest)

########################## DiscretizedRadonMatrixUnitTest #######################
add_executable (DiscretizedRadonMatrixUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	DiscretizedRadonMatrixUnitTest.cpp) 
target_link_libraries (DiscretizedRadonMatrixUnitTest 
	ComputerTomographyBase
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES}) 

add_test (DiscretizedRadonMatrixUnitTest 
	${PROJECT_BINARY_DIR}/src/ComputerTomographyBase/DiscretizedRadon/unittests/DiscretizedRadonMatrixUnitTest)

########################## point_tUnitTest #######################
add_executable (point_tUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	point_tUnitTest.cpp) 
target_link_libraries (point_tUnitTest 
	ComputerTomographyBase
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES}) 

add_test (point_tUnitTest 
	${PROJECT_BINARY_DIR}/src/ComputerTomographyBase/DiscretizedRadon/unittests/point_tUnitTest)

