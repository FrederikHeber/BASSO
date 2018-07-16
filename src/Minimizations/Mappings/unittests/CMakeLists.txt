### TESTS ###

########################## DualityMappingFactoryUnitTest #######################
add_executable (DualityMappingFactoryUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	DualityMappingFactoryUnitTest.cpp) 
target_link_libraries (DualityMappingFactoryUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES})

add_test (DualityMappingFactoryUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Mappings/unittests/DualityMappingFactoryUnitTest)

########################## IllegalDualityMappingUnitTest #######################
add_executable (IllegalDualityMappingUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	IllegalDualityMappingUnitTest.cpp) 
target_link_libraries (IllegalDualityMappingUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES})

add_test (IllegalDualityMappingUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Mappings/unittests/IllegalDualityMappingUnitTest)

########################## LpDualityMappingUnitTest #######################
add_executable (LpDualityMappingUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	LpDualityMappingUnitTest.cpp) 
target_link_libraries (LpDualityMappingUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES})

add_test (LpDualityMappingUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Mappings/unittests/LpDualityMappingUnitTest)

########################## L1DualityMappingUnitTest #######################
add_executable (L1DualityMappingUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	L1DualityMappingUnitTest.cpp) 
target_link_libraries (L1DualityMappingUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES})

add_test (L1DualityMappingUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Mappings/unittests/L1DualityMappingUnitTest)

########################## LInfinityDualityMappingUnitTest #######################
add_executable (LInfinityDualityMappingUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	LInfinityDualityMappingUnitTest.cpp) 
target_link_libraries (LInfinityDualityMappingUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES})

add_test (LInfinityDualityMappingUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Mappings/unittests/LInfinityDualityMappingUnitTest)

########################## RelativeShrinkageMappingUnitTest #######################
add_executable (RelativeShrinkageMappingUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	RelativeShrinkageMappingUnitTest.cpp) 
target_link_libraries (RelativeShrinkageMappingUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES})

add_test (RelativeShrinkageMappingUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Mappings/unittests/RelativeShrinkageMappingUnitTest)

########################## LinearMappingUnitTest #######################
add_executable (LinearMappingUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	LinearMappingUnitTest.cpp) 
target_link_libraries (LinearMappingUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES})

add_test (LinearMappingUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Mappings/unittests/LinearMappingUnitTest)

########################## NonLinearMappingUnitTest #######################
add_executable (NonLinearMappingUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	NonLinearMappingUnitTest.cpp) 
target_link_libraries (NonLinearMappingUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES})

add_test (NonLinearMappingUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Mappings/unittests/NonLinearMappingUnitTest)

