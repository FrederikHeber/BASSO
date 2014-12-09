### TESTS ###

########################## IllegalDualityMappingUnitTest #######################
add_executable (IllegalDualityMappingUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	IllegalDualityMappingUnitTest.cpp) 
target_link_libraries (IllegalDualityMappingUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}")

add_test (IllegalDualityMappingUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Mappings/unittests/IllegalDualityMappingUnitTest)

########################## LpDualityMappingUnitTest #######################
add_executable (LpDualityMappingUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	LpDualityMappingUnitTest.cpp) 
target_link_libraries (LpDualityMappingUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}")

add_test (LpDualityMappingUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Mappings/unittests/LpDualityMappingUnitTest)

########################## L1DualityMappingUnitTest #######################
add_executable (L1DualityMappingUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	L1DualityMappingUnitTest.cpp) 
target_link_libraries (L1DualityMappingUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}")

add_test (L1DualityMappingUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Mappings/unittests/L1DualityMappingUnitTest)

########################## LInfinityDualityMappingUnitTest #######################
add_executable (LInfinityDualityMappingUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	LInfinityDualityMappingUnitTest.cpp) 
target_link_libraries (LInfinityDualityMappingUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}")

add_test (LInfinityDualityMappingUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Mappings/unittests/LInfinityDualityMappingUnitTest)

########################## SoftThresholdingMappingUnitTest #######################
add_executable (SoftThresholdingMappingUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	SoftThresholdingMappingUnitTest.cpp) 
target_link_libraries (SoftThresholdingMappingUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}")

add_test (SoftThresholdingMappingUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Mappings/unittests/SoftThresholdingMappingUnitTest)
