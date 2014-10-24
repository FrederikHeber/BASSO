### TESTS ###

########################## LpDualityMappingUnitTest #######################
add_executable (LpDualityMappingUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	LpDualityMappingUnitTest.cpp) 
target_link_libraries (LpDualityMappingUnitTest 
	${EXTRA_LIBS}
	"${CPPUNIT_LIBRARY}")

add_test (LpDualityMappingUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/DualityMappings/unittests/LpDualityMappingUnitTest)

########################## IllegalDualityMappingUnitTest #######################
add_executable (IllegalDualityMappingUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	IllegalDualityMappingUnitTest.cpp) 
target_link_libraries (IllegalDualityMappingUnitTest 
	${EXTRA_LIBS}
	"${CPPUNIT_LIBRARY}")

add_test (IllegalDualityMappingUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/DualityMappings/unittests/IllegalDualityMappingUnitTest)

########################## SoftThresholdingOperatorUnitTest #######################
add_executable (SoftThresholdingOperatorUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	SoftThresholdingOperatorUnitTest.cpp) 
target_link_libraries (SoftThresholdingOperatorUnitTest 
	${EXTRA_LIBS}
	"${CPPUNIT_LIBRARY}")

add_test (SoftThresholdingOperatorUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/DualityMappings/unittests/SoftThresholdingOperatorUnitTest)
