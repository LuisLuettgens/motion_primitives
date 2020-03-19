find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../..)
add_executable(unittest main.cpp TWdimensionTest.cpp TWproblemTest.cpp)

target_include_directories(
	unittest
	PUBLIC
		"${CMAKE_CURRENT_SOURCE_DIR}"
		$<TARGET_PROPERTY:xmlio,INTERFACE_INCLUDE_DIRECTORIES>
		$<TARGET_PROPERTY:TransWORHP_glplot,INTERFACE_INCLUDE_DIRECTORIES>
		$<TARGET_PROPERTY:TransWORHP_core,INTERFACE_INCLUDE_DIRECTORIES>
		$<TARGET_PROPERTY:TransWORHP_base,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(
	unittest
	PUBLIC
		-DBOOST_ALL_DYN_LINK
)

target_compile_options(
	unittest
	PUBLIC
		-std=c++11
	PRIVATE
		-Wall
		-Wextra
		-pedantic
)

target_link_libraries(
	unittest
	PUBLIC
	${Boost_FILESYSTEM_LIBRARY}
	${Boost_SYSTEM_LIBRARY}
	${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
	TransWORHP
)
