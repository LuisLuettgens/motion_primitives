project(TransWORHP_tool)

add_library(
	${PROJECT_NAME} OBJECT
		toolregister.cpp
)

target_include_directories(
	${PROJECT_NAME}
	PUBLIC
		"${CMAKE_CURRENT_SOURCE_DIR}"
)

target_compile_options(
	${PROJECT_NAME}
	PUBLIC
		-std=c++11
	PRIVATE
		-Wall
		-Wextra
		-pedantic
)

set_target_properties(
	${PROJECT_NAME}
		PROPERTIES POSITION_INDEPENDENT_CODE ON
)

