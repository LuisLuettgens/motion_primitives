project(TransWORHP_glbase)

add_library(
	${PROJECT_NAME} OBJECT
		globject.cpp
		light.cpp
		smoothmovement.cpp
		viewport.cpp
		font.cpp
		model.cpp
		texture.cpp
		glfont.cpp
)

target_include_directories(
	${PROJECT_NAME}
	PUBLIC
		"${CMAKE_CURRENT_SOURCE_DIR}"
		$<TARGET_PROPERTY:xmlio,INTERFACE_INCLUDE_DIRECTORIES>
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

