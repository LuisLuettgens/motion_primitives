project(TransWORHP_toolbase)

add_library(
	${PROJECT_NAME} OBJECT
		tool.cpp 
		toolstatus.cpp
		toolmenu.cpp
		toolrect.cpp 
		toolmenuentry.cpp
		toolmenukey.cpp
		toolwindow.cpp
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

