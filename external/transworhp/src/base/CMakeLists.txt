project(TransWORHP_base)

add_library(
	${PROJECT_NAME} OBJECT
#		matrix.cpp # nicht genutzt
		point.cpp
		vectortools.cpp
		vektor.cpp
		color4.cpp
		exception.cpp
#		language.cpp # nicht genutzt
#		sharedmem.cpp # nicht genutzt
#		progressbar.cpp # nicht genutzt
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

if(CLANG_TIDY_EXE)
	set_target_properties(${PROJECT_NAME} PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
endif()
