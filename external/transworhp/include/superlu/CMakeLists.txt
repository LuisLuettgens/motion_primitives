project(superlu)

add_library(
	${PROJECT_NAME}
	SHARED
	IMPORTED
	GLOBAL
)

set_target_properties(
	${PROJECT_NAME}
	PROPERTIES
		IMPORTED_LOCATION
			"${CMAKE_SOURCE_DIR}/lib/libsuperlu.so"
		INTERFACE_INCLUDE_DIRECTORIES
			$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

