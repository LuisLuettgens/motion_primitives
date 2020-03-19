project(TransWORHP_core)

if(NOT BUILD-OPTIONS_CORE_LGS_SOLVER_FOR_ROW)
	set(
		BUILD-OPTIONS_CORE_LGS_SOLVER_FOR_ROW
			"NONE"
		CACHE
			STRING
			"If the library should be build with X. Only for ROW-methode in explTW!"
		FORCE
	)
	set_property(
		CACHE
			BUILD-OPTIONS_CORE_LGS_SOLVER_FOR_ROW
		PROPERTY
			STRINGS
				"lapack"
				"SuperLU"
				"NONE"
	)
endif()

set(${PROJECT_NAME}_SRC
	diffstructure.cpp
	newdouble.cpp
	TransWORHP.cpp
	FullDisTransWORHP.cpp
#	MagicTransWORHP.cpp
	ExplTransWORHP.cpp
	worhp_info.cpp
	butcher.cpp
	TWsharedmemory.cpp
	TWstrings.cpp
	TWparameter.cpp
	TWcount.cpp
	TWdebug.cpp
	TWfolder.cpp
	TWconsole.cpp
	MyStatus.cpp
	timing.cpp
	TWspline.cpp
	rkf45.cpp
	PmTransWORHP.cpp
	LobattoPmTransWORHP.cpp
	GaussPmTransWORHP.cpp
	TWproblem.cpp
)

if(${BUILD-OPTIONS_TransWORHP_GRAPHICS})
	set(${PROJECT_NAME}_SRC ${${PROJECT_NAME}_SRC} Viewer.cpp)
endif()

add_library(
	${PROJECT_NAME} OBJECT
		${${PROJECT_NAME}_SRC}
)

target_include_directories(
	${PROJECT_NAME}
	PRIVATE
		"${CMAKE_BINARY_DIR}"
	PUBLIC
		"${CMAKE_CURRENT_SOURCE_DIR}"
		$<TARGET_PROPERTY:worhp,INTERFACE_INCLUDE_DIRECTORIES>
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
		-fopenmp
)

if (${BUILD-OPTIONS_TransWORHP_BOOST})
	target_compile_definitions(
		${PROJECT_NAME}
		PRIVATE
			-DTW_WITH_BOOST
	)
endif()

set_target_properties(
	${PROJECT_NAME}
		PROPERTIES POSITION_INDEPENDENT_CODE ON
)

if(${BUILD-OPTIONS_CORE_LGS_SOLVER_FOR_ROW} STREQUAL "lapack")
	target_compile_definitions(
		${PROJECT_NAME}
		PRIVATE
			-DWITH_LAPACK
	)
	target_link_libraries(
		${PROJECT_NAME}
		PRIVATE
			lapack
			blas
	)
elseif(${BUILD-OPTIONS_CORE_LGS_SOLVER_FOR_ROW} STREQUAL "SuperLU")
	target_include_directories(
		${PROJECT_NAME}
		PRIVATE
			$<TARGET_PROPERTY:superlu,INTERFACE_INCLUDE_DIRECTORIES>
	)
	target_compile_definitions(
		${PROJECT_NAME}
		PRIVATE
			-DTW_WITH_SUPERLU
	)
	target_link_libraries(
		${PROJECT_NAME}
		PRIVATE
			superlu
			blas
	)
elseif(${BUILD-OPTIONS_CORE_LGS_SOLVER_FOR_ROW} STREQUAL "NONE")
	# Everything is ok!
else()
	message(FATAL_ERROR "Not supported Math backend")
endif()

