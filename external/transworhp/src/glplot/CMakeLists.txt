project(TransWORHP_glplot)

add_library(
	${PROJECT_NAME} OBJECT
		dataplot.cpp
		matrixplot.cpp
		baseplot.cpp
		sparseplot.cpp
		threedplot.cpp
		xopt_data.cpp
		generalplot.cpp
		phaseplot.cpp
		tabularplot.cpp
		userplot.cpp
		xopt_eps.cpp
		punktPlot.cpp
		gitterPlot.cpp
		lambdaPlot.cpp
		adjPlot.cpp
)

target_include_directories(
	${PROJECT_NAME}
	PUBLIC
		"${CMAKE_CURRENT_SOURCE_DIR}"
		$<TARGET_PROPERTY:TransWORHP_core,INTERFACE_INCLUDE_DIRECTORIES>
		$<TARGET_PROPERTY:xmlio,INTERFACE_INCLUDE_DIRECTORIES>
		$<TARGET_PROPERTY:worhp,INTERFACE_INCLUDE_DIRECTORIES>
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

