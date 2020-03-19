function(add_tutorial target)

	set(targetname tutorial_${target})

	add_executable(
		${targetname}
		${ARGN}
	)

	target_include_directories(
		${targetname}
		PUBLIC
			"${CMAKE_CURRENT_SOURCE_DIR}"
			$<TARGET_PROPERTY:TransWORHP_core,INTERFACE_INCLUDE_DIRECTORIES>
	)

	if(${BUILD-OPTIONS_TransWORHP_GRAPHICS})
		target_link_libraries(
			${targetname}
			PUBLIC
				TransWORHP
		)
	else()
		target_link_libraries(
			${targetname}
			PUBLIC
				TransWORHPmini
		)
	endif()

	target_compile_options(
		${targetname}
		PUBLIC
			-std=c++11
		PRIVATE
#			-Wall
#			-Wextra
#			-pedantic
	)

	set_target_properties(${targetname} PROPERTIES OUTPUT_NAME ${target})

endfunction()

add_tutorial(spline0 spline0.cpp)
add_tutorial(spline1 spline1.cpp)
add_tutorial(spline2 spline2.cpp)
add_tutorial(spline3 spline3.cpp)
add_tutorial(spline4 spline4.cpp)
add_tutorial(spline5 spline5.cpp)
add_tutorial(spline_load spline_load.cpp)
add_tutorial(spline_ruku spline_ruku.cpp)
#add_tutorial(spline_ad spline_ad.cpp)
add_tutorial(spline_int spline_int.cpp)
add_tutorial(spline_zen spline_zen.cpp)
add_tutorial(spline_phase spline_phase.cpp)

