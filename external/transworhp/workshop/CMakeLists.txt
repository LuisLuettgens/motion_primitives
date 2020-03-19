function(add_workshop target)

	set(targetname workshop_${target})

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

	target_link_libraries(
		${targetname}
		PUBLIC
			TransWORHP
	)

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

add_workshop(laufkatze laufkatze.cpp laufkatze_base.cpp)
#add_workshop(laufkatze_regelung laufkatze_regelung.cpp laufkatze_base.cpp)
#add_workshop(laufkatze_mpc laufkatze_mpc.cpp laufkatze_base.cpp)
add_workshop(raketenwagen1 raketenwagen1.cpp)
add_workshop(raketenwagen2 raketenwagen2.cpp)
add_workshop(raketenwagen3 raketenwagen3.cpp)
add_workshop(raketenwagen4 raketenwagen4.cpp)
add_workshop(raketenwagen5 raketenwagen5.cpp)
add_workshop(raketenwagen_constraints raketenwagen_constraints.cpp)
add_workshop(raketenwagen_load raketenwagen_load.cpp)
#add_workshop(raketenwagen_ad raketenwagen_ad.cpp)
add_workshop(raketenwagen_phase raketenwagen_phase.cpp)

