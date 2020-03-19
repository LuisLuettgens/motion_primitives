function(add_example target)

	set(targetname example_${target})

	add_executable(
		${targetname}
		${ARGN}
	)

	target_include_directories(
		${targetname}
		PUBLIC
			"${CMAKE_CURRENT_SOURCE_DIR}"
			$<TARGET_PROPERTY:xmlio,INTERFACE_INCLUDE_DIRECTORIES>
			$<TARGET_PROPERTY:TransWORHP_glplot,INTERFACE_INCLUDE_DIRECTORIES>
			$<TARGET_PROPERTY:TransWORHP_core,INTERFACE_INCLUDE_DIRECTORIES>
			$<TARGET_PROPERTY:TransWORHP_base,INTERFACE_INCLUDE_DIRECTORIES>
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


# checked
add_example(chemie chemie.cpp)
add_example(spline spline.cpp)
add_example(erzentlader erzentlader.cpp)
add_example(oszillator oszillator.cpp)
#
add_example(knickstab knickstab.cpp)
add_example(knickstab_int knickstab_int.cpp)
#
add_example(reentry reentry.cpp)
add_example(reentry_neben reentry_neben.cpp)
#
add_example(rakete rakete.cpp)
add_example(rakete_rand rakete_rand.cpp)
add_example(rakete_rand_neben rakete_rand_neben.cpp)
add_example(rakete_neben rakete_neben.cpp)
#
#
## unchecked
add_example(alpRider3_meshRef alpRider3_meshRef.cpp)
add_example(chemie_meshRef chemie_meshRef.cpp)
add_example(laufkatze_meshRef laufkatze_meshRef.cpp laufkatze_base.cpp)
add_example(rakete_meshRef rakete_meshRef.cpp)
add_example(spline_meshRef spline_meshRef.cpp)
add_example(rayleigh_meshRef rayleighOCP_meshRef.cpp)

add_example(rayleigh rayleighOCP.cpp)

add_example(emergencyLander emergencyLander.cpp)
add_example(emergencyLanderCtrl emergencyLanderCtrl.cpp)
add_example(goddard goddard.cpp)
add_example(einparken einparken.cpp)
add_example(einparken_LKW einparken_LKW.cpp)
add_example(ausweichen ausweichen.cpp)
add_example(brachistochrone brachistochrone.cpp)
add_example(brachistochrone2 brachistochrone2.cpp)
add_example(brachistochrone3 brachistochrone3.cpp)
add_example(brachistochrone4 brachistochrone4.cpp)
add_example(hangGlider hangGlider.cpp)
add_example(fischen fischen.cpp)
add_example(eismann eismann.cpp)
#add_example(orbit orbit.cpp)
#
add_example(turbolader turbolader_parId.cpp)
#
add_example(oszillator_neben oszillator_neben.cpp)
#
#
add_example(laufkatze laufkatze.cpp laufkatze_base.cpp)
#add_example(laufkatze_ad laufkatze_ad.cpp laufkatze_base.cpp)
add_example(laufkatze_int laufkatze_int.cpp laufkatze_base.cpp)
#
#add_example(laufkatze_scale laufkatze_scale.cpp)
add_example(laufkatze_phase laufkatze_phase.cpp laufkatze_base.cpp)
#add_example(laufkatze_zen laufkatze_zen.cpp)
add_example(flugzeug flugzeug.cpp)
#add_example(laufkatze_regelung laufkatze_regelung.cpp laufkatze_base.cpp)
add_example(laufkatze_mpc laufkatze_mpc.cpp laufkatze_base.cpp)
#
add_example(spline_phase spline_phase.cpp)
#add_example(spline_phase_meshRef spline_phase_meshRef.cpp)
#
#add_example(constrspline constrspline.cpp)

add_example(roboter roboter.cpp)
#add_example(geoleo geoleo.cpp)
#add_example(magicdiff magicdiff.cpp)
#
add_example(showResult showResult.cpp)
add_example(showResult2 showResult2.cpp laufkatze_base.cpp)

add_example(patterson_errorCalc patterson_errorCalc.cpp)
add_example(spline_errorCalc spline_errorCalc.cpp)

add_example(patterson patterson.cpp)

#add_example(orbitRaising_pm orbitRaising_pm.cpp)
#
add_example(LetterA LetterA3D.cpp)

