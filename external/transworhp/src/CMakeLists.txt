add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/base")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/core")

add_library(xmlio
	xmlio/src/xmlio/conversion.cpp
	xmlio/src/xmlio/textout.cpp
	xmlio/src/xmlio/xmlerror.cpp
	xmlio/src/xmlio/xmlnode.cpp
	xmlio/src/xmlio/xmlparser.cpp
)
target_include_directories(xmlio PUBLIC xmlio/src/xmlio)
set_target_properties(xmlio PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(${BUILD-OPTIONS_TransWORHP_GRAPHICS})
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/glbase")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/glplot")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/gui")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/imaging")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/tool")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/toolbase")
endif()

if(${BUILD-OPTIONS_TransWORHP_GRAPHICS})

find_package(GLEW REQUIRED)
find_package(OpenGL REQUIRED)

set(SDL2_LIBRARIES SDL2)
find_package(SDL2 QUIET)
string(STRIP ${SDL2_LIBRARIES} SDL2_LIBRARIES)

add_library(TransWORHP SHARED
	$<TARGET_OBJECTS:TransWORHP_base>
	$<TARGET_OBJECTS:TransWORHP_core>
	$<TARGET_OBJECTS:TransWORHP_glbase>
	$<TARGET_OBJECTS:TransWORHP_glplot>
	$<TARGET_OBJECTS:TransWORHP_gui>
	$<TARGET_OBJECTS:TransWORHP_imaging>
	$<TARGET_OBJECTS:TransWORHP_tool>
	$<TARGET_OBJECTS:TransWORHP_toolbase>
)

target_link_libraries(
	TransWORHP
	PUBLIC
		xmlio
		worhp
		${GLEW_LIBRARIES}
		${OPENGL_LIBRARIES}
		${SDL2_LIBRARIES}
		-fopenmp
)

else()

add_library(TransWORHPmini SHARED
	$<TARGET_OBJECTS:TransWORHP_base>
	$<TARGET_OBJECTS:TransWORHP_core>
)

target_link_libraries(
	TransWORHPmini
	PUBLIC
		xmlio
		worhp
		-fopenmp
)

endif()

if (${BUILD-OPTIONS_TransWORHP_UNITTEST})
	if (${BUILD-OPTIONS_TransWORHP_BOOST})
		add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/unittest")
	else()
		message("You have to activate boost!")
	endif()
endif ()

if(${BUILD-OPTIONS_TransWORHP_GRAPHICS})
	install(TARGETS TransWORHP LIBRARY DESTINATION lib)
else()
	install(TARGETS TransWORHPmini LIBRARY DESTINATION lib)
endif()

set(install_include_prefix include/TransWORHP)
install(FILES
	core/TransWORHP.h
	core/Viewer.h
	core/butcher.h
	core/timing.h
	core/MagicTransWORHP.h
	core/ExplTransWORHP.h
	core/PmTransWORHP.h
	core/LobattoPmTransWORHP.h
	core/FullDisTransWORHP.h
	core/diffstructure.h
	core/newdouble.h
	core/worhp_info.h
	core/TWsharedmemory.h
	core/TWparameter.h
	core/TWcount.h
	core/TWfolder.h
	core/TWconsole.h
	core/TWdebug.h
	core/TWproblem.h
	core/twstatus.h
	core/TWGUIconfig.h
		DESTINATION ${install_include_prefix}/core)

install(FILES
	base/defines.h
	base/color4.h
	base/vektor.h
	base/point.h
		DESTINATION ${install_include_prefix}/base)

install(FILES
	gui/sdlframe.h
	gui/sdlscreen.h
	gui/sdlcursor.h
	gui/sdlthread.h
		DESTINATION ${install_include_prefix}/gui)

install(FILES
	glbase/smoothmovement.h
	glbase/light.h
	glbase/font.h
	glbase/joystick.h
	glbase/texture.h
	glbase/viewport.h
	glbase/globject.h
	glbase/model.h
		DESTINATION ${install_include_prefix}/glbase)

install(FILES
	toolbase/toolmenu.h
	toolbase/toolmenukey.h
	toolbase/toolstatus.h
		DESTINATION ${install_include_prefix}/toolbase)

install(FILES
	glplot/plot.h
	glplot/baseplot.h
	glplot/threedplot.h
	glplot/xopt_data.h
	glplot/functions.h
		DESTINATION ${install_include_prefix}/glplot)

install(FILES
	xmlio/src/xmlio/textout.h
	xmlio/src/xmlio/xmlio.h
	xmlio/src/xmlio/conversion.h
		DESTINATION ${install_include_prefix}/xmlio)

SET(CPACK_GENERATOR "DEB")
set(CPACK_PACKAGE_VENDOR "Uni Bremen")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Matthias Rick <mrick@math.uni-bremen.de>")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "solver for optimal control problems")
set(CPACK_DEBIAN_PACKAGE_SECTION "math")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
#set(CPACK_DEBIAN_PACKAGE_DEPENDS "libsdl2-dev (>=2.0.1)")

set(CPACK_PACKAGE_VERSION_MAJOR "${${PROJECT_NAME}_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${${PROJECT_NAME}_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${${PROJECT_NAME}_VERSION_PATCH}")

include(CPack)
