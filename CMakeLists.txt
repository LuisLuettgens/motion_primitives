# Der Projektname
project(TransWORHP-Gnautic)

# Minimal benötigte CMake-Version wenn z.B. bestimmte CMake-Kommandos benutzt werden
cmake_minimum_required(VERSION 2.8)
set( CMAKE_VERBOSE_MAKEFILE on )


# Zusätzliche include-Verzeichnisse

IF( CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
ELSE()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
ENDIF()
# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# [TRANSWORHP]
set (BUILD-OPTIONS_TransWORHP_EXAMPLES OFF CACHE BOOL "don't make example files")
set (BUILD-OPTIONS_TransWORHP_TUTORIALS OFF CACHE BOOL "don't make tutorial files")
set (BUILD-OPTIONS_TransWORHP_WORKSHOP OFF CACHE BOOL "don't make workshop files")
set (BUILD-OPTIONS_TransWORHP_BETTS OFF CACHE BOOL "don't make examples from Betts")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/transworhp")

set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)


# set(targetname1 main_udp)
# set(targetname2 udp_server)
set(targetname1 main_motion_prims)
set(targetname2 main_arc)
set(targetname3 create_motion_prims)


# Die ausführbare Datei. Anstatt der Variablen können auch noch zusätzliche Dateien angegeben werden

file(GLOB_RECURSE DSL_SRC
    "${SRC_DIR}/external/dsl/lib/*.cc"
    "${SRC_DIR}/external/dsl/lib/*.c"
    "${SRC_DIR}/external/dsl/*.cpp"
)

# add_executable (${targetname1} ${SRC_DIR}/main_udp.cpp   ${SRC_DIR}/CourseCorrection.cpp ${SRC_DIR}/optimal_traj_Nships.cpp ${SRC_DIR}/Models/Model1v5MeckPomm.cpp ${SRC_DIR}/Models/Buggy_v1.cpp ${SRC_DIR}/shipviewer_base_Nships.cpp ${SRC_DIR}/geometric_constr.cpp ${SRC_DIR}/mapgrid_lla.cpp ${SRC_DIR}/utils/functions.cpp ${SRC_DIR}/DP/udp_client_server.cpp ${SRC_DIR}/DP/dp_base.cpp ${SRC_DIR}/DP/dp_0xA0.cpp ${SRC_DIR}/DP/dp_0xA1.cpp ${SRC_DIR}/DP/dp_0xA2.cpp ${SRC_DIR}/DP/dp_0xA3.cpp ${SRC_DIR}/DP/dp_0xA4.cpp ${SRC_DIR}/DP/dp_0xAC.cpp ${SRC_DIR}/DP/dp_0xAD.cpp ${SRC_DIR}/DP/dp_0xAE.cpp ${SRC_DIR}/DP/dp_0xAF.cpp ${SRC_DIR}/global.cpp ${DSL_SRC})
# add_executable (${targetname2} ${SRC_DIR}/udp_server.cpp ${SRC_DIR}/CourseCorrection.cpp ${SRC_DIR}/optimal_traj_Nships.cpp ${SRC_DIR}/Models/Model1v5MeckPomm.cpp ${SRC_DIR}/Models/Buggy_v1.cpp ${SRC_DIR}/shipviewer_base_Nships.cpp ${SRC_DIR}/geometric_constr.cpp ${SRC_DIR}/mapgrid_lla.cpp ${SRC_DIR}/utils/functions.cpp ${SRC_DIR}/DP/udp_client_server.cpp ${SRC_DIR}/DP/dp_base.cpp ${SRC_DIR}/DP/dp_0xA0.cpp ${SRC_DIR}/DP/dp_0xA1.cpp ${SRC_DIR}/DP/dp_0xA2.cpp ${SRC_DIR}/DP/dp_0xA3.cpp ${SRC_DIR}/DP/dp_0xA4.cpp ${SRC_DIR}/DP/dp_0xAC.cpp ${SRC_DIR}/DP/dp_0xAD.cpp ${SRC_DIR}/DP/dp_0xAE.cpp ${SRC_DIR}/DP/dp_0xAF.cpp ${SRC_DIR}/global.cpp ${DSL_SRC})
add_executable (${targetname1} ${SRC_DIR}/main_motion_prims.cpp   ${SRC_DIR}/CourseCorrection.cpp ${SRC_DIR}/optimal_traj_Nships.cpp ${SRC_DIR}/Models/Model1v5MeckPomm.cpp ${SRC_DIR}/Models/Buggy_v1.cpp ${SRC_DIR}/shipviewer_base_Nships.cpp ${SRC_DIR}/geometric_constr.cpp ${SRC_DIR}/mapgrid_lla.cpp ${SRC_DIR}/utils/functions.cpp ${SRC_DIR}/DP/udp_client_server.cpp ${SRC_DIR}/DP/dp_base.cpp ${SRC_DIR}/DP/dp_0xA0.cpp ${SRC_DIR}/DP/dp_0xA1.cpp ${SRC_DIR}/DP/dp_0xA2.cpp ${SRC_DIR}/DP/dp_0xA3.cpp ${SRC_DIR}/DP/dp_0xA4.cpp ${SRC_DIR}/DP/dp_0xAC.cpp ${SRC_DIR}/DP/dp_0xAD.cpp ${SRC_DIR}/DP/dp_0xAE.cpp ${SRC_DIR}/DP/dp_0xAF.cpp ${SRC_DIR}/global.cpp ${DSL_SRC})
add_executable (${targetname2} ${SRC_DIR}/main_arc.cpp   ${SRC_DIR}/CourseCorrection.cpp ${SRC_DIR}/optimal_traj_Nships.cpp ${SRC_DIR}/Models/Model1v5MeckPomm.cpp ${SRC_DIR}/Models/Buggy_v1.cpp ${SRC_DIR}/shipviewer_base_Nships.cpp ${SRC_DIR}/geometric_constr.cpp ${SRC_DIR}/mapgrid_lla.cpp ${SRC_DIR}/utils/functions.cpp ${SRC_DIR}/DP/udp_client_server.cpp ${SRC_DIR}/DP/dp_base.cpp ${SRC_DIR}/DP/dp_0xA0.cpp ${SRC_DIR}/DP/dp_0xA1.cpp ${SRC_DIR}/DP/dp_0xA2.cpp ${SRC_DIR}/DP/dp_0xA3.cpp ${SRC_DIR}/DP/dp_0xA4.cpp ${SRC_DIR}/DP/dp_0xAC.cpp ${SRC_DIR}/DP/dp_0xAD.cpp ${SRC_DIR}/DP/dp_0xAE.cpp ${SRC_DIR}/DP/dp_0xAF.cpp ${SRC_DIR}/global.cpp ${SRC_DIR}/external/dsl/dubins.cpp ${DSL_SRC})
add_executable (${targetname3} ${SRC_DIR}/create_motion_prims.cpp   ${SRC_DIR}/CourseCorrection.cpp ${SRC_DIR}/optimal_traj_Nships.cpp ${SRC_DIR}/Models/Model1v5MeckPomm.cpp ${SRC_DIR}/Models/Buggy_v1.cpp ${SRC_DIR}/shipviewer_base_Nships.cpp ${SRC_DIR}/geometric_constr.cpp ${SRC_DIR}/mapgrid_lla.cpp ${SRC_DIR}/utils/functions.cpp ${SRC_DIR}/DP/udp_client_server.cpp ${SRC_DIR}/DP/dp_base.cpp ${SRC_DIR}/DP/dp_0xA0.cpp ${SRC_DIR}/DP/dp_0xA1.cpp ${SRC_DIR}/DP/dp_0xA2.cpp ${SRC_DIR}/DP/dp_0xA3.cpp ${SRC_DIR}/DP/dp_0xA4.cpp ${SRC_DIR}/DP/dp_0xAC.cpp ${SRC_DIR}/DP/dp_0xAD.cpp ${SRC_DIR}/DP/dp_0xAE.cpp ${SRC_DIR}/DP/dp_0xAF.cpp ${SRC_DIR}/global.cpp ${DSL_SRC})

find_package(Eigen3 REQUIRED)
include_directories(${SRC_DIR}/external/eigen3/)

include_directories(${SRC_DIR}/external/dsl/lib)
include_directories(${SRC_DIR}/external/dsl/)
include_directories(/home/luis/worhp/include)

include_directories(${SRC_DIR})

target_include_directories(
		${targetname1}
		PUBLIC
            		"${CMAKE_CURRENT_SOURCE_DIR}"
            		$<TARGET_PROPERTY:xmlio,INTERFACE_INCLUDE_DIRECTORIES>
            		$<TARGET_PROPERTY:TransWORHP_glplot,INTERFACE_INCLUDE_DIRECTORIES>
            		$<TARGET_PROPERTY:TransWORHP_core,INTERFACE_INCLUDE_DIRECTORIES>
            		$<TARGET_PROPERTY:TransWORHP_base,INTERFACE_INCLUDE_DIRECTORIES>
    )

target_include_directories(
 		${targetname2}
 		PUBLIC
            		"${CMAKE_CURRENT_SOURCE_DIR}"
             		$<TARGET_PROPERTY:xmlio,INTERFACE_INCLUDE_DIRECTORIES>
            		$<TARGET_PROPERTY:TransWORHP_glplot,INTERFACE_INCLUDE_DIRECTORIES>
             		$<TARGET_PROPERTY:TransWORHP_core,INTERFACE_INCLUDE_DIRECTORIES>
             		$<TARGET_PROPERTY:TransWORHP_base,INTERFACE_INCLUDE_DIRECTORIES>
     )

target_include_directories(
 		${targetname3}
 		PUBLIC
            		"${CMAKE_CURRENT_SOURCE_DIR}"
             		$<TARGET_PROPERTY:xmlio,INTERFACE_INCLUDE_DIRECTORIES>
            		$<TARGET_PROPERTY:TransWORHP_glplot,INTERFACE_INCLUDE_DIRECTORIES>
             		$<TARGET_PROPERTY:TransWORHP_core,INTERFACE_INCLUDE_DIRECTORIES>
             		$<TARGET_PROPERTY:TransWORHP_base,INTERFACE_INCLUDE_DIRECTORIES>
     )


target_link_libraries(
        ${targetname1}
	PUBLIC
            TransWORHP
    )
target_link_libraries(
         ${targetname2}
	PUBLIC
           TransWORHP
    )
target_link_libraries(
         ${targetname3}
	PUBLIC
           TransWORHP
    )

#target_compile_options(
#        ${targetname}
#        PUBLIC
#            -std=c++11
#        PRIVATE
#            -Wall
#            -Wextra
#            -pedantic
#    )

set_target_properties(${targetname1} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set_target_properties(${targetname2} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set_target_properties(${targetname3} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

