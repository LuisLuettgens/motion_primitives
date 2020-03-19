# Minimal benoetigte CMake-Version wenn z.B. bestimmte CMake-Kommandos benutzt werden
cmake_minimum_required(VERSION 2.8.12)

# Für ältere C++-Compiler wie z. B. unter Debian 8:
set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

# Der Projektname
project(TransWORHP)

# TransWORHP Version setzen
set(${PROJECT_NAME}_VERSION_MAJOR 0)
set(${PROJECT_NAME}_VERSION_MINOR 5)
set(${PROJECT_NAME}_VERSION_PATCH 0)
set(${PROJECT_NAME}_VERSION_TWEAK 0)

execute_process(
	COMMAND
		git rev-parse --short HEAD
	WORKING_DIRECTORY
		"${CMAKE_CURRENT_SOURCE_DIR}"
	OUTPUT_VARIABLE
		${PROJECT_NAME}_VERSION_TWEAK
	RESULT_VARIABLE
		TWEAK_COMMAND_RESULT
)

find_program(CLANG_TIDY_EXE NAMES "clang-tidy" DOC "Path to clang-tidy executable")
if(NOT CLANG_TIDY_EXE)
	message(STATUS "clang-tidy not found.")
else()
	message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
	set(DO_CLANG_TIDY "${CLANG_TIDY_EXE}" "-checks=*")
endif()

if (NOT(${TWEAK_COMMAND_RESULT} EQUAL 0))
	message(FATAL_ERROR "Failed to generate tweak number")
else ()
	string(REPLACE "\n" "" ${PROJECT_NAME}_VERSION_TWEAK ${${PROJECT_NAME}_VERSION_TWEAK})
	set(${PROJECT_NAME}_VERSION "\"${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}.${${PROJECT_NAME}_VERSION_TWEAK}\"")
	configure_file(
		"${PROJECT_SOURCE_DIR}/src/transworhp_version.h.in"
		"${CMAKE_BINARY_DIR}/transworhp_version.h"
	)
endif ()

option(BUILD-OPTIONS_TransWORHP_GRAPHICS "Build TransWORHP with GUI" ON)
option(BUILD-OPTIONS_TransWORHP_EXAMPLES "Build TransWORHP with examples" ON)
option(BUILD-OPTIONS_TransWORHP_TUTORIALS "Build TransWORHP with tutorials" ON)
option(BUILD-OPTIONS_TransWORHP_WORKSHOP "Build TransWORHP with workshop folder" OFF)
option(BUILD-OPTIONS_TransWORHP_BETTS "Build TransWORHP with Betts-Testset" OFF)
option(BUILD-OPTIONS_TransWORHP_UNITTEST "Build TransWORHP UnitTests" OFF)
option(BUILD-OPTIONS_TransWORHP_BOOST "Build TransWORHP with boost" OFF)

if (${BUILD-OPTIONS_TransWORHP_GRAPHICS})
	set(TRANSWORHP_GRAPHICS 1)
endif ()
configure_file(
	"${PROJECT_SOURCE_DIR}/src/TWGUIconfig.h.in"
	"${PROJECT_SOURCE_DIR}/src/core/TWGUIconfig.h"
)

# Unit-Tests
if (${BUILD-OPTIONS_TransWORHP_UNITTEST})
	enable_testing()
	add_test(unittest unittest)
endif ()

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/include")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/src")
if (${BUILD-OPTIONS_TransWORHP_TUTORIALS})
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/tutorial")
endif ()
if (${BUILD-OPTIONS_TransWORHP_GRAPHICS})
if (${BUILD-OPTIONS_TransWORHP_EXAMPLES})
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/example")
endif()
if (${BUILD-OPTIONS_TransWORHP_WORKSHOP})
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/workshop")
endif ()
if (${BUILD-OPTIONS_TransWORHP_BETTS})
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Betts")
endif ()
endif ()

