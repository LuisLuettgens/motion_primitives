# Der Projektname
project(TransWORHP-XMLPARSER)

# Minimal benötigte CMake-Version wenn z.B. bestimmte CMake-Kommandos benutzt werden
cmake_minimum_required(VERSION 2.6)
#set( CMAKE_VERBOSE_MAKEFILE on )

# Benötigte oder optionale Pakete. Hier können z.B. Qt4 oder Boost stehen
#SET(wxWidgets_CONFIG_OPTIONS --version=2.8)
#FIND_PACKAGE(wxWidgets REQUIRED)


# Zusätzliche include-Verzeichnisse
include_directories(../include/ )
include_directories(../src/xmlio/ )
#include_directories(/localdata/opt/local/include/freetype2 )
include_directories(/localdata/opt/local/include/ )
include_directories(../src/glplot/ )
include_directories(../src/core/ )
include_directories(../src/ )



#SET(BUILD_SHARED_LIBS ON)

#link_directories(${PROJECT_SOURCE_DIR}/worhp/lib ${PROJECT_SOURCE_DIR}/glplot/lib)
link_directories(${PROJECT_SOURCE_DIR}/../lib )
link_directories(/localdata/opt/local/lib )




link_libraries(TransWORHP SDL2 GL glplot xmlio base glbase imaging gui tool toolbase worhp SDL2 imaging)

#if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
#   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -m32 -L/usr/lib32")
#   set(CMAKE_C_FLAGS "-I/usr/lib32/wx/include/gtk2-unicode-release-2.8 ${CMAKE_C_FLAGS} -m32 ")
#   set(CMAKE_CXX_FLAGS "-I/usr/lib32/wx/include/gtk2-unicode-release-2.8 ${CMAKE_CXX_FLAGS} -m32")
#endif( CMAKE_SIZEOF_VOID_P EQUAL 8 )

# Die ausführbare Datei. Anstatt der Variablen können auch noch zusätzliche Dateien angegeben werden

# checked
add_executable (xmlparser Main.cpp Variable.cpp Node.cpp XMLTransWorhp.cpp)

	set_target_properties(
	xmlparser
	PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/../xmlparser_release
)

