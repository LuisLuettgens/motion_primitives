# Install script for directory: /home/luis/Backup2018Seafile/galileonautic/external/transworhp/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libTransWORHP.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libTransWORHP.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libTransWORHP.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/libTransWORHP.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libTransWORHP.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libTransWORHP.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libTransWORHP.so"
         OLD_RPATH "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/include/worhp/../../lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libTransWORHP.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/TransWORHP/core" TYPE FILE FILES
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/TransWORHP.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/Viewer.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/butcher.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/timing.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/MagicTransWORHP.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/ExplTransWORHP.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/PmTransWORHP.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/LobattoPmTransWORHP.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/FullDisTransWORHP.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/diffstructure.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/newdouble.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/worhp_info.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/TWsharedmemory.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/TWparameter.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/TWcount.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/TWfolder.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/TWconsole.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/TWdebug.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/TWproblem.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/twstatus.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/core/TWGUIconfig.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/TransWORHP/base" TYPE FILE FILES
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/base/defines.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/base/color4.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/base/vektor.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/base/point.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/TransWORHP/gui" TYPE FILE FILES
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/gui/sdlframe.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/gui/sdlscreen.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/gui/sdlcursor.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/gui/sdlthread.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/TransWORHP/glbase" TYPE FILE FILES
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/glbase/smoothmovement.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/glbase/light.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/glbase/font.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/glbase/joystick.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/glbase/texture.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/glbase/viewport.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/glbase/globject.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/glbase/model.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/TransWORHP/toolbase" TYPE FILE FILES
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolmenu.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolmenukey.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/toolbase/toolstatus.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/TransWORHP/glplot" TYPE FILE FILES
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/glplot/plot.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/glplot/baseplot.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/glplot/threedplot.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/glplot/xopt_data.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/glplot/functions.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/TransWORHP/xmlio" TYPE FILE FILES
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/xmlio/src/xmlio/textout.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/xmlio/src/xmlio/xmlio.h"
    "/home/luis/Backup2018Seafile/galileonautic/external/transworhp/src/xmlio/src/xmlio/conversion.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/base/cmake_install.cmake")
  include("/home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/core/cmake_install.cmake")
  include("/home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/glbase/cmake_install.cmake")
  include("/home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/glplot/cmake_install.cmake")
  include("/home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/gui/cmake_install.cmake")
  include("/home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/imaging/cmake_install.cmake")
  include("/home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/tool/cmake_install.cmake")
  include("/home/luis/Backup2018Seafile/galileonautic/build/external/transworhp/src/toolbase/cmake_install.cmake")

endif()

