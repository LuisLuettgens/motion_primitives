set PATH=C:\mingw-w64\i686-4.9.2-posix-dwarf-rt_v4-rev2\mingw32\bin
set GCC=g++

set CFLAGS=-I../src/ -Iinclude -I../src/xmlio -DMINGW -m32

%GCC% %CFLAGS% -c ../src/core/C_cs0.cpp
%GCC% %CFLAGS% -c ../src/core/diffstructure.cpp
%GCC% %CFLAGS% -c ../src/core/newdouble.cpp
%GCC% %CFLAGS% -c ../src/core/TransWORHP.cpp
%GCC% %CFLAGS% -c ../src/core/ExplTransWORHP.cpp
%GCC% %CFLAGS% -c ../src/core/worhp_info.cpp
%GCC% %CFLAGS% -c ../src/core/butcher.cpp
%GCC% %CFLAGS% -c ../src/core/TWsharedmemory.cpp
%GCC% %CFLAGS% -c ../src/core/TWstrings.cpp
%GCC% %CFLAGS% -c ../src/core/TWparameter.cpp
%GCC% %CFLAGS% -c ../src/core/TWcount.cpp
%GCC% %CFLAGS% -c ../src/core/TWdebug.cpp
%GCC% %CFLAGS% -c ../src/core/TWspline.cpp
%GCC% %CFLAGS% -c ../src/core/TWfolder.cpp
%GCC% %CFLAGS% -c ../src/core/TWconsole.cpp
%GCC% %CFLAGS% -c ../src/core/MyStatus.cpp
%GCC% %CFLAGS% -c ../src/core/timing.cpp
%GCC% %CFLAGS% -c ../src/xmlio/conversion.cpp
%GCC% %CFLAGS% -c ../src/xmlio/xmlnode.cpp
%GCC% %CFLAGS% -c ../src/xmlio/xmlparser.cpp
%GCC% %CFLAGS% -c ../src/xmlio/xmlerror.cpp
%GCC% %CFLAGS% -c ../src/xmlio/textout.cpp
%GCC% %CFLAGS% -c ../src/base/point.cpp
%GCC% %CFLAGS% -c ../src/core/Viewer.cpp


%GCC% %CFLAGS% -c ../src/gui/sdlcursor.cpp
%GCC% %CFLAGS% -c ../src/gui/sdlframe.cpp
%GCC% %CFLAGS% -c ../src/gui/sdlscreen.cpp
%GCC% %CFLAGS% -c ../src/gui/sdlthread.cpp
%GCC% %CFLAGS% -c ../src/toolbase/toolmenu.cpp
%GCC% %CFLAGS% -c ../src/toolbase/toolwindow.cpp
%GCC% %CFLAGS% -c ../src/glplot/adjPlot.cpp
%GCC% %CFLAGS% -c ../src/glplot/dataplot.cpp
%GCC% %CFLAGS% -c ../src/glplot/generalplot.cpp
%GCC% %CFLAGS% -c ../src/glplot/lambdaPlot.cpp
%GCC% %CFLAGS% -c ../src/glplot/matrixplot.cpp
%GCC% %CFLAGS% -c ../src/glplot/gitterPlot.cpp
%GCC% %CFLAGS% -c ../src/glplot/phasePlot.cpp
%GCC% %CFLAGS% -c ../src/glplot/plot.cpp
%GCC% %CFLAGS% -c ../src/glplot/punktPlot.cpp
%GCC% %CFLAGS% -c ../src/glplot/sparseplot.cpp
%GCC% %CFLAGS% -c ../src/glplot/tabularplot.cpp
%GCC% %CFLAGS% -c ../src/glplot/threedplot.cpp
%GCC% %CFLAGS% -c ../src/glplot/userplot.cpp
%GCC% %CFLAGS% -c ../src/glplot/xopt_data.cpp
%GCC% %CFLAGS% -c ../src/glplot/xopt_eps.cpp
%GCC% %CFLAGS% -c ../src/base/color4.cpp
%GCC% %CFLAGS% -c ../src/glbase/font.cpp
%GCC% %CFLAGS% -c ../src/glbase/glfont.cpp
%GCC% %CFLAGS% -c ../src/glbase/globject.cpp
%GCC% %CFLAGS% -c ../src/imaging/imageformat.cpp
%GCC% %CFLAGS% -c ../src/imaging/imagewriter.cpp
%GCC% %CFLAGS% -c ../src/glbase/joystick.cpp
%GCC% %CFLAGS% -c ../src/base/language.cpp
%GCC% %CFLAGS% -c ../src/glbase/light.cpp
%GCC% %CFLAGS% -c ../src/base/matrix.cpp
%GCC% %CFLAGS% -c ../src/glbase/model.cpp

%GCC% %CFLAGS% -c ../src/base/point.cpp
%GCC% %CFLAGS% -c ../src/glbase/smoothmovement.cpp
%GCC% %CFLAGS% -c ../src/glbase/texture.cpp
%GCC% %CFLAGS% -c ../src/toolbase/tool.cpp
%GCC% %CFLAGS% -c ../src/toolbase/tooledit.cpp
%GCC% %CFLAGS% -c ../src/toolbase/toolfloating.cpp
%GCC% %CFLAGS% -c ../src/toolbase/toolmenuentry.cpp

%GCC% %CFLAGS% -c ../src/toolbase/toolmenukey.cpp
%GCC% %CFLAGS% -c ../src/toolbase/toolrect.cpp
%GCC% %CFLAGS% -c ../src/tool/toolregister.cpp
%GCC% %CFLAGS% -c ../src/toolbase/toolstatus.cpp
%GCC% %CFLAGS% -c ../src/base/vectortools.cpp
%GCC% %CFLAGS% -c ../src/base/vectortools2.cpp
%GCC% %CFLAGS% -c ../src/base/vektor.cpp
%GCC% %CFLAGS% -c ../src/glbase/viewport.cpp
%GCC% %CFLAGS% -c ../src/base/exception.cpp

%GCC% -shared  -o libtransworhp.dll *.o -Wl,--out-implib,libtransworhp.a -L. -llibworhp -lmingw32 libSDL2.dll.a libSDL2main.a -lopengl32

