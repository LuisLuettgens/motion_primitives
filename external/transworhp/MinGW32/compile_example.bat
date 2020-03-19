set PATH=C:\mingw-w64\i686-4.9.2-posix-dwarf-rt_v4-rev2\mingw32\bin
set GCC=g++

set CFLAGS=-I../src/ -I../src/core -I"..\Visual Studio 2012\TransWORHP\include" -I../src/xmlio -DMINGW -m32

%GCC% %CFLAGS% -DNOGRAPHICS ../tutorial/spline4.cpp  -L. -ltransworhp-mini -o spline4mini
%GCC% %CFLAGS% ../tutorial/spline4.cpp -L. -ltransworhp -o spline4 -lmingw32 libSDL2main.a libSDL2.dll.a

