set PATH=C:\mingw-w64\i686-4.9.2-posix-dwarf-rt_v4-rev2\mingw32\bin
set GCC=g++

set CFLAGS=-I../src/ -Iinclude -I../src/xmlio -DNOGRAPHICS -DMINGW -m32

del *.o
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

%GCC% -shared  -o libtransworhp-mini.dll *.o -Wl,--out-implib,libtransworhp-mini.a -L. -llibworhp

