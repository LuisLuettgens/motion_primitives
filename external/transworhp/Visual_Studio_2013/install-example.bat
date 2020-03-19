
set DIRC=TransWORHP_Example
set INC=%DIRC%\include

robocopy ..\example_release\ TransWORHP\Release64 transworhp.xml
REM worhp.xml

robocopy ..\src\ %INC%\ defines.h


rmdir %INC%\worhp
mkdir %INC%\worhp
copy TransWORHP\include\worhp\* %INC%\worhp

rmdir %INC%\core
mkdir %INC%\core
robocopy ..\src\core\ %INC%\core\ C_cs0.h defines.h TWfolder.h TWcount.h TWParameter.h TWConsole.h TransWORHP.h Viewer.h butcher.h diffstructure.h newdouble.h worhp_info.h transworhpsharedmem.h ExplTransWorhp.h timing.h MagicTransWorhp.h

rmdir %INC%\gui
mkdir %INC%\gui
robocopy ..\src\gui\ %INC%\gui\ sdlframe.h sdlscreen.h sdlcursor.h sdlthread.h

rmdir %INC%\base
mkdir %INC%\base
robocopy ..\src\base\ %INC%\base\ color4.h vektor.h point.h status.h defines.h

rmdir %INC%\glbase
mkdir %INC%\glbase
robocopy ..\src\glbase\ %INC%\glbase\ smoothmovement.h light.h font.h joystick.h texture.h viewport.h globject.h model.h

rmdir %INC%\gui
mkdir %INC%\gui

rmdir %INC%\glplot
mkdir %INC%\glplot
robocopy ..\src\glplot\ %INC%\glplot\ plot.h xopt_data.h functions.h

rmdir %INC%\toolbase
mkdir %INC%\toolbase
robocopy ..\src\toolbase\ %INC%\toolbase\ toolmenu.h toolmenukey.h toolstatus.h

rmdir %INC%\xmlio
mkdir %INC%\xmlio
robocopy ..\..\xmlio\src\xmlio\ %INC%\xmlio textout.h xmlio.h conversion.h



rmdir %DIRC%\tutorial
mkdir %DIRC%\tutorial
copy ..\tutorial\*.cpp %DIRC%\tutorial


rmdir %DIRC%\example
mkdir %DIRC%\example
copy ..\example\*.cpp %DIRC%\example
copy ..\example\*.c %DIRC%\example
copy ..\example\*.h %DIRC%\example

mkdir %DIRC%\lib32\ 
mkdir %DIRC%\lib64\ 
mkdir %DIRC%\release32\
mkdir %DIRC%\release64\
mkdir %DIRC%\debug32\
mkdir %DIRC%\debug64\

copy TransWORHP\Release64\transworhp.xml %DIRC%\release32\
copy TransWORHP\Release64\worhp.xml %DIRC%\release32\
copy TransWORHP\Release64\transworhp.xml %DIRC%\release64\
copy TransWORHP\Release64\worhp.xml %DIRC%\release64\

copy TransWORHP\Release64\transworhp.xml %DIRC%\debug32\
copy TransWORHP\Release64\worhp.xml %DIRC%\debug32\
copy TransWORHP\Release64\transworhp.xml %DIRC%\debug64\
copy TransWORHP\Release64\worhp.xml %DIRC%\debug64\


robocopy TransWORHP\Release64\ %DIRC%\release64\ TransWORHP.dll TransWORHP_mini.dll libworhp.dll SDL2.dll libiomp5md.dll
robocopy TransWORHP\Debug64\ %DIRC%\debug64\ TransWORHP.dll TransWORHP_mini_d.dll libworhp.dll SDL2.dll libiomp5md.dll

robocopy TransWORHP\Release\ %DIRC%\release32\ TransWORHP.dll TransWORHP_mini.dll libworhp.dll SDL2.dll libiomp5md.dll
robocopy TransWORHP\Debug\ %DIRC%\debug32\ TransWORHP.dll TransWORHP_mini_d.dll libworhp.dll SDL2.dll libiomp5md.dll

rem transworhp.xml worhp.xml flugzeug.xml laufkatze.xml

robocopy TransWORHP\Release64\ %DIRC%\lib64\ TransWORHP.lib TransWORHP_mini.lib
robocopy TransWORHP\Release\ %DIRC%\lib32\ TransWORHP.lib TransWORHP_mini.lib

robocopy TransWORHP\lib64\ %DIRC%\lib64 SDL2.lib SDL2main.lib
robocopy TransWORHP\lib32\ %DIRC%\lib32 SDL2.lib SDL2main.lib

robocopy TransWORHP\include\SDL2 %DIRC%\include\SDL2 *

robocopy ..\example_release\modell %DIRC%\release64\modell *
robocopy ..\example_release\ %DIRC%\release64\ flugzeug.xml laufkatze.xml

