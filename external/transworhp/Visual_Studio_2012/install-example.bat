
set DIRC=TransWORHP_Example
set INC=%DIRC%\include

robocopy ..\example_release\ TransWORHP\Release64 transworhp.xml worhp.xml

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
robocopy ..\src\xmlio\ %INC%\xmlio textout.h xmlio.h conversion.h



rmdir %DIRC%\tutorial
mkdir %DIRC%\tutorial
copy ..\tutorial\*.cpp %DIRC%\tutorial


rmdir %DIRC%\example
mkdir %DIRC%\example
copy ..\example\*.cpp %DIRC%\example
copy ..\example\*.c %DIRC%\example
copy ..\example\*.h %DIRC%\example

mkdir %DIRC%\lib\
mkdir %DIRC%\release\
rem copy TransWORHP\Release64\TransWORHP.lib %DIRC%\lib\
rem copy TransWORHP\Release64\TransWORHP.dll %DIRC%\release\
rem copy TransWORHP\Release64\TransWORHP_mini.lib %DIRC%\lib\
rem copy TransWORHP\Release64\TransWORHP_mini.dll %DIRC%\release\
copy TransWORHP\Release64\transworhp.xml %DIRC%\release\
copy TransWORHP\Release64\worhp.xml %DIRC%\release\
rem : copy TransWORHP\lib64\libworhp.lib %DIRC%\lib\
rem copy TransWORHP\Release64\libworhp.dll %DIRC%\release\

robocopy TransWORHP\Release64\ %DIRC%\release\ TransWORHP.dll TransWORHP_mini.dll libworhp.dll SDL2.dll libpng16.dll libiomp5md.dll

rem transworhp.xml worhp.xml flugzeug.xml laufkatze.xml

robocopy TransWORHP\Release64\ %DIRC%\lib\ TransWORHP.lib TransWORHP_mini.lib

robocopy TransWORHP\lib64\ %DIRC%\lib SDL2.lib SDL2main.lib

robocopy TransWORHP\include\SDL2 %DIRC%\include\SDL2 *

robocopy ..\example_release\modell %DIRC%\release\modell *
robocopy ..\example_release\ %DIRC%\release\ flugzeug.xml laufkatze.xml

