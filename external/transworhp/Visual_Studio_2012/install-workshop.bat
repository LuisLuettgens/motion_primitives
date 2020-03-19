

set DIRC=TransWORHP_Workshop
set INC=%DIRC%\include

robocopy ..\workshop_release\ TransWORHP\release transworhp.xml worhp.xml


rmdir %INC%\worhp
mkdir %INC%\worhp
copy TransWORHP\include\worhp\* %INC%\worhp


rmdir %INC%\src
mkdir %INC%\src
robocopy ..\src\ %INC%\src\ C_cs0.h defines.h TWfolder.h TWcount.h TWParameter.h TWConsole.h TransWORHP.h Viewer.h butcher.h diffstructure.h newdouble.h worhp_info.h transworhpsharedmem.h 

rmdir %INC%\gui
mkdir %INC%\gui
mkdir %INC%\gui\base
mkdir %INC%\gui\glbase
mkdir %INC%\gui\gui
mkdir %INC%\gui\src_glplot
mkdir %INC%\gui\toolbase
robocopy ..\gui\base\ %INC%\gui\base\ color4.h vektor.h point.h
robocopy ..\gui\gui\ %INC%\gui\gui\ sdlframe.h sdlscreen.h sdlcursor.h sdlthread.h
robocopy ..\gui\glbase\ %INC%\gui\glbase\ smoothmovement.h light.h font.h joystick.h texture.h viewport.h globject.h model.h

robocopy ..\gui\toolbase\ %INC%\gui\toolbase\ toolmenu.h toolmenukey.h toolstatus.h

robocopy ..\gui\src_glplot\ %INC%\gui\src_glplot\ plot.h xopt_data.h functions.h

robocopy ..\gui\src_xmlio\ %INC%\ textout.h xmlio.h conversion.h




rmdir %DIRC%\workshop
mkdir %DIRC%\workshop
copy ..\workshop\*.cpp %DIRC%\workshop

mkdir %DIRC%\lib\
mkdir %DIRC%\release\
copy TransWORHP\Release\TransWORHP.lib %DIRC%\lib\
copy TransWORHP\Release\TransWORHP.dll %DIRC%\release\
copy TransWORHP\Release\TransWORHP_mini.lib %DIRC%\lib\
copy TransWORHP\Release\TransWORHP_mini.dll %DIRC%\release\
copy TransWORHP\Release\transworhp.xml %DIRC%\release\
copy TransWORHP\Release\worhp.xml %DIRC%\release\
copy TransWORHP\lib\libworhp.lib %DIRC%\lib\
copy TransWORHP\Release\libworhp.dll %DIRC%\release\
