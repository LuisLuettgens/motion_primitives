
rmdir /S /Q TransWORHP_Example\x64

del TransWORHP_Example.zip
7z.exe a TransWORHP_Example.zip TransWORHP_Example\ -x!TransWORHP_Example\TransWORHP_Example.sdf -x!TransWORHP_Example\TransWORHP_Example.v11.suo -x!TransWORHP_Example\TransWORHP_Example.suo -x!"TransWORHP_Example\Visual Studio\*\x64\*" -x!"TransWORHP_Example\Visual Studio\*\Release\*" -x!"TransWORHP_Example\Visual Studio\*\x64" -x!"TransWORHP_Example\Visual Studio\*\Release"


