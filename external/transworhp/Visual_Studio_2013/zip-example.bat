:: Ausgabe-Ordner
set DDD=%date:~6,4%_%date:~3,2%_%date:~0,2%
set ZIP=TransWORHP_Example_VS2013_%DDD%.zip
set DIR=TransWORHP_Example

del %ZIP%

del %DIR%\debug32\*.lib
del %DIR%\debug32\*.ilk
del %DIR%\debug32\*.exp
del %DIR%\debug32\*.pdb
del %DIR%\release32\*.lib
del %DIR%\release32\*.ilk
del %DIR%\release32\*.exp
del %DIR%\release32\*.pdb
del %DIR%\debug64\*.lib
del %DIR%\debug64\*.ilk
del %DIR%\debug64\*.exp
del %DIR%\debug64\*.pdb
del %DIR%\release64\*.lib
del %DIR%\release64\*.ilk
del %DIR%\release64\*.exp
del %DIR%\release64\*.pdb

"C:\Program Files\7-Zip\7z.exe" a %ZIP% %DIR%\ -x!%DIR%\TransWORHP_Example.sdf -x!%DIR%\TransWORHP_Example.v11.suo -x!%DIR%\TransWORHP_Example.suo -x!"%DIR%\Visual Studio\*\x64\*" -x!"%DIR%\Visual Studio\*\Release\*" -x!"%DIR%\Visual Studio\*\Debug\*" -x!"%DIR%\Visual Studio\*\x64" -x!"%DIR%\Visual Studio\*\Release" -x!"%DIR%\Visual Studio\*\Debug"


