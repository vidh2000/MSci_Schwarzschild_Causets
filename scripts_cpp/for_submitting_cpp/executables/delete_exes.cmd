@ECHO OFF
:start

cd scripts_cpp

ECHO LISTING EXECUTABLE FILES
dir /b /s *.exe | findstr /e .exe 

SET choice=
SET /p choice=DO YOU WANT TO DELETE THESE [y/n]: 
IF NOT '%choice%'=='' SET choice=%choice:~0,1%
IF '%choice%'=='Y' GOTO yes
IF '%choice%'=='y' GOTO yes
IF '%choice%'=='N' GOTO no
IF '%choice%'=='n' GOTO no
IF '%choice%'=='' GOTO no
ECHO "%choice%" is not valid
ECHO.
GOTO start

:no
PAUSE
EXIT /b

:yes
del /s *.exe
PAUSE
EXIT /b

