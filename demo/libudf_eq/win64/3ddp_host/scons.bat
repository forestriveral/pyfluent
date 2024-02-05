@echo off
echo.
echo Copyright 1987-2023 ANSYS, Inc. All Rights Reserved.

set PYTHONHOME=%FLUENT_INC%\..\commonfiles\CPython\3_7\winx64\Release\python
set PYTHONPATH=%FLUENT_INC%\..\commonfiles\CPython\3_7\winx64\Release\python

echo Compiler and linker: Clang (builtin)
set FLUENT_UDF_COMPILER=clang
set FLUENT_UDF_CLANG=builtin
"%PYTHONPATH%\python.exe" "%PYTHONPATH%\Scripts\scons.exe" -s
