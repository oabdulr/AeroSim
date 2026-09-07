@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0"

if /i "%~1"=="msvc" goto :find_msvc

rem Prefer a compiler already in PATH.
where gcc >nul 2>nul && goto :gcc
where cl  >nul 2>nul && goto :cl
where clang >nul 2>nul && goto :clang

:find_msvc
rem Locate MSVC through vswhere.
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist "%VSWHERE%" goto :nocc
for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VSDIR=%%i"
if not defined VSDIR goto :nocc
call "%VSDIR%\VC\Auxiliary\Build\vcvars64.bat" >nul
goto :cl

:gcc
echo Building with GCC...
gcc -O3 -march=native -ffast-math -funroll-loops -fopenmp -Wall -Wextra ^
    src\main.c src\sim.c src\airfoil.c src\render.c ^
    -o AeroSim.exe -lgdi32 -luser32 -mwindows -lm
if errorlevel 1 exit /b 1
gcc -O3 -march=native -ffast-math -funroll-loops -fopenmp -Wall -Wextra ^
    src3d\main3d.c src3d\sim3d.c src3d\wing3d.c src3d\gl3d.c src3d\gpu3d.c src3d\model3d.c ^
    -o AeroSim3D.exe -lopengl32 -lgdi32 -luser32 -lcomdlg32 -mwindows -lm
if errorlevel 1 exit /b 1
goto :done

:clang
echo Building with clang...
clang -O3 -ffast-math -Wall -Wextra ^
    src\main.c src\sim.c src\airfoil.c src\render.c ^
    -o AeroSim.exe -lgdi32 -luser32 -Wl,/SUBSYSTEM:WINDOWS
if errorlevel 1 exit /b 1
clang -O3 -ffast-math -Wall -Wextra ^
    src3d\main3d.c src3d\sim3d.c src3d\wing3d.c src3d\gl3d.c src3d\gpu3d.c src3d\model3d.c ^
    -o AeroSim3D.exe -lopengl32 -lgdi32 -luser32 -lcomdlg32 -Wl,/SUBSYSTEM:WINDOWS
if errorlevel 1 exit /b 1
goto :done

:cl
echo Building with MSVC...
cl /nologo /O2 /fp:fast /openmp /std:c11 /W3 /D_CRT_SECURE_NO_WARNINGS ^
   src\main.c src\sim.c src\airfoil.c src\render.c ^
   /Fe:AeroSim.exe /link /SUBSYSTEM:WINDOWS user32.lib gdi32.lib
if errorlevel 1 exit /b 1
cl /nologo /O2 /fp:fast /openmp /std:c11 /W3 /D_CRT_SECURE_NO_WARNINGS ^
   src3d\main3d.c src3d\sim3d.c src3d\wing3d.c src3d\gl3d.c src3d\gpu3d.c src3d\model3d.c ^
   /Fe:AeroSim3D.exe /link /SUBSYSTEM:WINDOWS user32.lib gdi32.lib opengl32.lib comdlg32.lib
if errorlevel 1 exit /b 1
del /q *.obj >nul 2>nul
goto :done

:nocc
echo ERROR: no C compiler found (tried gcc, cl, clang, vswhere).
exit /b 1

:done
echo.
echo Built AeroSim.exe and AeroSim3D.exe
exit /b 0
