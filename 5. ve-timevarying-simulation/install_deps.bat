@echo off
REM ============================================================================
REM  install_deps.bat - install every R package + the CmdStan backend this
REM  project needs. Run this ONCE before run_full.bat (re-running is safe).
REM
REM  Double-click, or from a terminal:  install_deps.bat
REM
REM  Optional:  set VE_FORCE=1 & install_deps.bat   (reinstall CRAN packages)
REM
REM  Note: M4 also needs the Rtools C++ toolchain. If it is missing, the
REM  installer prints the official CRAN download link (do NOT use winget for
REM  Rtools on this machine - the winget package is broken here).
REM ============================================================================

setlocal enabledelayedexpansion
cd /d "%~dp0"

REM --- locate Rscript.exe -----------------------------------------------------
set "RSCRIPT="
for /f "delims=" %%i in ('where Rscript.exe 2^>nul') do (
    if not defined RSCRIPT set "RSCRIPT=%%i"
)
if not defined RSCRIPT (
    for /f "delims=" %%d in ('dir /b /ad /o-n "%ProgramFiles%\R\R-*" 2^>nul') do (
        if not defined RSCRIPT if exist "%ProgramFiles%\R\%%d\bin\Rscript.exe" (
            set "RSCRIPT=%ProgramFiles%\R\%%d\bin\Rscript.exe"
        )
    )
)
if not defined RSCRIPT (
    for /f "delims=" %%d in ('dir /b /ad /o-n "%ProgramFiles(x86)%\R\R-*" 2^>nul') do (
        if not defined RSCRIPT if exist "%ProgramFiles(x86)%\R\%%d\bin\Rscript.exe" (
            set "RSCRIPT=%ProgramFiles(x86)%\R\%%d\bin\Rscript.exe"
        )
    )
)
if not defined RSCRIPT (
    echo [ERROR] Could not find Rscript.exe.
    echo         Install R from https://cran.r-project.org/ first, then re-run.
    pause
    exit /b 1
)

echo Using R: "%RSCRIPT%"
echo Installing project dependencies...
echo.

"%RSCRIPT%" --vanilla install_deps.R
set "RC=%ERRORLEVEL%"

echo.
if "%RC%"=="0" (
    echo Dependency setup finished. Review the verification table above.
) else (
    echo [ERROR] install_deps.R exited with code %RC%.
)

pause
exit /b %RC%
