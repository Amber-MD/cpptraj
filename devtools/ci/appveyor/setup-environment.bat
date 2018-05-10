@echo on

rem We need printf exponents to have 2 digits like the C standard says,
rem or tests will fail. Set environment variable for mingw-w64:
rem https://github.com/Alexpux/mingw-w64/blob/master/mingw-w64-crt/stdio/mingw_pformat.c#L223
set PRINTF_EXPONENT_DIGITS=2

set HOME=.
set MSYSTEM=MINGW64
set "PATH=C:/msys64/usr/bin;C:/msys64/mingw64/bin;%PATH%"

set MINGWPREFIX=x86_64-w64-mingw32
set CC=%MINGWPREFIX%-gcc.exe
set CXX=%MINGWPREFIX%-g++.exe
set FC=%MINGWPREFIX%-gfortran.exe
sh -lc "pacman -S --noconfirm --needed mingw-w64-x86_64-openblas mingw-w64-x86_64-arpack mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-ncurses mingw-w64-x86_64-readline mingw64/mingw-w64-x86_64-netcdf diffutils"

if %USE_CMAKE equ 1 (
		
	rem path from here: https://www.appveyor.com/docs/build-environment/#mingw-msys-cygwin
	set MINGWDIR=C:\mingw-w64\x86_64-6.3.0-posix-seh-rt_v5-rev1
	
	rem add mingw32-make and makensis to the PATH
	set "PATH=%PATH%;%MINGWDIR%\bin;C:\Program Files (x86)\NSIS"
	
) 
