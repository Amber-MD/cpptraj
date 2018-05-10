@echo on

rem add makensis to the PATH
set "PATH=%PATH%;C:\Program Files (x86)\NSIS"

set MINGWPREFIX=x86_64-w64-mingw32

if %BUILD_TYPE% equ cmake-vs (
	
	rem download prebuilt NetCDF and FFTW (hosted by Jamie Smith)
	powershell -Command "(new-object System.Net.WebClient).DownloadFile('https://app.box.com/shared/static/cc04s69672igsbfb4n3efbmv4hxlycu2.7z','cpptraj-msvc-prebuilts.7z')
	7z x cpptraj-msvc-prebuilts.7z
	set PREBUILTS_DIR=%cd%\cpptraj-msvc-prebuilts
	
	rem add prebuilt DLLs and Dumpbin to the PATH
	set "PATH=%PATH%;%PREBUILTS_DIR%\bin;C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin"
	
) else (
	rem set up for MinGW:
	rem We need printf exponents to have 2 digits like the C standard says,
	rem or tests will fail. Set environment variable for mingw-w64:
	rem https://github.com/Alexpux/mingw-w64/blob/master/mingw-w64-crt/stdio/mingw_pformat.c#L223
	set PRINTF_EXPONENT_DIGITS=2

	set HOME=.
	set MSYSTEM=MINGW64
	set "PATH=C:/msys64/usr/bin;C:/msys64/mingw64/bin;%PATH%"

	set CC=%MINGWPREFIX%-gcc.exe
	set CXX=%MINGWPREFIX%-g++.exe
	set FC=%MINGWPREFIX%-gfortran.exe
	sh -lc "pacman -S --noconfirm --needed mingw-w64-x86_64-openblas mingw-w64-x86_64-arpack mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-ncurses mingw-w64-x86_64-readline mingw64/mingw-w64-x86_64-netcdf diffutils"

	rem path from here: https://www.appveyor.com/docs/build-environment/#mingw-msys-cygwin
	set MINGWDIR=C:\mingw-w64\x86_64-6.3.0-posix-seh-rt_v5-rev1
)

if %BUILD_TYPE% equ cmake-mingw (
			
	rem add mingw32-make to the PATH
	set "PATH=%PATH%;%MINGWDIR%\bin"
	
) 

