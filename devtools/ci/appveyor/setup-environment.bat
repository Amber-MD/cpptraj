@echo on
set HOME=.
set MSYSTEM=MINGW64
set "PATH=C:/msys64/usr/bin;C:/msys64/mingw64/bin;%PATH%"

set MINGWPREFIX=x86_64-w64-mingw32
set CC=%MINGWPREFIX%-gcc.exe
set CXX=%MINGWPREFIX%-g++.exe
set FC=%MINGWPREFIX%-gfortran.exe
sh -lc "pacman -S --noconfirm --needed mingw-w64-x86_64-openblas mingw-w64-x86_64-arpack mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-ncurses mingw-w64-x86_64-readline diffutils"
sh -lc "^
 if [ ! -f /usr/local/lib/libnetcdf.a ]; then^
   curl -fsS -o netcdf-4.3.3.zip ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.3.zip &&^
   7z x netcdf-4.3.3.zip -o$APPVEYOR_BUILD_FOLDER > /dev/null &&^
   cd $APPVEYOR_BUILD_FOLDER/netcdf-4.3.3 &&^
   exec 0</dev/null &&^
   ./configure --enable-static --disable-netcdf-4 --prefix=/usr/local/ --disable-dap &&^
   make -r install;^
 else^
   echo 'Have Cached NetCDF';^
 fi"

if %USE_CMAKE equ 1 (
		
	rem path from here: https://www.appveyor.com/docs/build-environment/#mingw-msys-cygwin
	set MINGWDIR=C:\mingw-w64\x86_64-6.3.0-posix-seh-rt_v5-rev1
	
	rem add mingw32-make and makensis to the PATH
	set "PATH=%PATH%;%MINGWDIR%\bin;C:\Program Files (x86)\NSIS"
	
) 
