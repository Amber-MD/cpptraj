@echo on

rem must have forward slashes
set SRCDIR=C:/projects/cpptraj

if %BUILD_TYPE% equ cmake-vs (
	
	mkdir build
	cd build
	
	cmake .. "-GVisual Studio 14 2015 Win64" -DCMAKE_PREFIX_PATH=%PREBUILTS_DIR% -DPRINT_PACKAGING_REPORT=TRUE -DARCHIVE_FORMAT=ZIP -DINSTALL_HEADERS=TRUE -DCOMPILER=MSVC -DCMAKE_INSTALL_PREFIX=%SRCDIR%
	msbuild /m cpptraj.sln
	msbuild /m INSTALL.vcxproj
	msbuild /m PACKAGE.vcxproj
)

if %BUILD_TYPE% equ cmake-mingw (

	mkdir build
		
	cd build
	
	rem make sure to pick up netcdf in /usr/local
	cmake .. "-GMinGW Makefiles" -DCMAKE_LIBRARY_PATH=C:/msys64/mingw64/lib;C:/msys64/usr/local/lib -DCMAKE_INCLUDE_PATH=C:/msys64/mingw64/include;C:/msys64/usr/local/include "-DCMAKE_SH=" -DPACKAGE_TYPE=ARCHIVE -DPRINT_PACKAGING_REPORT=TRUE -DARCHIVE_FORMAT=ZIP -DINSTALL_HEADERS=TRUE -DCOMPILER=MANUAL -DCMAKE_INSTALL_PREFIX=%SRCDIR% || exit /b
	mingw32-make -j2 install || exit /b
	mingw32-make -j2 package || exit /b
	cd ..
)

if %BUILD_TYPE% equ configure-mingw (
	sh -lc "./configure --with-netcdf=/usr/local/ --with-blas=/mingw64/ -openblas --with-bzlib=/mingw64/ --with-zlib=/mingw64 --with-arpack=/mingw64 --with-readline=/mingw64/ -shared -windows gnu" || exit /b
	make libcpptraj -j2 || exit /b
	make install -j2 || exit /b
)

