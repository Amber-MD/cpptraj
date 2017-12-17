@echo on

rem must have forward slashes
set SRCDIR=C:/projects/cpptraj

if %USE_CMAKE% equ 1 (

	mkdir build
		
	cd build
	
	rem make sure to pick up netcdf in /usr/local
	cmake .. "-GMinGW Makefiles" -DCMAKE_LIBRARY_PATH=C:/msys64/mingw64/lib;C:/msys64/usr/local/lib -DFORCE_DISABLE_LIBS=arpack -DCMAKE_INCLUDE_PATH=C:/msys64/mingw64/include;C:/msys64/usr/local/include "-DCMAKE_SH=" -DPACKAGE_TYPE=ARCHIVE -DPRINT_PACKAGING_REPORT=TRUE -DARCHIVE_FORMAT=ZIP -DINSTALL_HEADERS=TRUE -DCOMPILER=manual -DCMAKE_INSTALL_PREFIX=%SRCDIR% || exit /b
	mingw32-make -j2 install || exit /b
	mingw32-make -j2 package || exit /b
	cd ..
) else (
	"sh configure --with-netcdf=/usr/local/ --with-blas=/mingw64/ -openblas --with-bzlib=/mingw64/ --with-zlib=/mingw64 --with-arpack=/mingw64 --with-readline=/mingw64/ -shared -windows gnu" || exit /b
	make libcpptraj -j2 || exit /b
	make install -j2 || exit /b
)

