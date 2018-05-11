@echo on

set CPPTRAJ_TEST_OS=windows

rem debugging: print library dependencies of cpptraj
dumpbin /dependents bin/cpptraj.exe

if %BUILD_TYPE% equ cmake-vs (
	rem make sure MSYS's sh inherits our full path
	rem (from https://github.com/cisco/ChezScheme/issues/19 )
	set MSYS2_PATH_TYPE=inherit
	
	set "PATH=%cd%\bin;%PATH%"
	
	cd C:/msys64/usr/bin
	
	cpptraj --version
	sh -lc "cd /c/projects/cpptraj/test; make test.showerrors"
	
) else (
	sh -lc "cd test; make test.showerrors" || exit /b
)