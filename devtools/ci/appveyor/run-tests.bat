@echo on

set CPPTRAJ_TEST_OS=windows

rem debugging: print library dependencies of cpptraj
dumpbin /dependents bin/cpptraj.exe

rem currently tests can't be run under Visual Studio due to some strange library error
if %BUILD_TYPE% neq cmake-vs (
	sh -lc "cd test; make test.showerrors"
)