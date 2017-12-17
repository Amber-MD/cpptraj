
echo "hello"

if %USE_CMAKE equ 1 (
	
	rem move and rename already created package file
	xcopy /Y /F build\cpptraj.zip .\cpptraj-%APPVEYOR_BUILD_ID%.zip
	
	rem move executables to bin folder
	mkdir bin
	xcopy /Y /F /I build\src\*.exe bin\
	
) else (
	7z a cpptraj-%APPVEYOR_BUILD_ID%.zip bin/ambpdb.exe bin/cpptraj.exe lib/libcpptraj.dll.a src/*.h
)

find C:/msys64/usr/local