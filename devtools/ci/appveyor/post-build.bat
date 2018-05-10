
echo "hello"

if %BUILD_TYPE% equ "configure-mingw" (
	
	7z a cpptraj-%APPVEYOR_BUILD_ID%.zip bin/ambpdb.exe bin/cpptraj.exe lib/libcpptraj.dll.a src/*.h
	
) else (
	rem move and rename already created package file
	echo F | xcopy /Y /F build\cpptraj.zip .\cpptraj-%APPVEYOR_BUILD_ID%.zip
)
