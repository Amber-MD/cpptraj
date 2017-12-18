
echo "hello"

if %USE_CMAKE% equ 1 (
	
	rem move and rename already created package file
	echo F | xcopy /Y /F build\cpptraj.zip .\cpptraj-%APPVEYOR_BUILD_ID%.zip
	
) else (
	7z a cpptraj-%APPVEYOR_BUILD_ID%.zip bin/ambpdb.exe bin/cpptraj.exe lib/libcpptraj.dll.a src/*.h
)

find C:/msys64/usr/local