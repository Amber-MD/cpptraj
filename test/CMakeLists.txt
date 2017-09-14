
if(INSTALL_TESTS)
	if(INSIDE_AMBER)
		
		# install to full path inside AmberTools directory
		install(DIRECTORY . USE_SOURCE_PERMISSIONS DESTINATION AmberTools/src/cpptraj/test/ COMPONENT Tests)
	else()
		# install to top-level dir
		install(DIRECTORY . USE_SOURCE_PERMISSIONS DESTINATION test/ COMPONENT Tests)
	endif()
endif()