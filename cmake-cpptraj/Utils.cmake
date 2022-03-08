# this macro removes the lib prefix on each of the library targets provided so that this doesn't happen
macro(remove_prefix) #LIBRARIES
        set_target_properties(${ARGN} PROPERTIES PREFIX "")
        set_target_properties(${ARGN} PROPERTIES IMPORT_PREFIX "")
endmacro(remove_prefix)

#make the provided object library position independent if shared libraries are turned on
function(make_pic_if_needed OBJECT_LIBRARY)
        set_property(TARGET ${OBJECT_LIBRARY} PROPERTY POSITION_INDEPENDENT_CODE ${SHARED})
endfunction(make_pic_if_needed)

