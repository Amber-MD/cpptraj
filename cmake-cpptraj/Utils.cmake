# this macro removes the lib prefix on each of the library targets provided so that this doesn't happen
macro(remove_prefix) #LIBRARIES
        set_target_properties(${ARGN} PROPERTIES PREFIX "")
        set_target_properties(${ARGN} PROPERTIES IMPORT_PREFIX "")
endmacro(remove_prefix)

#make the provided object library position independent if shared libraries are turned on
function(make_pic_if_needed OBJECT_LIBRARY)
        set_property(TARGET ${OBJECT_LIBRARY} PROPERTY POSITION_INDEPENDENT_CODE ${SHARED})
endfunction(make_pic_if_needed)

#Checks that the cache variable VARIABLE is set to one of VALID_VALUES and prints an error if is not.
#Also creates a pull-down menu for the variable in the GUI containing these choices
macro(validate_configuration_enum VARIABLE) #2nd argument: VALID_VALUES...

        list_contains(VALID ${${VARIABLE}} ${ARGN})

        if(NOT VALID)
                list_to_space_separated(VALID_VALUES_STRING ${ARGN})

                message(FATAL_ERROR "${${VARIABLE}} is not a valid value for ${VARIABLE} -- must be one of: ${VALID_VALUES_STRING}")
        endif()

          set_property(CACHE ${VARIABLE} PROPERTY STRINGS ${ARGN})
endmacro(validate_configuration_enum)

#converts a list into a string with each of its elements seperated by a space
macro(list_to_space_separated OUTPUT_VAR)# 2nd arg: LIST...
        string(REPLACE ";" " " ${OUTPUT_VAR} "${ARGN}")
endmacro(list_to_space_separated)

