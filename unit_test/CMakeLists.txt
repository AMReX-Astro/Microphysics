
if(BUILD_UNIT_TEST_PC)
  #add primordial_chem test
  add_subdirectory(burn_cell_primordial_chem)
endif()

if(BUILD_UNIT_TEST_MC)
  #add metal_chem test
  add_subdirectory(burn_cell_metal_chem)
endif()

#uncomment below if you want to see the values of all CMAKE variables
#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
