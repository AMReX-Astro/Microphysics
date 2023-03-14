add_subdirectory(burn_cell_primordial_chem)

#add_dependencies(test_burn_cell_primordial_chem ${CMAKE_CURRENT_SOURCE_DIR}/extern_parameters.H)

set(depfile "${MICROPHYSICS_AUTO_SOURCE_DIR}/extern_parameters.H")

#add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/extern_parameters.H"
#    COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/util/build_scripts/write_probin.py"
#    DEPFILE ${depfile} --argument
#    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/extern_parameters.H")

add_custom_command(OUTPUT ${depfile}
    COMMAND python3 "${CMAKE_CURRENT_SOURCE_DIR}/../util/build_scripts/write_probin.py"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/../util/build_scripts/write_probin.py"
    COMMENT "writing extern_parameters.H ${PYTHON_EXECUTABLE}")

add_custom_target(generate_extern_parameters
    DEPENDS ${depfile})

add_dependencies(test_burn_cell_primordial_chem generate_extern_parameters)

get_cmake_property(_variableNames VARIABLES)
list (SORT _variableNames)
foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
endforeach()
