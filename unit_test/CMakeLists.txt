include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR}/util ${CMAKE_SOURCE_DIR}/util/gcem/include ${CMAKE_SOURCE_DIR}/integration/VODE ${CMAKE_SOURCE_DIR}/integration/utils ${CMAKE_SOURCE_DIR}/integration
                    ${CMAKE_SOURCE_DIR}/EOS ${CMAKE_SOURCE_DIR}/EOS/primordial_chem ${CMAKE_SOURCE_DIR}/networks/primordial_chem ${CMAKE_SOURCE_DIR}/networks ${CMAKE_SOURCE_DIR}/interfaces
                    ${CMAKE_SOURCE_DIR}/constants)

add_subdirectory(burn_cell_primordial_chem)

set(externparamfile "${CMAKE_CURRENT_BINARY_DIR}/extern_parameters.H")
set(paramfile "${CMAKE_SOURCE_DIR}/networks/primordial_chem/_parameters")
set(networkpropfile "${CMAKE_CURRENT_BINARY_DIR}/network_properties.H")
set(networkfile "${CMAKE_SOURCE_DIR}/networks/primordial_chem/pynucastro.net")
set(networkpropertiesfile "${CMAKE_SOURCE_DIR}/networks/primordial_chem/NETWORK_PROPERTIES")

add_custom_command(OUTPUT ${externparamfile}
    COMMAND python3 "${CMAKE_SOURCE_DIR}/util/build_scripts/write_probin.py" --pa "${paramfile}"
    DEPENDS "${CMAKE_SOURCE_DIR}/util/build_scripts/write_probin.py"
    DEPENDS "${CMAKE_SOURCE_DIR}/util/build_scripts/runtime_parameters.py"
    DEPENDS "${paramfile}"
    COMMENT "writing extern_parameters.H")

add_custom_target(generate_extern_parameters
    DEPENDS ${externparamfile}
    )

add_dependencies(test_burn_cell_primordial_chem generate_extern_parameters)


add_custom_command(OUTPUT ${networkpropfile}
    COMMAND python3 "${CMAKE_SOURCE_DIR}/networks/general_null/network_param_file.py"
    COMMAND python3 "${CMAKE_SOURCE_DIR}/networks/general_null/write_network.py" -s "${networkfile}" --other_properties "${networkpropertiesfile}"
    COMMAND python3 "${CMAKE_SOURCE_DIR}/networks/update_headers.py" --net "primordial_chem"
    DEPENDS "${CMAKE_SOURCE_DIR}/networks/update_headers.py"
    DEPENDS "${CMAKE_SOURCE_DIR}/networks/general_null/write_network.py"
    DEPENDS "${CMAKE_SOURCE_DIR}/networks/general_null/network_param_file.py"
    COMMENT "writing network_properties.H")

add_custom_target(generate_network_properties
    DEPENDS ${networkpropfile}
    )

add_dependencies(test_burn_cell_primordial_chem generate_network_properties)


get_cmake_property(_variableNames VARIABLES)
list (SORT _variableNames)
foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
endforeach()
