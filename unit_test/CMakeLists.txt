#include these directories based on the output of the traditional make
include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR}/util ${CMAKE_SOURCE_DIR}/util/gcem/include ${CMAKE_SOURCE_DIR}/integration/VODE ${CMAKE_SOURCE_DIR}/integration/utils ${CMAKE_SOURCE_DIR}/integration
                    ${CMAKE_SOURCE_DIR}/EOS ${CMAKE_SOURCE_DIR}/EOS/primordial_chem ${CMAKE_SOURCE_DIR}/networks/primordial_chem ${CMAKE_SOURCE_DIR}/networks ${CMAKE_SOURCE_DIR}/interfaces
                    ${CMAKE_SOURCE_DIR}/constants ${CMAKE_SOURCE_DIR}/networks/general_null)


#add primordial_chem test
add_subdirectory(${CMAKE_SOURCE_DIR}/unit_test/burn_cell_primordial_chem)

set(networkpropfile "${CMAKE_CURRENT_BINARY_DIR}/network_properties.H")
set(networkfile "${CMAKE_SOURCE_DIR}/networks/primordial_chem/pynucastro.net")
set(networkdir "${CMAKE_SOURCE_DIR}/networks/primordial_chem/")
set(networkpropertiesfile "${CMAKE_SOURCE_DIR}/networks/primordial_chem/NETWORK_PROPERTIES")
set(networkheadertemplatefile "${CMAKE_SOURCE_DIR}/networks/general_null/network_header.template")


add_custom_command(OUTPUT ${networkpropfile}
    COMMAND python3 "${CMAKE_SOURCE_DIR}/networks/general_null/write_network.py" --header_template "${networkheadertemplatefile}" --header_output "${networkpropfile}" -s "${networkfile}" --other_properties "${networkpropertiesfile}"
    COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${PYTHONPATH}:${CMAKE_SOURCE_DIR}/networks/general_null" python3 "${CMAKE_SOURCE_DIR}/networks/get_naux.py" --net "${networkdir}" --microphysics_path "${CMAKE_SOURCE_DIR}/"
    COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${PYTHONPATH}:${CMAKE_SOURCE_DIR}/networks/general_null" python3 "${CMAKE_SOURCE_DIR}/networks/update_headers.py" --net "${networkdir}" --microphysics_path "${CMAKE_SOURCE_DIR}/" --odir "${CMAKE_CURRENT_BINARY_DIR}"
    DEPENDS "${CMAKE_SOURCE_DIR}/networks/get_naux.py"
    DEPENDS "${CMAKE_SOURCE_DIR}/networks/update_headers.py"
    DEPENDS "${CMAKE_SOURCE_DIR}/networks/general_null/write_network.py"
    DEPENDS "${CMAKE_SOURCE_DIR}/networks/general_null/network_param_file.py"
    DEPENDS "${CMAKE_SOURCE_DIR}/networks/general_null/network_header.template"
    DEPENDS "${networkfile}"
    COMMENT "writing network_properties.H"
    )

add_custom_target(generate_network_properties
    DEPENDS ${networkpropfile}
    )

#the test depends on this custom target
add_dependencies(test_burn_cell_primordial_chem generate_network_properties)

#uncomment below if you want to see the values of all CMAKE variables
#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
