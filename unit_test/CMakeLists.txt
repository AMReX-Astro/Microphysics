#include these directories based on the output of the traditional make
include_directories(${primordial_chem_dirs})

#add primordial_chem test
add_subdirectory(burn_cell_primordial_chem)

set(networkpropfile "${CMAKE_BINARY_DIR}/network_properties.H")
set(networkfile "${CMAKE_CURRENT_SOURCE_DIR}/../networks/primordial_chem/pynucastro.net")
set(networkdir "${CMAKE_CURRENT_SOURCE_DIR}/../networks/primordial_chem/")
set(networkpropertiesfile "${CMAKE_CURRENT_SOURCE_DIR}/../networks/primordial_chem/NETWORK_PROPERTIES")
set(networkheadertemplatefile "${CMAKE_CURRENT_SOURCE_DIR}/../networks/general_null/network_header.template")


add_custom_command(OUTPUT ${networkpropfile}
    COMMAND python3 "${CMAKE_CURRENT_SOURCE_DIR}/../networks/general_null/write_network.py" --header_template "${networkheadertemplatefile}" --header_output "${networkpropfile}" -s "${networkfile}" --other_properties "${networkpropertiesfile}"
    COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${PYTHONPATH}:${CMAKE_CURRENT_SOURCE_DIR}/../networks/general_null" python3 "${CMAKE_CURRENT_SOURCE_DIR}/../networks/get_naux.py" --net "${networkdir}" --microphysics_path "${CMAKE_CURRENT_SOURCE_DIR}/../"
    COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${PYTHONPATH}:${CMAKE_CURRENT_SOURCE_DIR}/../networks/general_null" python3 "${CMAKE_CURRENT_SOURCE_DIR}/../networks/update_headers.py" --net "${networkdir}" --microphysics_path "${CMAKE_CURRENT_SOURCE_DIR}/../" --odir "${CMAKE_BINARY_DIR}"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/../networks/get_naux.py"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/../networks/update_headers.py"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/../networks/general_null/write_network.py"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/../networks/general_null/network_param_file.py"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/../networks/general_null/network_header.template"
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
