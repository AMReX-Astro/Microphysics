include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR}/util ${CMAKE_SOURCE_DIR}/util/gcem/include ${CMAKE_SOURCE_DIR}/integration/VODE ${CMAKE_SOURCE_DIR}/integration/utils ${CMAKE_SOURCE_DIR}/integration
                    ${CMAKE_SOURCE_DIR}/EOS ${CMAKE_SOURCE_DIR}/EOS/primordial_chem ${CMAKE_SOURCE_DIR}/networks/primordial_chem ${CMAKE_SOURCE_DIR}/networks ${CMAKE_SOURCE_DIR}/interfaces
                    ${CMAKE_SOURCE_DIR}/constants ${CMAKE_SOURCE_DIR}/networks/general_null)

add_subdirectory(burn_cell_primordial_chem)

set(externparamfile "${CMAKE_CURRENT_BINARY_DIR}/extern_parameters.H")

set(paramfile "${CMAKE_SOURCE_DIR}/unit_test/burn_cell_primordial_chem/_parameters")
set(EOSparamfile "${CMAKE_SOURCE_DIR}/EOS/primordial_chem/_parameters")
set(networkpcparamfile "${CMAKE_SOURCE_DIR}/networks/primordial_chem/_parameters")
set(networkparamfile "${CMAKE_SOURCE_DIR}/networks/_parameters")
set(VODEparamfile "${CMAKE_SOURCE_DIR}/integration/VODE/_parameters")
set(integrationparamfile "${CMAKE_SOURCE_DIR}/integration/_parameters")
set(unittestparamfile "${CMAKE_SOURCE_DIR}/unit_test/_parameters")


set(networkpropfile "${CMAKE_CURRENT_BINARY_DIR}/network_properties.H")
set(networkfile "${CMAKE_SOURCE_DIR}/networks/primordial_chem/pynucastro.net")
set(networkdir "${CMAKE_SOURCE_DIR}/networks/primordial_chem/")
set(networkpropertiesfile "${CMAKE_SOURCE_DIR}/networks/primordial_chem/NETWORK_PROPERTIES")
set(networkheadertemplatefile "${CMAKE_SOURCE_DIR}/networks/general_null/network_header.template")


add_custom_command(OUTPUT ${externparamfile}
    COMMAND python3 "${CMAKE_SOURCE_DIR}/util/build_scripts/write_probin.py" --pa "${paramfile} ${EOSparamfile} ${networkpcparamfile} ${networkparamfile} ${VODEparamfile} ${integrationparamfile} ${unittestparamfile}"
    DEPENDS "${CMAKE_SOURCE_DIR}/util/build_scripts/write_probin.py"
    DEPENDS "${CMAKE_SOURCE_DIR}/util/build_scripts/runtime_parameters.py"
    DEPENDS "${paramfile}"
    COMMENT "writing extern_parameters.H")

add_custom_target(generate_extern_parameters
    DEPENDS ${externparamfile}
    )

add_dependencies(test_burn_cell_primordial_chem generate_extern_parameters)


#PYTHONPATH=../../networks/general_null ../../networks/update_headers.py \
#           --microphysics_path ../.. \
#           --net primordial_chem \
#           --odir tmp_build_dir/microphysics_sources/1d.gnu.DEBUG.EXE \
#           --defines "-DCXX_REACTIONS -DMICROPHYSICS_DEBUG -DBL_SPACEDIM=1 -DAMREX_SPACEDIM=1 -DBL_FORT_USE_UNDERSCORE -DAMREX_FORT_USE_UNDERSCORE -DBL_Linux -DAMREX_Linux -DBL_USE_ASSERTION -DAMREX_USE_ASSERTION -DCRSEGRNDOMP -DREACTIONS -DSTRANG -DSCREEN_METHOD=0 -DNETWORK_HAS_CXX_IMPLEMENTATION -DALLOW_JACOBIAN_CACHING "-DNAUX_NET=0""
#write_network.py: working on network file ../../networks/primordial_chem/pynucastro.net ...
#writing tmp_build_dir/microphysics_sources/1d.gnu.DEBUG.EXE/network_properties.H

#COMMAND python3 "${CMAKE_SOURCE_DIR}/networks/general_null/network_param_file.py"

#set(ENV{PYTHONPATH} "${CMAKE_SOURCE_DIR}/networks/general_null")
#set(PYTHON_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/networks/general_null")

#Python_add_library("${CMAKE_SOURCE_DIR}/networks/general_null")

#export ${PYTHONPATH}="${PYTHONPATH}:${CMAKE_SOURCE_DIR}/networks/general_null"


#COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_SOURCE_DIR}/networks/general_null"
#COMMAND set "PYTHONPATH=${CMAKE_SOURCE_DIR}/networks/general_null;%PYTHONPATH%"

add_custom_command(OUTPUT ${networkpropfile}
    COMMAND python3 "${CMAKE_SOURCE_DIR}/networks/general_null/write_network.py" --header_template "${networkheadertemplatefile}" --header_output "${networkpropfile}" -s "${networkfile}" --other_properties "${networkpropertiesfile}"
    #COMMAND python3 "${CMAKE_SOURCE_DIR}/networks/general_null/write_network.py" --header_template "${networkheadertemplatefile}" --header_output "${networkpropfile}" -s "${networkfile}" --other_properties "${networkpropertiesfile}"
    COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${PYTHONPATH}:${CMAKE_SOURCE_DIR}/networks/general_null" python3 "${CMAKE_SOURCE_DIR}/networks/update_headers.py" --net "${networkdir}" --microphysics_path "${CMAKE_SOURCE_DIR}/" --odir "${CMAKE_CURRENT_BINARY_DIR}"
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

add_dependencies(test_burn_cell_primordial_chem generate_network_properties)


get_cmake_property(_variableNames VARIABLES)
list (SORT _variableNames)
foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
endforeach()
