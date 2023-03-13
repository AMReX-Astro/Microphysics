add_subdirectory(burn_cell_primordial_chem)

#add_dependencies(burn_cell_primordial_chem ${CMAKE_CURRENT_BINARY_DIR}/extern_parameters.h)

set(depfile "${CMAKE_CURRENT_SOURCE_DIR}/extern_parameters.H" ABSOLUTE)

add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/extern_parameters.H"
    COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/util/build_scripts/write_probin.py"
    DEPFILE ${depfile} --argument
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/extern_parameters.H")

#add_custom_target(generate_header ALL
#    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/extern_parameters.h")
