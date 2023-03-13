add_subdirectory(burn_cell_primordial_chem)

add_dependencies(burn_cell_primordial_chem ${CMAKE_CURRENT_BINARY_DIR}/extern_parameters.h)

set(depfile "${CMAKE_CURRENT_SOURCE_DIR}/NETWORK_DIR??" ABSOLUTE)

add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/extern_parameters.h"
    COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/util/build_scripts/write_probin.py"
    DEPFILE ${depfile} --argument
    DEPENDS some_target)

add_custom_target(generate_header ALL
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/extern_parameters.h")
