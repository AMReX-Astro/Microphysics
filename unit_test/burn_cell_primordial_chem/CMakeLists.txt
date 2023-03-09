add_executable(test_burn_cell_primordial_chem main.cpp)

if(AMReX_GPU_BACKEND MATCHES "CUDA")
    setup_target_for_cuda_compilation(test_burn_cell_primordial_chem)
endif(AMReX_GPU_BACKEND MATCHES "CUDA")

add_test(NAME burn_cell_primordial_chem COMMAND test_burn_cell_primordial_chem inputs_primordial_chem WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/unit_test/burn_cell_primordial_chem/)

add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/extern_parameters.h"
    COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/util/build_scripts/write_probin.py" --argument
    DEPENDS some_target)

add_custom_target(generate_header ALL
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/extern_parameters.h")

