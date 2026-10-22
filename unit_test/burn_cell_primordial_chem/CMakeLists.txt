#all .cpp files that need to be compiled for the test go here
add_executable(test_burn_cell_primordial_chem main.cpp ${primordial_chem_sources})

#test needs to find burn_cell.H, unit_test.H, and extern_parameters.H
target_include_directories(test_burn_cell_primordial_chem PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} "includes/burn_cell.H"
                           ${CMAKE_CURRENT_SOURCE_DIR}/../ "includes/unit_test.H" ${CMAKE_BINARY_DIR} "includes/extern_parameters.H")

#below line will add NAUX_NET and STRANG compiler flags
#these are needed for preprocessor variables
target_compile_definitions(test_burn_cell_primordial_chem PUBLIC NAUX_NET STRANG)


if(AMReX_GPU_BACKEND MATCHES "CUDA")
    setup_target_for_cuda_compilation(test_burn_cell_primordial_chem)
endif(AMReX_GPU_BACKEND MATCHES "CUDA")

add_test(NAME burn_cell_primordial_chem COMMAND test_burn_cell_primordial_chem inputs_primordial_chem WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

if(AMReX_GPU_BACKEND MATCHES "CUDA" AND MICROPHYSICS_INTEGRATOR STREQUAL "Rosenbrock")
    enable_language(CUDA)
    add_executable(test_primordial_rodas5p_hopper
                   hopper_compare.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/../../integration/Rosenbrock/primordial_rodas5p_hopper.cu
                   ${primordial_chem_sources})
    target_include_directories(test_primordial_rodas5p_hopper PUBLIC
                               ${CMAKE_CURRENT_SOURCE_DIR}
                               ${CMAKE_CURRENT_SOURCE_DIR}/../
                               ${CMAKE_BINARY_DIR})
    target_compile_definitions(test_primordial_rodas5p_hopper PUBLIC NAUX_NET STRANG)
    setup_target_for_cuda_compilation(test_primordial_rodas5p_hopper)
    add_test(NAME primordial_rodas5p_hopper
             COMMAND test_primordial_rodas5p_hopper inputs_primordial_chem integrator.rosenbrock_tableau=0
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
