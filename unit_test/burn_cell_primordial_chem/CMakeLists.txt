#all .cpp files that need to be compiled for the test go here
add_executable(test_burn_cell_primordial_chem main.cpp ${primordial_chem_sources})

#test needs to find burn_cell.H, unit_test.H, and extern_parameters.H
target_include_directories(test_burn_cell_primordial_chem PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} "includes/burn_cell.H"
                           ${CMAKE_CURRENT_SOURCE_DIR}/../ "includes/unit_test.H" ${CMAKE_BINARY_DIR} "includes/extern_parameters.H")
#test needs to be linked to all of Microphysics
target_link_libraries(test_burn_cell_primordial_chem Microphysics)


if(AMReX_GPU_BACKEND MATCHES "CUDA")
    setup_target_for_cuda_compilation(test_burn_cell_primordial_chem)
endif(AMReX_GPU_BACKEND MATCHES "CUDA")

add_test(NAME burn_cell_primordial_chem COMMAND test_burn_cell_primordial_chem inputs_primordial_chem WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
