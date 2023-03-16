add_executable(test_burn_cell_primordial_chem main.cpp)
include_directories(${CMAKE_SOURCE_DIR}/unit_test/burn_cell_primordial_chem "includes/burn_cell.H" ${CMAKE_SOURCE_DIR}/unit_test "includes/unit_test.H")


if(AMReX_GPU_BACKEND MATCHES "CUDA")
    setup_target_for_cuda_compilation(test_burn_cell_primordial_chem)
endif(AMReX_GPU_BACKEND MATCHES "CUDA")

add_test(NAME burn_cell_primordial_chem COMMAND test_burn_cell_primordial_chem inputs_primordial_chem WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/unit_test/burn_cell_primordial_chem/)
