cmake_minimum_required(VERSION 3.20)

project(Microphysics
    VERSION 1.0.0
    DESCRIPTION "building primordial_chemistry network in Microphysics with CMake"
    LANGUAGES CXX C)

#----------------------------------------------------------------------------------------------------------------------
# general settings and options
#----------------------------------------------------------------------------------------------------------------------

find_package (Python COMPONENTS Interpreter Development)
find_package(Threads)

include(CTest)

add_library(Microphysics SHARED) # initialized below

#----------------------------------------------------------------------------------------------------------------------
# mylib sources
#----------------------------------------------------------------------------------------------------------------------

set(sources
    EOS/eos_composition.H
    integration/integrator.H
    integration/integrator_data.H
    integration/integrator_rhs_strang.H
    integration/integrator_type_strang.H
    interfaces/ArrayUtilities.H
    interfaces/burner.H
    interfaces/eos.H
    interfaces/eos_data.H
    interfaces/eos_type.H
    interfaces/eos_override.H
    interfaces/network.H
    interfaces/rhs_type.H
    networks/rate_type.H
    networks/rhs.H
    unit_test/unit_test.H)
#source_group(TREE "${CMAKE_SOURCE_DIR}" FILES ${sources})

#----------------------------------------------------------------------------------------------------------------------
# mylib target
#----------------------------------------------------------------------------------------------------------------------

include(CMakePackageConfigHelpers)

target_sources(Microphysics PUBLIC ${sources})
#below line will add NAUX_NET and STRANG compiler flags
#these are needed for preprocessor variables
target_compile_definitions(Microphysics PUBLIC NAUX_NET STRANG)

set_target_properties(Microphysics PROPERTIES
    SOVERSION ${PROJECT_VERSION_MAJOR}
    VERSION ${PROJECT_VERSION}
    LINKER_LANGUAGE CXX)

#we need to have extern_parameters.cpp be available at configure time
#the script write_probin.py writes this .cpp file so we call it here
set(paramfile "${CMAKE_SOURCE_DIR}/unit_test/burn_cell_primordial_chem/_parameters")
set(EOSparamfile "${CMAKE_SOURCE_DIR}/EOS/primordial_chem/_parameters")
set(networkpcparamfile "${CMAKE_SOURCE_DIR}/networks/primordial_chem/_parameters")
set(networkparamfile "${CMAKE_SOURCE_DIR}/networks/_parameters")
set(VODEparamfile "${CMAKE_SOURCE_DIR}/integration/VODE/_parameters")
set(integrationparamfile "${CMAKE_SOURCE_DIR}/integration/_parameters")
set(unittestparamfile "${CMAKE_SOURCE_DIR}/unit_test/_parameters")

execute_process(COMMAND python3 "${CMAKE_SOURCE_DIR}/util/build_scripts/write_probin.py" --pa "${paramfile} ${EOSparamfile} ${networkpcparamfile} ${networkparamfile} ${VODEparamfile} ${integrationparamfile} ${unittestparamfile}")

#fetching amrex
set (AMREX_TAG "development")
include(FetchContent)
FetchContent_Declare(
  amrex
  GIT_REPOSITORY https://github.com/AMReX-Codes/amrex
  GIT_TAG "${AMREX_TAG}"
)

# CMake 3.14+
FetchContent_MakeAvailable(amrex)

include_directories(${amrex_INCLUDE_DIRS_RET})
link_libraries(AMReX::amrex)

#adding unit_test as subdirectories
add_subdirectory(${CMAKE_SOURCE_DIR}/unit_test)


