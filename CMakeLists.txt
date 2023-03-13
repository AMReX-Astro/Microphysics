cmake_minimum_required(VERSION 3.16)

project(Microphysics
    VERSION 1.0.0
    DESCRIPTION "some description of the project; built with CMake"
    LANGUAGES CXX C)

#----------------------------------------------------------------------------------------------------------------------
# general settings and options
#----------------------------------------------------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_COMPILER "clang" CACHE PATH "")
set(CMAKE_CXX_COMPILER "clang++" CACHE PATH "")
include(CTest)

add_library(Microphysics) # initialized below
add_library(Microphysics::Microphysics ALIAS Microphysics)

# Search for your dependencies here

#----------------------------------------------------------------------------------------------------------------------
# mylib sources
#----------------------------------------------------------------------------------------------------------------------

set(sources
    EOS/eos_composition.H
    integration/integrator.H
    integration/integrator_data.H
    interfaces/ArrayUtilities.H
    interfaces/eos.H
    interfaces/eos_data.H
    interfaces/eos_type.H
    interfaces/eos_override.H
    interfaces/network.H
    interfaces/rhs_type.H
    networks/rate_type.H
    networks/rhs.H)
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${sources})

#----------------------------------------------------------------------------------------------------------------------
# mylib target
#----------------------------------------------------------------------------------------------------------------------

include(CMakePackageConfigHelpers)

target_sources(Microphysics PRIVATE ${sources})
#target_compile_definitions(Microphysics PUBLIC "$<$<NOT:$<BOOL:${BUILD_SHARED_LIBS}>>:MYLIB_STATIC_DEFINE>")

#target_include_directories(Microphysics
#    PUBLIC
#        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
#        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>")

set_target_properties(Microphysics PROPERTIES
    SOVERSION ${PROJECT_VERSION_MAJOR}
    VERSION ${PROJECT_VERSION}
    LINKER_LANGUAGE CXX)


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


#adding unit_test as subdirectories
add_subdirectory(${Microphysics_SOURCE_DIR}/unit_test)


