project(Microphysics
    VERSION 1.0.0
    DESCRIPTION "some description of the project; built with CMake"
    LANGUAGES CXX)

#----------------------------------------------------------------------------------------------------------------------
# general settings and options
#----------------------------------------------------------------------------------------------------------------------

add_library(Microphysics) # initialized below
add_library(Microphysics::Microphysics ALIAS Microphysics)

# Search for your dependencies here

#----------------------------------------------------------------------------------------------------------------------
# mylib sources
#----------------------------------------------------------------------------------------------------------------------

set(sources
    EOS/eos_composition.h
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

#fetching amrex
FetchContent_Declare(
  amrex
  GIT_REPOSITORY https://github.com/AMReX-Codes/amrex.git
  GIT_TAG        v23.03
)

# CMake 3.14+
FetchContent_MakeAvailable(amrex)


#adding unit_test as subdirectories
add_subdirectory(${Microphysics_SOURCE_DIR}/unit_test)


#----------------------------------------------------------------------------------------------------------------------
# mylib target
#----------------------------------------------------------------------------------------------------------------------

