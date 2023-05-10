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
    interfaces/burn_type.H
    interfaces/eos.H
    interfaces/eos_data.H
    interfaces/eos_type.H
    interfaces/eos_override.H
    interfaces/network.H
    interfaces/rhs_type.H
    networks/rate_type.H
    networks/rhs.H
    unit_test/unit_test.H)

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


#set a cache variable that controls whether
#we want to build the unit test or not
set(BUILD_UNIT_TEST false CACHE BOOL "Do you want to build the unit test? (true/false)")
message("Building unit test -- ${BUILD_UNIT_TEST}")

#setting sourcefiles and directories needed to make the test here
#so that they are accessible to codes using
#Microphysics as a submodule
#include parent_scope for sourcefiles and directories
#when built as a submodule in an external code
#otherwise these are not visible to the external code
#this will generate a warning but it will build successfully
#do not need unit_test paramfiles when unit_test is not built


#use this function to build Microphysics using different EOS for different test problems in external codes
function(setup_target_for_microphysics_compilation network_name output_dir)

  set(networkparamfile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/_parameters")
  set(VODEparamfile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/integration/VODE/_parameters")
  set(integrationparamfile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/integration/_parameters")
  set(unittestparamfile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/unit_test/_parameters")

  set(networkpropfile "${output_dir}/network_properties.H")
  set(networkheadertemplatefile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/general_null/network_header.template")

  if (${network_name} STREQUAL "gamma_law")
    set(EOSparamfile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/EOS/gamma_law/_parameters")

    set(networkfile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/general_null/gammalaw.net")
    set(networkdir "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/general_null/")
    set(networkpropertiesfile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/general_null/NETWORK_PROPERTIES")

    set (gamma_law_dirs ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/util ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/util/gcem/include 
                             ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/interfaces
                             ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/EOS ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/EOS/gamma_law 
                             ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/constants PARENT_SCOPE)

    execute_process(COMMAND python3 "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/util/build_scripts/write_probin.py" --pa "${EOSparamfile} ${networkparamfile} 
                    ${VODEparamfile} ${integrationparamfile} " WORKING_DIRECTORY ${output_dir}/)

    set(gamma_law_sources ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/interfaces/eos_data.cpp
                          ${output_dir}/extern_parameters.cpp PARENT_SCOPE)

  elseif (${network_name} STREQUAL "primordial_chem")
    #need these to write extern_parameters.H
    set(paramfile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/unit_test/burn_cell_primordial_chem/_parameters")
    set(EOSparamfile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/EOS/primordial_chem/_parameters")
    set(networkpcparamfile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/primordial_chem/_parameters")

    #similarly, we want network_properties.H
    set(networkfile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/primordial_chem/pynucastro.net")
    set(networkdir "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/primordial_chem/")
    set(networkpropertiesfile "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/primordial_chem/NETWORK_PROPERTIES")

    set (primordial_chem_dirs ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/util ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/util/gcem/include 
                             ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/integration/VODE ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/integration/utils 
                             ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/integration ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/interfaces
                             ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/EOS ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/EOS/primordial_chem 
                             ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/primordial_chem ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks 
                             ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/constants 
                             ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/general_null PARENT_SCOPE)

    #we need to have extern_parameters.cpp be available at configure time
    #the script write_probin.py writes this .cpp file so we call it here
    #note, execute_process only works on 'cmake' and not 'make'
    #so, if any of the _parameter files are changed, one needs to re-run 'cmake'
    #to generate updated header files

    if(BUILD_UNIT_TEST)
      execute_process(COMMAND python3 "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/util/build_scripts/write_probin.py" --pa "${paramfile} ${EOSparamfile} 
                      ${networkpcparamfile} ${networkparamfile} ${VODEparamfile} ${integrationparamfile} ${unittestparamfile}" WORKING_DIRECTORY ${output_dir}/)
    else()
      #do not need paramfile and unittestparamfile
      execute_process(COMMAND python3 "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/util/build_scripts/write_probin.py" --pa "${EOSparamfile} ${networkpcparamfile} 
                      ${networkparamfile} ${VODEparamfile} ${integrationparamfile} " WORKING_DIRECTORY ${output_dir}/)
    endif()

    set(primordial_chem_sources ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/interfaces/eos_data.cpp
                               ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/interfaces/network_initialization.cpp
                               ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/EOS/primordial_chem/actual_eos_data.cpp
                               ${output_dir}/extern_parameters.cpp PARENT_SCOPE)

    execute_process(COMMAND python3 "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/general_null/write_network.py" --header_template "${networkheadertemplatefile}" --header_output "${networkpropfile}" -s "${networkfile}" --other_properties "${networkpropertiesfile}" WORKING_DIRECTORY ${output_dir}/)

    #below to update headers
    execute_process(COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${PYTHONPATH}:${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/general_null" python3 "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/update_headers.py" --net "${networkdir}" --microphysics_path "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/" --odir "${output_dir}" WORKING_DIRECTORY ${output_dir}/)
    #below for NAUX
    execute_process(COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${PYTHONPATH}:${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/general_null" python3 "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/networks/get_naux.py" --net "${networkdir}" --microphysics_path "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/" WORKING_DIRECTORY ${output_dir}/)

  endif()


endfunction()



#set a cache variable that controls whether
#we want to fetch and link amrex or not
set(BUILD_AMREX false CACHE BOOL "Do you want to build and link amrex? (true/false)")
message("Building and linking amrex -- ${BUILD_AMREX}")

if(BUILD_AMREX)
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
   message("Built and fetched amrex in Microphysics")
else()
   message("Not building and fetching amrex in Microphysics")
endif()


if(BUILD_UNIT_TEST)
  setup_target_for_microphysics_compilation("primordial_chem" ${CMAKE_BINARY_DIR})

  include_directories(${primordial_chem_dirs})

  #adding unit_test as subdirectories
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/unit_test)
  message("Building primordial chem unit test")

else()

  message("Not building primordial chem unit test")
endif()



