cmake_minimum_required( VERSION 3.11 FATAL_ERROR )
project( demotrack_cuda_cxx LANGUAGES C CXX )

option( USE_BLOCKSIZE_CALCULATOR "Use blocksize calculator" OFF )
set( DEMOTRACK_DEFAULT_BLOCK_SIZE "512" CACHE STRING "Default grid block size" )

if( USE_BLOCKSIZE_CALCULATOR )
    set( DEMOTRACK_USE_BLOCKSIZE_CALCULATOR "1" )
else()
    set( DEMOTRACK_USE_BLOCKSIZE_CALCULATOR "0" )
endif()

configure_file( config.h.template ${CMAKE_SOURCE_DIR}/config.h @ONLY )

# -------------------------------------------------------------------------------
# Enable CUDA support

# Set environment variable CUDA_BIN_PATH to the directory containing nvcc or
# HIP_ROOT_DIR variable to hint the CUDA install location

if(NOT DEFINED HIP_PATH)
    if(NOT DEFINED ENV{HIP_PATH})
        set(HIP_PATH "/opt/rocm/hip" CACHE PATH "Path to which HIP has been installed")
    else()
        set(HIP_PATH $ENV{HIP_PATH} CACHE PATH "Path to which HIP has been installed")
    endif()
endif()

set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" ${CMAKE_MODULE_PATH})
find_package( HIP REQUIRED )

# ******************************************************************************
# demo applications

set( CLANG_MCPU_ARCH "gfx900" ) #--amdgpu-target=${CLANG_MCPU_ARCH}

set( CMAKE_POSITION_INDEPENDENT_CODE ON )
set( CMAKE_CXX_STANDARD 11 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CUDA_STANDARD 11 )
set( CMAKE_CUDA_STANDARD_REQUIRED ON )

include_directories( ${CMAKE_SOURCE_DIR} )

# ------------------------------------------------------------------------------
# demo01: particles stored in global memory, lattice stored in global memory

set( COMMON_HEADERS definitions.h beam_elements.h beamfields.h config.h
                    lattice.h particle.h )

HIP_ADD_EXECUTABLE( demo01_sc0 demo01.cpp ${COMMON_HEADERS} )
target_compile_options( demo01_sc0 PRIVATE -fno-fast-math )
target_link_libraries( demo01_sc0 PRIVATE ${DEMOTRACK_LINK_LIBRARIES} )

HIP_ADD_EXECUTABLE( demo02_sc0 demo02.cpp ${COMMON_HEADERS} )
target_compile_options( demo02_sc0 PRIVATE -fno-fast-math )
target_link_libraries( demo02_sc0 PRIVATE ${DEMOTRACK_LINK_LIBRARIES} )

# ------------------------------------------------------------------------------
# demo02: particles stored in global memory, lattice stored in global memory
#         BUT: tracking with thread-local / private copy of the particles,
#              i.e. avoiding the frequent global memory accesses

set( MY_HIPCC_FLAGS "${MY_HIPCC_FLAGS};-DDEMOTRACK_ENABLE_BEAMFIELDS=1" )

HIP_ADD_EXECUTABLE( demo01_sc1 demo01.cpp ${COMMON_HEADERS} )
target_link_libraries( demo01_sc1 PRIVATE ${DEMOTRACK_LINK_LIBRARIES} )
target_compile_options( demo01_sc1 PRIVATE -fno-fast-math )
target_compile_definitions( demo01_sc1 PRIVATE -DDEMOTRACK_ENABLE_BEAMFIELDS=1 )

HIP_ADD_EXECUTABLE( demo02_sc1 demo02.cpp ${COMMON_HEADERS}
                    HIPCC_OPTIONS ${MY_HIPCC_FLAGS} )
target_link_libraries( demo02_sc1 PRIVATE ${DEMOTRACK_LINK_LIBRARIES} )
target_compile_options( demo02_sc1 PRIVATE -fno-fast-math )
target_compile_definitions( demo02_sc1 PRIVATE -DDEMOTRACK_ENABLE_BEAMFIELDS=1 )
