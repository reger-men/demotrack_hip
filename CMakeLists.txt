cmake_minimum_required( VERSION 3.11 FATAL_ERROR )
project( demotrack LANGUAGES C CXX )

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE "Release" CACHE STRING
       "Build type [Debug ReleaseRelWithDebInfo MinSizeRel]" FORCE )
endif()

if( EXISTS ${CMAKE_SOURCE_DIR}/Settings.cmake )
    message( STATUS "loading custom settings file ${CMAKE_SOURCE_DIR}/Settings.cmake" )
    include( ${CMAKE_SOURCE_DIR}/Settings.cmake )
else()
    message( STATUS "loading default settings file ${CMAKE_SOURCE_DIR}/Settings.cmake.original" )
    include( ${CMAKE_SOURCE_DIR}/Settings.cmake.original )
endif()

if( USE_BLOCKSIZE_CALCULATOR )
    set( DEMOTRACK_USE_BLOCKSIZE_CALCULATOR "1" )
else()
    set( DEMOTRACK_USE_BLOCKSIZE_CALCULATOR "0" )
endif()

configure_file( config.h.template ${CMAKE_SOURCE_DIR}/include/config.h @ONLY )

##### Detect enabled target arch ####
# This shhould passed as arg, so here it hardcoded. FIX-ME
set(ENABLE_HIP "ON")

#CUDA
if(ENABLE_CUDA)
  #Do CUDA stuff
endif()

#OpenCL
if(ENABLE_OPENCL)
  #Do OCL stuff
  find_package(OpenCL)
endif()

#HIP
if(ENABLE_HIP)
  # -------------------------------------------------------------------------------
  # Enable HIP support
  
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
  
   if(HIP_FOUND)
    set(HIP_ENABLED ON)
    message(STATUS "HIP Found (${hip_HIPCC_EXECUTABLE})")
    #set(HIP_CMAKE_CXX_FLAGS "-fcuda-flush-denormals-to-zero -Wno-invalid-command-line-argument -Wno-unused-command-line-argument -Wno-invalid-constexpr -Wno-ignored-optimization-argument -Wno-unused-private-field")

    if(DEFINED HIP_AMDGPUTARGET)
      if("${HIP_PLATFORM}" STREQUAL hcc)	    
        set(HIP_CMAKE_CXX_FLAGS "${HIP_CMAKE_CXX_FLAGS} --amdgpu-target=${HIP_AMDGPUTARGET}")
      elseif("${HIP_PLATFORM}" STREQUAL nvcc)
        set(HIP_CMAKE_CXX_FLAGS "${HIP_CMAKE_CXX_FLAGS} -arch=${HIP_AMDGPUTARGET}")
      endif()
    else()
      if("${HIP_PLATFORM}" STREQUAL hcc)	    
	set(HIP_AMDGPUTARGET "gfx906")
        set(HIP_CMAKE_CXX_FLAGS "${HIP_CMAKE_CXX_FLAGS} --amdgpu-target=gfx906")
      elseif("${HIP_PLATFORM}" STREQUAL nvcc)
	set(HIP_AMDGPUTARGET "compute_70")
        set(HIP_CMAKE_CXX_FLAGS "${HIP_CMAKE_CXX_FLAGS} -arch=compute_70")
      endif()	
    endif()
  endif()

  # ******************************************************************************
  # demo applications
  
  #set( CLANG_MCPU_ARCH "gfx900" ) #--amdgpu-target=${CLANG_MCPU_ARCH}
  
  set( CMAKE_CXX_STANDARD_REQUIRED ON )

  if("${HIP_PLATFORM}" STREQUAL hcc)
    set( CMAKE_POSITION_INDEPENDENT_CODE ON )
    set( CMAKE_CXX_STANDARD 11 )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${HIP_CMAKE_CXX_FLAGS} -fno-gpu-rdc -fno-fast-math")
    message(STATUS "Building demo with HIP support ${HIP_AMDGPUTARGET}")
  elseif("${HIP_PLATFORM}" STREQUAL nvcc)
    set(HIP_NVCC_FLAG "-x cu -arch=${HIP_AMDGPUTARGET}")
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${HIP_NVCC_FLAG}")  
    message(STATUS "Building demo CUDA with HIP support ${HIP_AMDGPUTARGET}")
  endif()

  add_subdirectory(src)
  include_directories( ${CMAKE_SOURCE_DIR}/include )
  add_subdirectory( tools )
  
endif()
