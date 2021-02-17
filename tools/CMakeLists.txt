add_executable( create_particle_data create_particle_data.cpp )

if("${HIP_PLATFORM}" STREQUAL hcc)
set_target_properties( create_particle_data PROPERTIES CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON POSITION_INDEPENDENT_CODE ON )
endif()


add_executable( compare_particle_data compare_particle_data.cpp )

if("${HIP_PLATFORM}" STREQUAL hcc)
set_target_properties( compare_particle_data PROPERTIES CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON POSITION_INDEPENDENT_CODE ON )
endif()
