include_directories( ${CMAKE_SOURCE_DIR}/include )


set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set( COMMON_HEADERS 
	${PROJECT_SOURCE_DIR}/include/definitions.h 
	${PROJECT_SOURCE_DIR}/include/beam_elements.h 
	${PROJECT_SOURCE_DIR}/include/beamfields.h 
	${PROJECT_SOURCE_DIR}/include/config.h
	${PROJECT_SOURCE_DIR}/include/lattice.h 
	${PROJECT_SOURCE_DIR}/include/particle.h 
	${PROJECT_SOURCE_DIR}/include/demo01_sc0.h 
	${PROJECT_SOURCE_DIR}/include/demo01_sc1.h 
	${PROJECT_SOURCE_DIR}/include/demo02_sc0.h 
	${PROJECT_SOURCE_DIR}/include/demo02_sc1.h)


# ------------------------------------------------------------------------------
# demo01: particles stored in global memory, lattice stored in global memory

HIP_ADD_EXECUTABLE( demotrack main.cpp ${COMMON_HEADERS} )
target_link_libraries( demotrack PRIVATE ${DEMOTRACK_LINK_LIBRARIES} )
target_compile_options( demotrack PRIVATE )
target_compile_definitions( demotrack PRIVATE -DDEMOTRACK_ENABLE_BEAMFIELDS=1 )
