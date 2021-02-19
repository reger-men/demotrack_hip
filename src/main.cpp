#include <hip/hip_runtime.h>
#ifndef __HIP_PLATFORM_NVCC__
#include "hip/hip_ext.h"
#endif

#include "definitions.h"
#include "config.h"
#include "particle.h"
#include "beam_elements.h"
#if defined( DEMOTRACK_ENABLE_BEAMFIELDS ) && ( DEMOTRACK_ENABLE_BEAMFIELDS == 1 )
#include "beamfields.h"
#endif /* DEMOTRACK_ENABLE_BEAMFIELDS == 1 */
#include "lattice.h"
#include "utils.h"
#include "demo01_sc0.h"
#include "demo01_sc1.h"
#include "demo02_sc0.h"
#include "demo02_sc1.h"



int main( int argc, char* argv[] )
{
    parseArguments(argc, argv);

    /********************************************************************** */
    /* Debug functions */
    hipDeviceProp_t hipDeviceProp;
    int mDeviceId = 0; //TODO fix-me
    int debugLevel = 1;

    GPUCheckFail(hipGetDeviceProperties(&hipDeviceProp, mDeviceId));


    if (debugLevel >= 1) {
      GPUInfo("Using HIP Device %s with Properties:", hipDeviceProp.name);
      GPUInfo("\ttotalGlobalMem = %lld", (unsigned long long int)hipDeviceProp.totalGlobalMem);
      GPUInfo("\tsharedMemPerBlock = %lld", (unsigned long long int)hipDeviceProp.sharedMemPerBlock);
      GPUInfo("\tregsPerBlock = %d", hipDeviceProp.regsPerBlock);
      GPUInfo("\twarpSize = %d", hipDeviceProp.warpSize);
      GPUInfo("\tmaxThreadsPerBlock = %d", hipDeviceProp.maxThreadsPerBlock);
      GPUInfo("\tmaxThreadsDim = %d %d %d", hipDeviceProp.maxThreadsDim[0], hipDeviceProp.maxThreadsDim[1], hipDeviceProp.maxThreadsDim[2]);
      GPUInfo("\tmaxGridSize = %d %d %d", hipDeviceProp.maxGridSize[0], hipDeviceProp.maxGridSize[1], hipDeviceProp.maxGridSize[2]);
      GPUInfo("\ttotalConstMem = %lld", (unsigned long long int)hipDeviceProp.totalConstMem);
      GPUInfo("\tmajor = %d", hipDeviceProp.major);
      GPUInfo("\tminor = %d", hipDeviceProp.minor);
      GPUInfo("\tclockRate = %d", hipDeviceProp.clockRate);
      GPUInfo("\tmemoryClockRate = %d", hipDeviceProp.memoryClockRate);
      GPUInfo("\tmultiProcessorCount = %d", hipDeviceProp.multiProcessorCount);
      GPUInfo("\tPCIBusID = %d", hipDeviceProp.pciBusID);
      GPUInfo(" ");
    }
    

    /* ********************************************************************* */
    /* Prepare particle data: */

    std::vector< dt::Particle > particles_host;
    dt::Particles_load( particles_host, NUM_PARTICLES, path_to_particle_data );

    /* ********************************************************************* */
    /* Prepare lattice / machine description: */

    std::vector< double > lattice_host;
    dt::uint64_type const LATTICE_SIZE = dt::load_lattice(
        lattice_host, path_to_lattice_data );

    /* ********************************************************************** */
    /* Allocate buffers on the device */

    dt::Particle* particles_dev = nullptr;
    double* lattice_dev = nullptr;

    GPUCheckFail(hipMalloc( ( void** )&particles_dev, sizeof( dt::Particle ) * NUM_PARTICLES ));
    GPUCheckFail(hipMalloc( ( void** )&lattice_dev, LATTICE_SIZE * sizeof( double )));

    /* Copy particle and lattice data to device */

    GPUCheckFail(hipMemcpy( lattice_dev, lattice_host.data(), LATTICE_SIZE * sizeof( double ), ::hipMemcpyHostToDevice ));
    GPUCheckFail(hipMemcpy( particles_dev, particles_host.data(), particles_host.size() * sizeof( dt::Particle ), hipMemcpyHostToDevice));

    /* ******************************************************************** */
    /* Estimate block size */

    int BLOCK_SIZE = 0;
    int MIN_GRID_SIZE = 0;
    int MAX_GRID_SIZE = 0;

#if defined( DEMOTRACK_HIP_CALCULATE_BLOCKSIZE ) && ( DEMOTRACK_HIP_CALCULATE_BLOCKSIZE == 1 )
    GPUFailedMsg(hipOccupancyMaxPotentialBlockSize(&MIN_GRID_SIZE, &BLOCK_SIZE, Track_particles_until_turn, 0, 0));
    GPUFailedMsg(hipOccupancyMaxActiveBlocksPerMultiprocessor(&MAX_GRID_SIZE, Track_particles_until_turn, BLOCK_SIZE, 0));
#elif defined( DEMOTRACK_DEFAULT_BLOCK_SIZE ) && ( DEMOTRACK_DEFAULT_BLOCK_SIZE > 0 )
    BLOCK_SIZE = DEMOTRACK_DEFAULT_BLOCK_SIZE;
#else
    BLOCK_SIZE = 1;
#endif /* DEMOTRACK_HIP_CALCULATE_BLOCKSIZE */

    assert( BLOCK_SIZE > 0 );
    int const GRID_SIZE = ( NUM_PARTICLES + BLOCK_SIZE - 1 ) / BLOCK_SIZE;

    GPUInfo("Kernel                : %20s Block size: %4d, Maximum active blocks: %3d, Suggested blocks: %3d", "Track_particles_until_turn",BLOCK_SIZE,MAX_GRID_SIZE,MIN_GRID_SIZE);

    /* ******************************************************************** */
    /* Run kernel: */

    GPUCheckFail(hipDeviceGetPCIBusId( pci_bus_id_str, 32, mDeviceId));

    GPUInfo("number of particles   : %8llu \nnumber of turns       : %8lld", NUM_PARTICLES, TRACK_UNTIL_TURN);

    if( !path_to_particle_data.empty() )
    {
        GPUInfo("particle data         : %50s", path_to_particle_data.c_str());
    }
    else
    {
        GPUInfo("particle data         : %s", "generated");
    }

    if( !path_to_lattice_data.empty() )
    {
        GPUInfo("lattice               : %50s", path_to_lattice_data.c_str());
    }
    else
    {
        GPUInfo("lattice               : %s", "generated fodo lattice");
    }

#if defined( DEMOTRACK_ENABLE_BEAMFIELDS ) && ( DEMOTRACK_ENABLE_BEAMFIELDS == 1 )
    GPUInfo("space-charge enabled  : %s", "true");
#else
    GPUInfo("space-charge enabled  :  %s", "false");
#endif /* SC emabled */

    if( !path_to_output_data.empty() )
    {
        GPUInfo("path to output data : %50s", path_to_output_data.c_str());
    }

    GPUInfo("DEVICE                : %15s (%10s)",  pci_bus_id_str, hipDeviceProp.name);
    GPUInfo("NUM_OF_BLOCKS         : %6d", GRID_SIZE);
    GPUInfo("MIN_GRID_SIZE         : %6d", MIN_GRID_SIZE);
    GPUInfo("THREADS_PER_BLOCK     : %6d", BLOCK_SIZE);

    /* Prepare cuda events to estimate the elapsed wall time */

    hipEvent_t start;
    GPUCheckFail(hipEventCreate( &start ));

    hipEvent_t stop;
    GPUCheckFail(hipEventCreate( &stop ));

    GPUInfo("KERNEL VERSION        : %6d", kernel_version);
    GPUInfo("Event Timing          : %6d", event_timing);

    /******************************* Run kernel **********************************/
    if(!event_timing){
      GPUCheckFail(hipEventRecord( start ));

      if(kernel_version == 0){
        hipLaunchKernelGGL(Track_particles_until_turn00_00, dim3(GRID_SIZE), dim3(BLOCK_SIZE ), 0, 0,
          particles_dev, NUM_PARTICLES, lattice_dev, LATTICE_SIZE, TRACK_UNTIL_TURN );
      }else if(kernel_version == 1){
        hipLaunchKernelGGL(Track_particles_until_turn00_01, dim3(GRID_SIZE), dim3(BLOCK_SIZE ), 0, 0,
          particles_dev, NUM_PARTICLES, lattice_dev, LATTICE_SIZE, TRACK_UNTIL_TURN );	     
      }else if(kernel_version == 2){
        hipLaunchKernelGGL(Track_particles_until_turn01_00, dim3(GRID_SIZE), dim3(BLOCK_SIZE ), 0, 0,
          particles_dev, NUM_PARTICLES, lattice_dev, LATTICE_SIZE, TRACK_UNTIL_TURN );
      }else if(kernel_version == 3){
        hipLaunchKernelGGL(Track_particles_until_turn01_01, dim3(GRID_SIZE), dim3(BLOCK_SIZE ), 0, 0,
          particles_dev, NUM_PARTICLES, lattice_dev, LATTICE_SIZE, TRACK_UNTIL_TURN );
      }else{
        GPUError("Error: Kernel version %d not found!", kernel_version);
        return 1;
      }

      GPUCheckFail(hipDeviceSynchronize());

      /* Estimate wall time */
      GPUCheckFail(hipEventRecord( stop ));
      GPUCheckFail(hipEventSynchronize( stop ));
    }else{
#ifndef __HIP_PLATFORM_NVCC__	    
      if(kernel_version == 0){
      hipExtLaunchKernelGGL(Track_particles_until_turn00_00, dim3(GRID_SIZE), dim3(BLOCK_SIZE ), 0, 0,  start, stop, 0,
			       particles_dev, NUM_PARTICLES, lattice_dev, LATTICE_SIZE, TRACK_UNTIL_TURN);
      }else if(kernel_version == 1){
         hipExtLaunchKernelGGL(Track_particles_until_turn00_01, dim3(GRID_SIZE), dim3(BLOCK_SIZE ), 0, 0,  start, stop, 0,
                               particles_dev, NUM_PARTICLES, lattice_dev, LATTICE_SIZE, TRACK_UNTIL_TURN);	      
      }else if(kernel_version == 2){
         hipExtLaunchKernelGGL(Track_particles_until_turn01_00, dim3(GRID_SIZE), dim3(BLOCK_SIZE ), 0, 0,  start, stop, 0,
                               particles_dev, NUM_PARTICLES, lattice_dev, LATTICE_SIZE, TRACK_UNTIL_TURN);
      }else if(kernel_version == 3){
         hipExtLaunchKernelGGL(Track_particles_until_turn01_01, dim3(GRID_SIZE), dim3(BLOCK_SIZE ), 0, 0,  start, stop, 0,
                               particles_dev, NUM_PARTICLES, lattice_dev, LATTICE_SIZE, TRACK_UNTIL_TURN);
      }else{
        GPUError("Error: Kernel version %d not found!", kernel_version);
	return 1;
      }

      GPUCheckFail(hipEventSynchronize(stop));
#endif      
    }


    float wtime = 0.0;
    GPUCheckFail(hipEventElapsedTime( &wtime, start, stop ));

    std::cout << "-------------------------------------------------------\r\n"
              << "Elapsed time          : " << wtime << " msec total \r\n"
              << "                      : " << wtime / (
                std::max( NUM_PARTICLES * TRACK_UNTIL_TURN,
                          dt::uint64_type{ 1 } ) )
              << " msec / particle / turn\r\n";

    /* Fetch data */

    GPUCheckFail(hipMemcpy( particles_host.data(), particles_dev,
                           particles_host.size() * sizeof( dt::Particle ),
                           hipMemcpyDeviceToHost ));

    /* ********************************************************************* */
    /* Verify tracking results */

    dt::uint64_type num_active_particles = 0u;
    dt::uint64_type num_lost_particles = 0u;

    for( auto& p : particles_host )
    {
        if( ( p.state == 1 ) && ( p.at_turn == TRACK_UNTIL_TURN ) )
        {
            ++num_active_particles;
        }
        else if( ( p.state == 0 ) && ( p.at_turn < TRACK_UNTIL_TURN ) )
        {
            ++num_lost_particles;
        }
        else
        {
            std::cerr << "illegal particle id = " << p.id
                      << ", at_turn = " << p.at_turn
                      << ", at_element = " << p.at_element
                      << ", state = " << p.state << std::endl;
        }
    }

    std::cout << "-------------------------------------------------------\r\n"
              << "num lost particles    : " << num_lost_particles << "\r\n"
              << "num active particles  : " << num_active_particles << "\r\n"
              << std::endl;

    if( !path_to_output_data.empty() )
    {
        FILE* fp = std::fopen( path_to_output_data.c_str(), "wb" );
        double const temp = static_cast< double >( particles_host.size() );
        auto ret = std::fwrite( &temp, sizeof( double ), 1u, fp );
        bool success = ( ret == 1 );

        for( auto const& p : particles_host )
        {
            ret = std::fwrite( &p, sizeof( dt::Particle ), 1u, fp );
            success &= ( ret == 1 );
        }

        if( success )
        {
            std::cout << "Written particle state to " << path_to_output_data
                      << "\r\n" << std::endl;
        }

        std::fflush( fp );
        std::fclose( fp );
    }

    /* ********************************************************************* */
    /* Cleaning up, Freeing resources */

    hipFree( lattice_dev );
    lattice_dev = nullptr;

    hipFree( particles_dev );
    particles_dev = nullptr;

    hipEventDestroy( start );
    hipEventDestroy( stop );

    return 0;
}
