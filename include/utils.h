#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <string>
#include <vector>

bool event_timing  = false;
int kernel_version = 0;

char pci_bus_id_str[] =
{
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0'
};

/********************* UTILS *******************************************/
#define GPUCheckFail(x) GPUCheckFailF(x, __FILE__, __LINE__)
#define GPUError(...) GPUInfo(__VA_ARGS__)
#define GPUInfo(string, ...)            \
{                                     \
  printf(string "\n", ##__VA_ARGS__); \
}

static int GPUCheckFailF(const long long int error, const char* file, int line)
{
  // Check for HIP Error and in the case of an error display the corresponding error string
  if (error == hipSuccess) {
    return (0);
  }
  GPUError("HIP Error: %lld / %s (%s:%d)", error, hipGetErrorString((hipError_t)error), file, line);
  return 1;
}

/* ***************************************************************** */
namespace dt = demotrack;
dt::uint64_type NUM_PARTICLES = 50 * 1024;
dt::int64_type  TRACK_UNTIL_TURN = 1000;
std::string path_to_lattice_data = std::string{};
std::string path_to_particle_data = std::string{};
std::string path_to_output_data = std::string{};

void parseArguments(int argc, char *argv[])
{
  // Extract optional args
  int new_argc = argc;
  for (int i = 1; i < argc; i++)
  { 
    if (!std::string("--event-timing").compare(argv[i]) ||
             !std::string("-e").compare(argv[i]))
    {
      event_timing = true;
      new_argc--;
    }else if(!std::string("-k").compare(argv[i])){
      kernel_version = atoi(argv[++i]);
      new_argc=new_argc-2;
    }
    
  }
  argc = new_argc;

  if( argc >= 2 )
  {
      NUM_PARTICLES = std::stoi( argv[ 1 ] );
      if( argc >= 3 )
      {
          TRACK_UNTIL_TURN = std::stoi( argv[ 2 ] );
          if( argc >= 4 )
          {
              path_to_particle_data = std::string{ argv[ 3 ] };
              if( path_to_particle_data.compare( "default" ) == 0 ) {
                  path_to_particle_data.clear(); }
              if( argc >= 5 )
              {
                  path_to_lattice_data = std::string{ argv[ 4 ] };
                  if( path_to_lattice_data.compare( "default" ) == 0 ) {
                      path_to_lattice_data.clear(); }
                  if( argc >= 6 )
                  {
                      path_to_output_data = std::string{ argv[ 5 ] };
                      if( path_to_output_data.compare( "none" ) == 0 ) {
                          path_to_output_data.clear(); }
                  }
              }
          }
      }
  }
  else
  {
      std::cout << "Usage : " << argv[ 0 ]
                << " [NUM_PARTICLES] [TRACK_UNTIL_TURN]"
                << " [PATH_TO_PARTICLE_DATA] [PATH_TO_LATTICE_DATA]"
                << " [PATH_TO_OUTPUT_DATA]\r\n"
                << std::endl;
  }

}
