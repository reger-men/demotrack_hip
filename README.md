# demotrack_hip
Hip Version of the demotrack repository

## Build instructions

Building for NVIDIA:
```
git clone https://github.com/reger-men/demotrack_hip
cd demotrack_hip
mkdir build
cd build
export HIP_PLATFORM=nvcc && CXX=hipcc cmake -DHIP_AMDGPUTARGET=compute_70  ..
make
```

Building for AMD:
```
git clone https://github.com/reger-men/demotrack_hip
cd demotrack_hip
mkdir build
cd build
export HIP_PLATFORM=hcc && CXX=hipcc cmake -DHIP_AMDGPUTARGET=gfx906 ..
make
Scanning dependencies of target demo02_sc1
[ 12%] Building CXX object CMakeFiles/demo02_sc1.dir/demo02.cpp.o
[ 25%] Linking HIP executable demo02_sc1
[ 25%] Built target demo02_sc1
Scanning dependencies of target demo01_sc1
[ 37%] Building CXX object CMakeFiles/demo01_sc1.dir/demo01.cpp.o
[ 50%] Linking HIP executable demo01_sc1
[ 50%] Built target demo01_sc1
Scanning dependencies of target demo02_sc0
[ 62%] Building CXX object CMakeFiles/demo02_sc0.dir/demo02.cpp.o
[ 75%] Linking HIP executable demo02_sc0
[ 75%] Built target demo02_sc0
Scanning dependencies of target demo01_sc0
[ 87%] Building CXX object CMakeFiles/demo01_sc0.dir/demo01.cpp.o
[100%] Linking HIP executable demo01_sc0
[100%] Built target demo01_sc0
``` 

Run the test:
``` 
# space charge enabled with globale memory access
./demotrack 102400 100 ../data/lhc02/lhc_no_bb_particles.bin ../data/lhc02/lhc_no_bb_elements.bin -k 1
./demotrack 102400 100 ../data/full_SIS_lattice/sis100_particles_flat.bin ../data/full_SIS_lattice/sis100_elements_flat.bin -k 1

# space charge enabled with local memory access
./demotrack 102400 100 ../data/lhc02/lhc_no_bb_particles.bin ../data/lhc02/lhc_no_bb_elements.bin -k 3
./demotrack 102400 100 ../data/full_SIS_lattice/sis100_particles_flat.bin ../data/full_SIS_lattice/sis100_elements_flat.bin -k 3
``` 

