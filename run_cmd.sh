##### Demo01 -k 1 ###########################################################

## fc1 21344.5 msec | fc0 29494.5 
HIP_VISIBLE_DEVICES=3 ./demotrack 102400 100 ../data/lhc02/lhc_no_bb_particles.bin ../data/lhc02/lhc_no_bb_elements.bin -k 1

## fc1 47117.8 msec | fc0 99784.2
HIP_VISIBLE_DEVICES=3 ./demotrack 102400 100 ../data/full_SIS_lattice/sis100_particles_flat.bin ../data/full_SIS_lattice/sis100_elements_flat.bin -k 1


##### Demo02 -k 3 ###########################################################

## fc1 9615.06 msec | fc0 9702.61
HIP_VISIBLE_DEVICES=3 ./demotrack 102400 100 ../data/lhc02/lhc_no_bb_particles.bin ../data/lhc02/lhc_no_bb_elements.bin -k 3

## fc1 31579 msec   | fc0 31755.2 
HIP_VISIBLE_DEVICES=3 ./demotrack 102400 100 ../data/full_SIS_lattice/sis100_particles_flat.bin ../data/full_SIS_lattice/sis100_elements_flat.bin -k 3
