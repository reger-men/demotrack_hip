
step=32

for i in {1..32}
do 
  threads=$((step*i))
  echo "Run $threads Threads" &>> tun.log

  #Replace threads
  sed -i "9s/DEMOTRACK_DEFAULT_BLOCK_SIZE.*/DEMOTRACK_DEFAULT_BLOCK_SIZE $threads/" ../include/config.h

  #Rebuild
  make -j

  #Run
  HSA_NO_SCRATCH_RECLAIM=0 HIP_VISIBLE_DEVICES=1 ./demotrack 102400 10 ../data/full_SIS_lattice/sis100_particles_flat.bin ../data/full_SIS_lattice/sis100_elements_flat.bin -k 3 -e | grep "Elapsed time          :" >> tun.log

  echo "-----------------------------------" &>> tun.log
done

