#!./bin/bash
filename="input.dat"
start=`date +%s`
shock_str=""
cfl_str=""
for imax in 32 #16 32 64 128 256 512
do
  for shock in 0 #0 1
  do
    for ramp in 0 #0 1
    do
      for prat in 0.1 #0.1 0.5 0.6
      do
        for k2 in 0.5 #0.5 0.4 0.3 0.25
        do
          for k4 in 0.03125 #0.03125 0.02 0.015625 
          do
            for cfl in 0.1 #0.1 0.5 0.9
            do
              if [ $shock -eq 0 ]
              then
                shock_str="Isentropic"
              else
                shock_str="Normal Shock"
              fi
              if [ $ramp -eq 1 ]
              then
                cfl_str="cfl-ramp"
              else
               cfl_str="cfl-const"
              fi
echo ""
echo "################################################################################"
echo " N=$imax | $shock_str | $cfl_str | P_rat=$prat | k2=$k2 | k4=$k4 | CFL=$cfl "
echo "################################################################################"
              echo "! Input Data" > $filename
              echo "----|----|----|----|----|----|" >> $filename
              echo "imax      $imax" >> $filename
              echo "shock     $shock" >> $filename
              echo "ramp      $ramp" >> $filename
              echo "p_rat     $prat" >> $filename
              echo "CFL       $cfl" >> $filename
              echo "k2        $k2" >> $filename
              echo "k4        $k4" >> $filename
              ./../build/bin/test_program
            done
          done
        done      
      done
    done        
  done
done
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
