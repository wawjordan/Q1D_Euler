#!./bin/bash
input="input.dat"
summary="../results/summary.dat"
temp="temp.txt"
start=`date +%s`
shock_str=""
cfl_str=""
touch $summary
if [ -f $summary ]; then
  rm -f "$summary"
fi
for shock in 1 #0 1
do
  for imax in 512 #16 32 64 128 256 512
  do
    for ramp in 0 #0 1
    do
      for prat in 0.4 #0.1 0.5 0.6
      do
        for k2 in 0.5 0.25 #0.5 0.4 0.3 0.25
        do
          for k4 in 0.03125 0.015625 #0.03125 0.02 0.015625 
          do
            for cfl in 0.1 #0.1 0.5 0.9
            do
              if [ $shock -eq 0 ]; then
                shock_str="Isentropic"
              else
                shock_str="Normal Shock"
              fi
              if [ $ramp -eq 1 ]; then
                cfl_str="cfl-ramp"
              else
               cfl_str="cfl-const"
              fi
echo ""
echo "################################################################################"
echo " N=$imax | $shock_str | $cfl_str | P_rat=$prat | k2=$k2 | k4=$k4 | CFL=$cfl "
echo "################################################################################"
              echo "! Input Data" > $input
              echo "----|----|----|----|----|----|" >> $input
              echo "imax      $imax" >> $input
              echo "shock     $shock" >> $input
              echo "ramp      $ramp" >> $input
              echo "p_rat     $prat" >> $input
              echo "CFL       $cfl" >> $input
              echo "k2        $k2" >> $input
              echo "k4        $k4" >> $input
              ./../build/bin/test_program
              #echo "N$imax" >> $summary
              cat "$temp" >> $summary
            done
          done
        done      
      done
    done        
  done
done
rm -f temp.txt
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
