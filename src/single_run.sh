#!./bin/bash
input="input.dat"
summary="../results/summary.dat"
temp="temp.txt"
start=`date +%s`
shock_str=""
cfl_str=""
#touch $summary
#if [ -f $summary ]; then
#  rm -f "$summary"
#fi
imax=513 #16 32 64 128 256 512
p0=300.0
T0=600.0
flux=2
shock=0 #0 1
ramp=0 #0 1
prat=0.4 #0.1 0.5 0.6
cfl=0.5 #0.1 0.5 0.9
maxk=100000
k2=0.5 #0.5 0.4 0.3 0.25
k4=0.03125 #0.03125 0.02 0.015625 
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
echo "p0        $p0" >> $input
echo "T0        $T0" >> $input
echo "flux      $flux" >> $input
echo "shock     $shock" >> $input
echo "ramp      $ramp" >> $input
echo "p_rat     $prat" >> $input
echo "CFL       $cfl" >> $input
echo "maxk      $maxk" >> $input
echo "k2        $k2" >> $input
echo "k4        $k4" >> $input
./../build/bin/test_program
#cat "$temp" >> $summary
rm -f temp.txt
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"