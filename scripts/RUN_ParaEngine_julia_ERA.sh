#!/bin/sh

module use /home/scykd/mayj/Soft_modulefiles
module load apps/julia/1.7.2
module load NWChem/7.0.0-GNU9-openmpi4.0.0

date1=$(date +%Y%m%d%H%M%S)
echo $date1

HOSTFILE=./hosts.$LSB_JOBID
touch $HOSTFILE
let inum=0
for HOST in $LSB_MCPU_HOSTS
do
  let inum=${inum}+1
  echo "HOST : " $HOST  "  inum : " $inum
#  if [ $ODD -eq 1 ]; then
#    echo -n $HOST >> $HOSTFILE
#    HOSTLIST="$HOSTLIST,$HOST"
#    ODD=0
#  else
  let imod=$inum%2
  if [ $imod == 0 ]
  then
    #echo "$HOST" >> $HOSTFILE
    echo $imod  " ==== "  $HOST 
  else
    echo "$HOST" >> $HOSTFILE
    echo $imod  " oooo "  $HOST 
  fi
done

julia --machine-file=$HOSTFILE ParaEngine.jl $1 > out-$date1

