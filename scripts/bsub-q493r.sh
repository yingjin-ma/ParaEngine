#!/bin/bash

for x in {7,8,9}
do
for i in {0..9}
do
    bsub -n 240 -q cpuII -e err_"$x"_$i.log -o out_"$x"_$i.log -R "select[hname!='m4510']" ./RUN_ParaEngine_julia_ERA.sh input-wb97xd-6311gss-Q493R_"$x"_$i-STDDLB-10nodes
    sleep 0.5
done
done
