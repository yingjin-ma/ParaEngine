#!/bin/bash
for x in {0..9}
do
    for i in {0..9}
    do
	echo "/work1/scquant/ParaEngine-git-COVID2019/Workdir-B2PLYP-631gss-GD3-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i"
	cd /work1/scquant/ParaEngine-git-COVID2019/Workdir-B2PLYP-631gss-GD3-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i
	ls Frag*.log > log.txt
	grep "Frag" log.txt >existlog.txt
	if [ -s existlog.txt ]
	then
	    echo "!!!!!!!!!!!!!!!!!!!There are already gjf tasks here.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!=======Workdir-B2PLYP-631gss-GD3-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i" 
	fi
    done
done

	
