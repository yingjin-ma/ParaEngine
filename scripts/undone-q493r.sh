#!/bin/bash
for x in {0..9}
do
    for i in {0..9}
    do
	echo "/work1/scquant/ParaEngine-git-COVID2019/Workdir-B2PLYP-631gss-GD3-12A-ACE2_Ab-Omicron_Q493_"$x"_$i"
	cd /work1/scquant/ParaEngine-git-COVID2019/Workdir-B2PLYP-631gss-GD3-12A-ACE2_Ab-Omicron_Q493_"$x"_$i
	ls Frag*.log > log.txt
	grep "Frag" log.txt >existlog.txt
	if [ -s existlog.txt ]
	then
	    echo "!!!!!!!!!!!!!!!!!!!There are already gjf tasks here.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	else
	    ls *.gjf > gjf.txt
	    grep ".gjf" gjf.txt > undone.txt
	    if [ -s undone.txt ]
	    then
	    	sed -i "s/.gjf//g" undone.txt
	    	head -n 4 undone.txt
	    	for j in `cat undone.txt`
	        do
	      	     bsub -n 24 -q cpuII -e err.log -o out.log g09 $j.gjf
	        done
	    else
	        echo "----------------------All the gif are done for the first time.-------------------------"
	    fi
	fi
    done
done

	
