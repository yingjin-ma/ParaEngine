#!/bin/bash
for x in {0..9}
do
for i in {0..9}
do
    FILE=/work1/scquant/ParaEngine-git-COVID2019/Workdir-wb97xd-6311gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i
    if [ -d $FILE ]
    then
	echo "===========================================================>Workdir-wb97xd-6311gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i<==============================================="
	cd /work1/scquant/ParaEngine-git-COVID2019/Workdir-wb97xd-6311gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i
	grep -B1 "^[[:space:]]*$" tail.txt > space.txt	
	grep "==>" space.txt > space_name.txt
	sed -i "s/==> //g" space_name.txt
	sed -i "s/.log <==//g" space_name.txt
	for j in `cat space_name.txt`
	do
	    ls -l $j.log
	    bsub -n 24 -q cpuII -e err_"$x"_$i.log -o out_"$x"_$i.log -R "select[hname!='m4510']" g09 $j.gjf
	done
    else
	echo "==============================================================$FILE not exist.========================================================="
    fi
done
done

