#!/bin/bash
for x in {0..9}
do
for i in {0..9}
do
    FILE=/work1/scquant/ParaEngine-git-COVID2019/Workdir-wb97xd-6311gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i
    if [ -d $FILE ]
    then
	echo "=========================================================Workdir-wb97xd-6311gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i================================="
	cd /work1/scquant/ParaEngine-git-COVID2019/Workdir-wb97xd-6311gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i
	cp err.txt err-left.txt
	sed -e "/File lengths/d" -e "/--/d" err-left.txt | grep -B1 -v "==>" | grep "==>" > file-left.txt
	sed -i "s/==> //g" file-left.txt
	sed -i "s/.log <==//g" file-left.txt
	head -n 4 file-left.txt
	for j in `cat file-left.txt`
	do
	    bsub -n 24 -q cpuII -e err.log -o out.log -R "select[hname!='m4510']" g09 $j.gjf
	done
    else
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Workdir-wb97xd-6311gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i not exist!!!!!!!!!!!!!!!!!!!!!!"
    fi
done
done



