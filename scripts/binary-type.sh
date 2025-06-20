#!/bin/bash

rm binary-type.txt
for x in {0..9}
do
for i in {0..9}
do
    FILE=/work1/scquant/ParaEngine-git-COVID2019/Workdir-wb97xd-6311gss-12A-ACE2_Ab-Wild-type_"$x"_$i
    if [ -d $FILE ]
    then
	echo "======================================>Workdir-wb97xd-6311gss-12A-ACE2_Ab-Wild-type_"$x"_$i" >> binary-type.txt
	grep "Binary" Workdir-wb97xd-6311gss-12A-ACE2_Ab-Wild-type_"$x"_$i/err.txt >> binary-type.txt

	head -n 8 binary-type.txt  
    else
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!/work1/scquant/ParaEngine-git-COVID2019/Workdir-wb97xd-6311gss-12A-ACE2_Ab-Wild-type_"$x"_$i doesn't exist"
    fi
done
done

