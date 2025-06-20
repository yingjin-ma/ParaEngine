#!/bin/bash
for x in {0..9}
do
for i in {0..9}
do
    FILE=/work1/scquant/ParaEngine-git-COVID2019/Workdir-wb97xd-6311gss-12A-ACE2_Ab-Wild-type_"$x"_$i
    if [ -d "$FILE" ]
	then
		echo "Workdir-wb97xd-6311gss-12A-ACE2_Ab-Wild-type_"$x"_$i exists."
	else
		echo "!!!!!!!!!!!!!!!!!!!!Workdir-wb97xd-6311gss-12A-ACE2_Ab-Wild-type_"$x"_$i not exist.!!!!!!!!!!!!!!!!!!!!!!!"
    fi
done
done



