#!/bin/bash
for x in {0..9}
  do
  for i in {0..9}
    do
	echo "-----/work1/scquant/ParaEngine-git-COVID2019/Workdir-wb97xd-6311gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i ----------------"
	cd /work1/scquant/ParaEngine-git-COVID2019/Workdir-wb97xd-6311gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i
	tail -n 1 */*.log > tail.txt
	sed -e "/Normal/d" -e "/^\s*$/d" tail.txt > tail1.txt
	grep -B1 -v "==>" tail1.txt > err.txt
	grep "==>" err.txt > err2.txt
	sed -e "s/\s//g" -e "s/==>//g" -e "s/.log<==//g" err2.txt  > file.txt
	head -n 3 err.txt
	echo "========================================================================================================================"
    done
  done
