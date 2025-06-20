#!/bin/bash
for x in {0..9}
  do
  for i in {0..9}
    do
	echo "=========================>/work1/scquant/ParaEngine-git-COVID2019/Workdir-wb97xd-6311gss-12A-ACE2_Ab-Wild-type_"$x"_$i"
	cd /work1/scquant/ParaEngine-git-COVID2019/Workdir-wb97xd-6311gss-12A-ACE2_Ab-Wild-type_"$x"_$i
	tail -n 1 */*.log > tail.txt
	sed -e "/Normal/d" -e "/^\s*$/d" tail.txt > tail1.txt                        #Just the err lines are left.
	grep -B1 -v "==>" tail1.txt > err.txt                                        #get the err lines and their gjf names.
	grep "==>" err.txt > err2.txt                                                #get the gjf names of the err lines.
	sed -e "s/\s//g" -e "s/==>//g" -e "s/.log<==//g" err2.txt  > file.txt        #get the gjf to submit a task
	head -n 3 err.txt
	echo "========================================================================================================================"
    done
  done
