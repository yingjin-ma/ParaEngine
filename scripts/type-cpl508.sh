#!/bin/bash
rm /work1/scquant/ParaEngine-git-COVID2019/Workdir-MP2-631gss-12A-ACE2_Ab-Wild-type*/checkerr.txt
for x in {0..9}
do
	for i in {0..9}
		do
		echo  "----------------Workdir-MP2-631gss-12A-ACE2_Ab-Wild-type_"$x"_$i-----------------------"
		FILE=/work1/scquant/ParaEngine-git-COVID2019/Workdir-MP2-631gss-12A-ACE2_Ab-Wild-type_"$x"_$i
		if [ -d "$FILE" ]
		then
    			#echo "*******************$FILE is a directory.*****************"
			cd /work1/scquant/ParaEngine-git-COVID2019/Workdir-MP2-631gss-12A-ACE2_Ab-Wild-type_"$x"_$i
			for k in `cat file.txt`
			do
                            echo "==>$k" >> checkerr.txt
			    tail -n 4 $k.log >> checkerr.txt
			done
			    grep "/soft/apps/Paid/g09/l508.exe" checkerr.txt > exist508.txt
			    grep "/soft/apps/Paid/g09/l202.exe" checkerr.txt > exist202.txt
			if [ -s exist508.txt ];then
			    grep -B2 "/soft/apps/Paid/g09/l508.exe" checkerr.txt > sub508.txt
			    grep "==>" sub508.txt > sub508name.txt
			    sed -i "s/==>//g" sub508name.txt
			    cp sub508name.txt /work1/scquant/ParaEngine-git-COVID2019/Workdir-MP2-631gss-12A-ACE2_Ab-Wild-type_"$x"_"$i"_l508.txt
			else
			    echo "======================================There are no l508 tasks.=============================="
			fi
			if [ -s exist202.txt ];then
			    grep -B2 "/soft/apps/Paid/g09/l202.exe" checkerr.txt > sub202.txt
                            grep "==>" sub202.txt > sub202name.txt
                            sed -i "s/==>//g" sub202name.txt
                            cp sub202name.txt /work1/scquant/ParaEngine-git-COVID2019/Workdir-MP2-631gss-12A-ACE2_Ab-Wild-type_"$x"_"$i"_l202.txt
			else
			    echo "*****************There are no l202 tasks.********************"
			fi
		else
			echo "!!!!!!!!!!!!!!!!!!!!!!$FILE not exist.!!!!!!!!!!!!!!!!!!!!!!!"
		fi
	done
done





