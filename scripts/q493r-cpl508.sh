#!/bin/bash
rm /work1/scquant/ParaEngine-git-COVID2019/Workdir-MP2-631gss-12A-ACE2_Ab-Omicron_Q493R*/checkerr.txt
for x in {0..9}
do
	for i in {0..9}
		do
		echo  "----------------Workdir-MP2-631gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i-----------------------"
		FILE=/work1/scquant/ParaEngine-git-COVID2019/Workdir-MP2-631gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i
		if [ -d "$FILE" ]
		then
    			#echo "*******************$FILE is a directory.*****************"
			cd /work1/scquant/ParaEngine-git-COVID2019/Workdir-MP2-631gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i
			for k in `cat file.txt`
			do
                            echo "==>$k" >> checkerr.txt
			    tail -n 4 $k.log >> checkerr.txt
			done
			    grep "/soft/apps/Paid/g09/l508.exe" checkerr.txt > exist508.txt
			    grep "/soft/apps/Paid/g09/l801.exe" checkerr.txt > exist801.txt
			    grep "/soft/apps/Paid/g09/l202.exe" checkerr.txt > exist202.txt
			if [ -s exist508.txt ];then
			    grep -B2 "/soft/apps/Paid/g09/l508.exe" checkerr.txt > sub508.txt
			    grep "==>" sub508.txt > sub508name.txt
			    sed -i "s/==>//g" sub508name.txt
			    head -n 4 sub508name.txt
			    cp sub508name.txt /work1/scquant/ParaEngine-git-COVID2019/Workdir-MP2-631gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_"$i"_l508.txt
			else
			    echo "==========================There are no l508 tasks.====================================="
			fi
			if [ -s exist202.txt ];then 
			    grep -B2 "/soft/apps/Paid/g09/l202.exe" checkerr.txt > sub202.txt
                            grep "==>" sub202.txt > sub202name.txt
                            sed -i "s/==>//g" sub202name.txt
                            head -n 4 sub202name.txt
                            cp sub202name.txt /work1/scquant/ParaEngine-git-COVID2019/Workdir-MP2-631gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_"$i"_l202.txt
			else
			    echo "==========================There are no l202 tasks.====================================="
			fi
			if [ -s exist801.txt ];then
			    grep -B2 "/soft/apps/Paid/g09/l801.exe" checkerr.txt > sub801.txt
                            grep "==>" sub801.txt > sub801name.txt
                            sed -i "s/==>//g" sub801name.txt
                            head -n 4 sub801name.txt
#			    for j in `cat sub801name.txt`
#			    do
#				sed -i "s/scf\(maxcyc\=100\)/scf\(maxcyc\=100\) iop\(8\/11\=1\)/g" $j.gjf
#				sed -i "s/scf\(maxcyc\=100\,qc\)/scf\(maxcyc\=100\,qc\) iop\(8\/11\=1\)/g" $j.gjf
#				head -n 4 $j.gjf
#				bsub -n 24 -q cpuII -e err.log -o out.log g09 $j.gjf
#			    done
                            cp sub801name.txt /work1/scquant/ParaEngine-git-COVID2019/Workdir-MP2-631gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_"$i"_l801.txt

			else
			    echo "*****************There are no l801 tasks.********************"
			fi
		else
			echo "!!!!!!!!!!!!!!!!!!!!!!$FILE not exist.!!!!!!!!!!!!!!!!!!!!!!!"
		fi
	done
done





