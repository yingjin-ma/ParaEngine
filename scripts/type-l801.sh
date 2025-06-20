#!/bin/bash
rm /work1/scquant/ParaEngine-git-COVID2019/Workdir-wb97xd-6311gss-12A-ACE2_Ab-Wild-type*/checkerr.txt
for x in {0..9}
do
	for i in {0..9}
		do
		echo  "----------------Workdir-wb97xd-6311gss-12A-ACE2_Ab-Wild-type_"$x"_$i-----------------------"
		FILE=/work1/scquant/ParaEngine-git-COVID2019/Workdir-wb97xd-6311gss-12A-ACE2_Ab-Wild-type_"$x"_$i
		if [ -d "$FILE" ]
		then
    			#echo "*******************$FILE is a directory.*****************"
			cd /work1/scquant/ParaEngine-git-COVID2019/Workdir-wb97xd-6311gss-12A-ACE2_Ab-Wild-type_"$x"_$i
			for k in `cat file.txt`
			do
                            echo "==>$k" >> checkerr.txt
			    tail -n 4 $k.log >> checkerr.txt
			done
			    grep "/soft/apps/Paid/g09/l801.exe" checkerr.txt > exist801.txt
			if [ -s exist801.txt ]
			then
			    grep -B2 "l801.exe" checkerr.txt > sub801.txt
                            grep "==>" sub801.txt > sub801name.txt
                            sed -i "s/==>//g" sub801name.txt
                            head -n 4 sub801name.txt
			    for j in `cat sub801name.txt`
			    do
				sed -i "s/scf(maxcyc=100,qc)/scf(maxcyc=100,qc) iop(8\/11=1)/g" $j.gjf
				sed -i "s/scf(maxcyc=100)/scf(maxcyc=100) iop(8\/11=1)/g" $j.gjf
				head -n 4 $j.gjf
				bsub -n 24 -q c_soft -R "select[hname!="m4510"]" -e err.log -o out.log g09 $j.gjf
			    done
			else
			    echo "*****************There are no l801 tasks.********************"
			fi
		else
			echo "!!!!!!!!!!!!!!!!!!!!!!$FILE not exist.!!!!!!!!!!!!!!!!!!!!!!!"
		fi
	done
done





