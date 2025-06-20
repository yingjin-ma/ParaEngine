#!/bin/bash
rm /work1/scquant/ParaEngine-git-COVID2019/Workdir-wb97xd-6311gss-12A-ACE2_Ab-Wild-type*/check502.txt
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
                            echo "==>$k" >> check502.txt
			    tail -n 4 $k.log >> check502.txt
			done
			    grep "/soft/apps/Paid/g09/l502.exe" check502.txt > exist502.txt
			
			if [ -s exist502.txt ]
			then
                            
			    grep -B2 "/soft/apps/Paid/g09/l502.exe" check502.txt > sub502.txt
			    grep "==>" sub502.txt > sub502name.txt
			    sed -i "s/==>//g" sub502name.txt
			    head -n 4 sub502name.txt
                            #exit 0
			    for j in `cat sub502name.txt`
			    do
			        sed -i "s/scf(maxcyc=100)/scf(maxcyc=100,qc)/g" $j.gjf
			        #grep "scf(maxcyc=100,qc)" $j.gjf > gjf.txt
			        #if [ -s gjf.txt ]
			        #then
			            #echo "adding scf=qc is ok"
			            bsub -n 24 -q cpuII -R "select[hname!='m4510']" -e err.log -o out.log  g09 $j.gjf
			        #fi
			    done
                            #exit 0
			else
			    echo "*****************There's no l502 task.********************"
			fi
		else
			echo "!!!!!!!!!!!!!!!!!!!!!!$FILE not exist.!!!!!!!!!!!!!!!!!!!!!!!"
		fi
	done
done





