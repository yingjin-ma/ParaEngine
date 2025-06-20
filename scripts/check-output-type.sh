#!/bin/bash
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
			grep "Out" err.txt > checkout.txt
			grep "write" err.txt >>checkout.txt
			grep "NtrErr" err.txt >>checkout.txt
			grep "Error" err.txt >>checkout.txt
			if [ -s checkout.txt ]
                	then
                           grep -B1 "Out" err.txt > newout.txt
		           grep -B1 "write" err.txt >> newout.txt
			   grep -B1 "Error" err.txt >> newout.txt
			   grep -B1 "NtrErr" err.txt >> newout.txt 
        	           sed -e "/Out/d" -e "/--/d" -e "s/[[:space:]]//g" -e "s/==>//g" -e "s/.log<==//g" -e "/g_write/d" -e "/NtrErr/d" -e "/Error/d" newout.txt > out1.txt
		           for j in `cat out1.txt`
		       	  	   do
		                   bsub -n 24 -q cpuII -e err.log -o out.log -R "select[hname!="m4510"]" g09 $j.gjf
		                   done
		           head -n 3 out1.txt
			fi
		else
			echo "!!!!!!!!!!!!!!!!!!!!!!$FILE not exist.!!!!!!!!!!!!!!!!!!!!!!!"
		fi
	        done
done	
