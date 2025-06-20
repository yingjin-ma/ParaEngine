#!/bin/bash
rm /work1/scquant/ParaEngine-git-COVID2019/Workdir-MP2-631gss-12A-ACE2_Ab-Wild-type*/checkwrite.txt
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
                            echo "==>$k" >> checkwrite.txt
			    tail -n 1 $k.log >> checkwrite.txt
			done
			    grep "g_write" checkwrite.txt > write.txt
			    grep "Error" checkwrite.txt >>write.txt
			    grep "JobTyp" checkwrite.txt >>write.txt
			
			if [ -s write.txt ]
			then
                            
			    grep -B1 "g_write" checkwrite.txt > subwrite.txt
			    grep -B1 "Error" checkwrite.txt >>subwrite.txt
			    grep -B1 "JobTyp" checkwrite.txt >>subwrite.txt
			    grep "==>" subwrite.txt > writename.txt
			    sed -i "s/==>//g" writename.txt
			    head -n 4 writename.txt
                            #exit 0
			    for j in `cat writename.txt`
			    do
			        sed -i "s/B2PLYP\/6-31G\*\*/B2PLYP=FullDirect\/6-31G\*\*/g" $j.gjf
			        #grep "scf(maxcyc=100,qc)" $j.gjf > gjf.txt
			        #if [ -s gjf.txt ]
			        #then
			            head -n 4 $j.gjf
				    echo "=========================The $j.gjf is submitted.============================"
			            bsub -n 24 -q c_soft -e err.log -o out.log -R "select[hname!="c5433"]" g09 $j.gjf
			        #fi
			    done
                            #exit 0
			else
			    echo "*****************There's no g_write task.********************"
			fi
		else
			echo "!!!!!!!!!!!!!!!!!!!!!!$FILE not exist.!!!!!!!!!!!!!!!!!!!!!!!"
		fi
	done
done





