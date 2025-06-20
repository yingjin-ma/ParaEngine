#!/bin/bash
rm /work1/scquant/ParaEngine-git-COVID2019/Workdir-MP2-631gss-12A-ACE2_Ab-Omicron_Q493R*/checkwrite.txt
rm /work1/scquant/ParaEngine-git-COVID2019/Workdir-MP2-631gss-12A-ACE2_Ab-Omicron_Q493R*/write.txt 
rm /work1/scquant/ParaEngine-git-COVID2019/Workdir-MP2-631gss-12A-ACE2_Ab-Omicron_Q493R*/subwrite.txt
for x in {0..9}
do
	for i in {0..9}
		do
		echo  "----------------Workdir-MP2-631gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i-----------------------"
		FILE=/work1/scquant/ParaEngine-git-COVID2019/Workdir-MP2-631gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i
		if [ -d "$FILE" ]
		then
			echo "====================/work1/scquant/ParaEngine-git-COVID2019/Workdir-MP2-631gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i========================="
			cd /work1/scquant/ParaEngine-git-COVID2019/Workdir-MP2-631gss-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i
			for k in `cat file.txt`
			do
                            echo "==>$k" >> checkwrite.txt
			    tail -n 1 $k.log >> checkwrite.txt
			done
			    grep "g_write" checkwrite.txt >> write.txt
			    grep "Error" checkwrite.txt >>write.txt
			    grep "JobTyp=" checkwrite.txt >>write.txt
			if [ -s write.txt ]
			then
			    grep -B1 "g_write" checkwrite.txt >> subwrite.txt
			    grep -B1 "Error" checkwrite.txt >>subwrite.txt
			    grep -B1 "JobTyp=" checkwrite.txt >>subwrite.txt
			    grep "==>" subwrite.txt > writename.txt
			    sed -i "s/==>//g" writename.txt
			    head -n 4 writename.txt
			    for j in `cat writename.txt`
			    do
			        sed -i "s/B2PLYP\/6-31G\*\*/B2PLYP=FullDirect\/6-31G\*\*/g" $j.gjf
			            head -n 4 $j.gjf
				    echo "=========================The $j.gjf is submitted.============================"
			            bsub -n 24 -q c_soft -e err.log -o out.log -R "select[hname!="m5610" hname!="m4516" hname!="c5433"]" g09 $j.gjf
			    done
			else
			    echo "*****************There are no g_write or Error tasks.********************"
			fi
		else
			echo "!!!!!!!!!!!!!!!!!!!!!!$FILE not exist.!!!!!!!!!!!!!!!!!!!!!!!"
		fi
	done
done





