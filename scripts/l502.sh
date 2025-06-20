#!/bin/bash
for x in {0..9}
do
	for i in {0..9}
		do
		echo  "----------------Workdir-M062X-631gss-GD3-12A-ACE2_Ab-Wild-type_"$x"_$i-----------------------"
		FILE=/work1/scquant/ParaEngine-git-COVID2019/Workdir-M062X-631gss-GD3-12A-ACE2_Ab-Wild-type_"$x"_$i
		if [ -d "$FILE" ]
		then
    			#echo "*******************$FILE is a directory.*****************"
			cd /work1/scquant/ParaEngine-git-COVID2019/Workdir-M062X-631gss-GD3-12A-ACE2_Ab-Wild-type_"$x"_$i
			for k in `cat file.txt`
			do
			    tail -n 4 $k.log > check502.txt
			    grep "/soft/apps/Paid/g09/l502.exe" check502.txt > exist502.txt
			done
			if [ -s exist502.txt ]
			then
			    grep -B2 "/soft/apps/Paid/g09/l502.exe" check502.txt > sub502.txt
			    grep "==>" sub502.txt > sub502name.txt
			    sed -i "s/==>//g" sub502name.txt
			    for j in `cat sub502name.txt`
			    do
			        sed -i "s/scf(maxcyc=100)/scf(maxcyc=100,qc)/g" $j.gjf
			       # grep "scf(maxcyc=100,qc)" $j.gjf > gjf.txt
			       # if [ -s gjf.txt ]
			       # then
			            echo "adding scf=qc is ok"
			           # bsub -n 24 -q cpuII -e err.log -o out.log -R "select[hname!="m5610" hname!="m4516" hname!="c5433"]" g09 $j.gjf
			       # fi
			    done
			else
			    echo "*****************There's no l502 task.********************"
			fi
		else
			echo "!!!!!!!!!!!!!!!!!!!!!!$FILE not exist.!!!!!!!!!!!!!!!!!!!!!!!"
		fi
	done
done





