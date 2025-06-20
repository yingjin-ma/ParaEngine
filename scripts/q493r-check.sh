#!/bin/bash
rm /work1/scquant/ParaEngine-git-COVID2019/Workdir-B2PLYP-631gss-GD3-12A-ACE2_Ab-Omicron_Q493R*/checkwrite.txt
for x in {0..9}
do
	for i in {0..9}
		do
		echo  "----------------Workdir-B2PLYP-631gss-GD3-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i-----------------------"
		FILE=/work1/scquant/ParaEngine-git-COVID2019/Workdir-B2PLYP-631gss-GD3-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i
		if [ -d "$FILE" ]
		then
    			#echo "*******************$FILE is a directory.*****************"
			cd /work1/scquant/ParaEngine-git-COVID2019/Workdir-B2PLYP-631gss-GD3-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i
			for k in `cat file.txt`
			do
                            echo "==>$k" >> checkwrite.txt
			    tail -n 4 $k.log >> checkwrite.txt
			    cat checkwrite.txt
			done
		else
			echo "!!!!!!!!!!!!!!!!!!!!!!$FILE not exist.!!!!!!!!!!!!!!!!!!!!!!!"
		fi
	done
done





