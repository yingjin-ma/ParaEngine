!#/bin/bash
for x in {0..5}
  do
  for i in {0..9}
    do
    echo "==================================>/work1/scquant/ParaEngine-git-COVID2019/Workdir-M062X-631gss-GD3-12A-ACE2_Ab-Wild-type_"$x"_$i"
    cd /work1/scquant/ParaEngine-git-COVID2019/Workdir-M062X-631gss-GD3-12A-ACE2_Ab-Wild-type_"$x"_$i
    tail -n 1 */*.log > tail.txt
    sed -e "/Normal/d" -e "/^\s*$/d" tail.txt > tail1.txt
    grep -B1 -v "==>" tail1.txt > err.txt
    grep "==>" err.txt > err2.txt
    sed -e "s/\s//g" -e "s/==>//g" -e "s/.log<==//g" err2.txt  > file.txt

    for j in `cat file.txt`
	do
		tail -n 4 $j.log > checktail.txt
	echo "-------------------------------------------------------------------------------------------------------------"
	grep "Output" checktail.txt
	if [ -s file.txt ]
		then
        		echo "file has content"
		else
        		echo "file is empty"
	fi 
	 
	done
    done
  done
