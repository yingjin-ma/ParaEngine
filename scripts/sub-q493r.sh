!#/bin/bash
for x in {0..9}
  do
  for i in {0..9}
    do
    echo "/work1/scquant/ParaEngine-git-COVID2019/Workdir-M062X-631gss-GD3-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i"
    cd /work1/scquant/ParaEngine-git-COVID2019/Workdir-M062X-631gss-GD3-12A-ACE2_Ab-Omicron_Q493R_"$x"_$i
    head -n 3 err.txt
    for j in `cat file.txt`
	tail -n 4 $j.log
	done
    done
  done
