!#/bin/bash

for i in 0
do
	for x in {2..9}
	do
		cp input-wb97xd-6311gss-Wild-type_0_0-STDDLB-10nodes input-wb97xd-6311gss-Wild-type_"$i"_$x-STDDLB-10nodes
	        sed -i "s/type_0_0/type_"$i"_$x/g" input-wb97xd-6311gss-Wild-type_"$i"_$x-STDDLB-10nodes
		echo "input-wb97xd-6311gss-Wild-type_"$i"_$x-STDDLB-10nodes"
		sed -n '4,9p' input-wb97xd-6311gss-Wild-type_"$i"_$x-STDDLB-10nodes
	done
done

