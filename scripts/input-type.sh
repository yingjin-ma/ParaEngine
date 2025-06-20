#!/bin/bash
for x in {1..9}
do
    for i in {0..9}
    do
	cp input-b2plyp-631gss-gd3-Wild-type_0_0-STDDLB-10nodes input-b2plyp-631gss-gd3-Wild-type_"$x"_$i-STDDLB-10nodes 
	sed -i "s/Wild-type_0_0/Wild-type_"$x"_$i/g" input-b2plyp-631gss-gd3-Wild-type_"$x"_$i-STDDLB-10nodes
	sed -n "4,9p" input-b2plyp-631gss-gd3-Wild-type_"$x"_$i-STDDLB-10nodes
    done
done
