#!/bin/bash

rho=( 0.84 0.845 )
l_s=( 0.2 )
l_n=( 0.3 )

noSteps=100000
stepsPerTime=10
noCells=10000
run=0
ID="10000"

rm -f input-prof.txt
g++ -pg active_jam_nbr_15.cpp -I boost_1_64_0/ -O3 -o a.out -std=c++11

for i in ${rho[@]}
do
	for j in ${l_n[@]}
	do
		for k in ${l_s[@]}
		do
			printf "$ID run$run $noCells $noSteps $stepsPerTime $k $j $i\n" >> input-prof.txt
			run=$((run+1))
		done
	done
done
exit

