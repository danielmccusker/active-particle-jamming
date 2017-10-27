#!/bin/bash

rho=( 0.84 0.845 0.86 0.9 )
l_s=( 0.003 0.01 0.03 0.1 0.3 )
l_n=( 0.01 0.1 0.3 0.4 0.5 0.6 0.7 1.0 )

noCells=1024
noSteps=100000000
stepsPerTime=10
run=0
ID="newCode1024"

rm -f input.txt
g++ -pg active_jam_nbr_15.cpp -I boost_1_64_0/ -O3 -o a.out -std=c++11
for i in ${rho[@]}
do
	for j in ${l_n[@]}
	do
		for k in ${l_s[@]}
		do
                        printf "$ID run$run $noCells $noSteps $stepsPerTime $k $j $i\n" >> input.txt
			run=$((run+1))
		done
	done
done
exit

