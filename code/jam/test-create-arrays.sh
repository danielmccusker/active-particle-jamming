#!/bin/bash

rho=( 0.84 0.845 )
l_s=( 0.003 )
l_n=( 0.01 )

noSteps=10000000
stepsPerTime=10
run=0
ID="test"

g++ test.cpp -I boost_1_64_0/ -O3 -o test.out -std=c++11
rm -f input.txt
#rm ../../output/test -rf
#rm -f jamming.sh.e* jamming.sh.o*

for i in ${rho[@]}
do
	for j in ${l_n[@]}
	do
		for k in ${l_s[@]}
		do
			printf "$ID run$run $noSteps $stepsPerTime $k $j $i\n" >> input.txt
			run=$((run+1))
		done
	done
done
exit

