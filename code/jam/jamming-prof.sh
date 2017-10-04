#!/bin/bash
#PBS -l nodes=1:ppn=2

parameters=$(sed -n -e "${PBS_ARRAYID}p" /home/dmccusker/remote/jamming_cluster/code/jam/input-prof.txt)
parameterArray=($parameters)

ID=${parameterArray[0]}
currentRun=${parameterArray[1]}
noSteps=${parameterArray[2]}
stepsPerTime=${parameterArray[3]}
lambda_s=${parameterArray[4]}
lambda_n=${parameterArray[5]}
rho=${parameterArray[6]}

/home/dmccusker/remote/jamming_cluster/code/jam/a-prof.out $ID $currentRun $noSteps $stepsPerTime $lambda_s $lambda_n $rho
