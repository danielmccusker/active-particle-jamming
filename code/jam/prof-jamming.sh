#!/bin/bash
#PBS -l nodes=1:ppn=1

parameters=$(sed -n -e "${PBS_ARRAYID}p" /home/dmccusker/remote/jamming-dynamics/code/jam/input-prof.txt)
parameterArray=($parameters)

ID=${parameterArray[0]}
currentRun=${parameterArray[1]}
noCells=${parameterArray[2]}
noSteps=${parameterArray[3]}
stepsPerTime=${parameterArray[4]}
lambda_s=${parameterArray[5]}
lambda_n=${parameterArray[6]}
rho=${parameterArray[7]}

/home/dmccusker/remote/jamming-dynamics/code/jam/a-prof.out $ID $currentRun $noCells $noSteps $stepsPerTime $lambda_s $lambda_n $rho
