#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=5:00:00
#$ -l h_vmem=4G # memory request
#$ -t 1-3

echo "Initializing assembly - time series generated after each invasion"

PARFILE=$HOME/LVMCM_src/LVMCM/parFiles/autonomous_turnover_example/autonomous_turnover_example_pars.txt
cd $HOME/LVMCM_src/LVMCM/build

./LVMCM $(sed -n "${SGE_TASK_ID}p" $PARFILE)

echo "Assembly complete"

exit 0
