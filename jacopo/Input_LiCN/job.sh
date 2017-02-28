#!/bin/sh
#PBS -N jacopo      
#PBS -j oe
#PBS -q s3par6c
#PBS -l nodes=1:ppn=12
#PBS -l walltime=72:00:00
##PBS -W  group_list=moe-sc
cd  $PBS_O_WORKDIR

JACOPOWT=/s3_home/jfregoni/WaveT-master/Input_LiCN

cd $JACOPOWT

module load composer_xe_2013.4.183
#./scan_seed.x

/s3_home/jfregoni/WaveT-master/embem.x < tdcis_vac.inp > primaprova.out
