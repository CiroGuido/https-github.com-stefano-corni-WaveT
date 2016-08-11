#! /bin/sh

#PBS -N gaussian
#PBS -j oe
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
#PBS -q s3par8c

export GAUSS_SCRDIR=/s3_home/spipolo/GAUSS_SCRDIR
export g09root=/scratch/spipolo/g09_parma/g09A02
. /unimore/prod/g09/bsd/g09.profile

cd  $PBS_O_WORKDIR

g09 < LiCN.gjf > LiCN.log   


