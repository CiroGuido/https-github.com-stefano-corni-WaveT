#!/bin/bash

file=$1
sts=$2
let ste=sts
let sts=sts+1
grep "States  " $file | awk 'NR<='$sts'*('$sts'+1)/2' > ci_mut.inp
awk '/Root      1 :/{for(i=1;i<='$sts';i++){print $0; getline}}' $file | head -n $ste > ci_energy.inp

