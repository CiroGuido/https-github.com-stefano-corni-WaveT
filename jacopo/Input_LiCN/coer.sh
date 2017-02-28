#!/bin/csh
#
set rhocol=`awk '{print $2,($7)**2+($8)**2}' rho_1.dat >> risu.out`
