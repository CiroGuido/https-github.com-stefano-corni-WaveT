#!/bin/csh
#
setenv fileinp tdcis_vac
foreach SEED (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)
cat <<EOF>  ${fileinp}_${SEED}.inp
$SEED ! numero aggiunto ai nomi dei file .dat, manca ancora medium_t 
10    ! numero stati eccitati ci
0.2  2000   1  ! dt(au) e numero step
dis  0.1  !dissipazione (diss), tipo di distribuzione (gaussiana per ora) coefficiente di attrito per random phase walk
pip
300 5 0.353259 ! centro dell'impulso (in a.u.), larghezza non usata (a.u.), frequenza angolare
0.000  0.000  1.000 ! valore massimo del campo (vettore)
non !rad aggiunge lo smorzamento radiativo
0.     0.     0.     ! molecular center of charge
2.0                  ! tau spectrum
0.                   ! start FT
1   ! direction for FT
vac!input for medium
1  ! stride between updating the interaction potential
ons ! coupling BEM (pcm) or dipolar (ons)
5.0 5.0 5.0
! input for dielectric function
deb ! Debye vs Drude-Lorentz (drl)
2.0 1.02 1000.  1. ! eps_0, eps_d, tau_deb(in au)  
! propagation scheme:dip evolve onsager field....
dip
! other options
NON     ! scf 
0       ! Include local Field ?  Binary(1=yes,0=no)
0       ! Debug ?  Binary(1=yes,0=no)
0       ! Equilibrium calculation? Binary(1=yes,0=no)
0       ! Equilibrium calculation? Binary(1=yes,0=no)

EOF

/scratch/corni/jacopo/WaveT-master/embem.x ${fileinp}_${SEED}.inp








