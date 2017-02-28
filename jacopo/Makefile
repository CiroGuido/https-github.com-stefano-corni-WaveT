#FC = gfortran
FC = ifort    
#FFLAGS_XEON= -O3 -ftz -align -arch pn4 -tune pn4 -xN -tpp7
#FFLAGS_PAR= -parallel -par_threshold75 -par_report3
#FFLAGS_DEB= -debug -C -g
#FFLAGS_DEB= -debug -traceback -C -g -diag-disable  warn -parallel -par_threshold75 -par_report3
#FFLAGS_DEB_PAR= -g -diag-disable  warn -parallel -par_threshold75 -par_report3
#FFLAGS_P3= -O3 -ftz -align -arch pn3 -tune pn3 -tpp6
#FFLAGS_P4= -O3 -ftz -align -arch pn4 -tune pn4 -tpp7 -xN
#FFLAGS_MIO= -O3 -ftz -align -fast   
#FFLAGS_HYD= -O3 -ftz -align -xHost 
#FFLAGS_HYD_PAR= -O3 -ftz -align -xHost -parallel -par_threshold75 -par_report3
#FFLAGS_WS= -O3 -ftz -align -xP -no-prec-div -ipo
#FFLAGS_GF= -O3
#LDLAGS= 
##FFLAGS= $(FFLAGS_HYD_PAR)
#LIBS = -L/usr/lib -llapack  -lblas 
LIBS = -lmkl_sequential -lmkl_core -lpthread -lmkl_intel_lp64
#LIBSF = -lm -L/usr/lib/i386-linux-gnu -lfftw3 -c 
INC = -I/usr/include/ -I/unimore/prod/fftw-3.3.4/include/

OBJ= random.o cav_types.o pedra_friends.o readio.o readio_medium.o rhomatrix.o spectra.o BEM_medium.o scf.o td_contmed.o QM_coupling.o propagate.o main.o

embem.x: $(OBJ) 
	$(FC) $(FFLAGS) $(LDLAGS) -o $@ $(OBJ)  $(LIBS) $(LIBSF)

clean:
	rm -f *.o *.mod *.il

random.o: random.f90
	$(FC) -c $(FFLAGS) $(INC)  $<

cav_types.o: cav_types.f90
	$(FC) -c $(FFLAGS) $(INC) $<

pedra_friends.o: pedra_friends.f90
	$(FC) -c $(FFLAGS) $(INC)  $<

readio.o: readio.f90
	$(FC) -c $(FFLAGS) $(INC) $<

readio_medium.o: readio_medium.f90
	$(FC) -c $(FFLAGS) $(INC) $<

scf.o: scf.f90
	$(FC) -c $(FFLAGS) $(INC) $<

rhomatrix.o: rhomatrix.f90
	$(FC) -c $(FFLAGS) $(INC)  $<

propagate.o: propagate.f90
	$(FC) -c $(FFLAGS) $(INC) $<

spectra.o: spectra.f90
	$(FC) -c $(FFLAGS)  $(INC)  $<

BEM_medium.o: BEM_medium.f90
	$(FC) -c $(FFLAGS) $(INC) $<

td_contmed.o: td_contmed.f90
	$(FC) -c $(FFLAGS) $(INC) $<

QM_coupling.o: QM_coupling.f90
	$(FC) -c $(FFLAGS) $(INC)  $<

main.o: main.f90
	$(FC) -c $(FFLAGS) $(INC) $<

