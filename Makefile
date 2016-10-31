#FC = gfortran
FC = ifort
FFLAGS_XEON= -O3 -ftz -align -arch pn4 -tune pn4 -xN -tpp7
FFLAGS_PAR= -parallel -par_threshold75 -par_report3
FFLAGS_DEB= -debug -C -g
FFLAGS_P3= -O3 -ftz -align -arch pn3 -tune pn3 -tpp6
FFLAGS_P4= -O3 -ftz -align -arch pn4 -tune pn4 -tpp7 -xN
FFLAGS_MIO= -O3 -ftz -align -arch pn4 -tune pn4 -tpp7 -xB -ipo
FFLAGS_WS= -O3 -ftz -align -xP -no-prec-div -ipo
FFLAGS_LAP= -O3 -ftz -align -xAVX -no-prec-div 
FFLAGS_GF= -O3
LDLAGS= 
#LDLAGS= --static --disable-shared --enable-static
FFLAGS= $(FFLAGS_LAP)
#LIBS = -L/usr/lib -llapack  -lblas 
#LIBS =  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread
LIBS =  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
#LIBS =  -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -Wl,--end-group
#LIBSF = -lm -L/usr/lib/i386-linux-gnu -lfftw3
INC = -I/usr/include/ 

OBJ= random.o cav_types.o pedra_friends.o readio.o readio_medium.o spectra.o BEM_medium.o scf.o td_contmed.o QM_coupling.o propagate.o main.o 
OBJ_FREQ= cav_types.o pedra_friends.o readio.o readio_medium.o spectra.o BEM_medium.o scf.o td_contmed.o main_freq.o 

embem.x: $(OBJ) 
	$(FC) $(FFLAGS) -o $@ $(OBJ)  $(LIBS) $(LIBSF)

freq_bem.x: $(OBJ_FREQ) 
	$(FC) $(FFLAGS) -o $@ $(OBJ_FREQ)  $(LIBS) $(LIBSF)

clean:
	rm -f *.o *.mod *.il

random.o: random.f90
	$(FC) -c $(FFLAGS) $<

cav_types.o: cav_types.f90
	$(FC) -c $(FFLAGS) $(INC) $<

pedra_friends.o: pedra_friends.f90
	$(FC) -c $(FFLAGS) $(INC) $<

readio.o: readio.f90
	$(FC) -c $(FFLAGS)  $<

readio_medium.o: readio_medium.f90
	$(FC) -c $(FFLAGS)  $<

scf.o: scf.f90
	$(FC) -c $(FFLAGS)  $<

propagate.o: propagate.f90
	$(FC) -c $(FFLAGS)  $<

spectra.o: spectra.f90
	$(FC) -c $(FFLAGS)  $<

BEM_medium.o: BEM_medium.f90
	$(FC) -c $(FFLAGS)  $<

td_contmed.o: td_contmed.f90
	$(FC) -c $(FFLAGS)  $<

QM_coupling.o: QM_coupling.f90
	$(FC) -c $(FFLAGS)  $<

main.o: main.f90
	$(FC) -c $(FFLAGS) $(INC) $<

main_freq.o: main_freq.f90
	$(FC) -c $(FFLAGS) $(INC) $<

