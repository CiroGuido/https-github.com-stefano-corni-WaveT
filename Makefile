#FC = gfortran
FC = ifort
FFLAGS_XEON= -O3 -ftz -align -arch pn4 -tune pn4 -xN -tpp7
FFLAGS_PAR= -parallel -par_threshold75 -par_report3
FFLAGS_HYD= -O3 -ftz -align -xHost
FFLAGS_HYD_PAR= -O3 -ftz -align -xHost -parallel -par_threshold75 -par_report3
FFLAGS_DEB= -debug -C -traceback
FFLAGS_P3= -O3 -ftz -align -arch pn3 -tune pn3 -tpp6
FFLAGS_P4= -O3 -ftz -align -arch pn4 -tune pn4 -tpp7 -xN
FFLAGS_MIO= -O3 -ftz -align -arch pn4 -tune pn4 -tpp7 -xB -ipo
FFLAGS_WS= -O3 -ftz -align -xP -no-prec-div -ipo
FFLAGS_LAP= -O3 -ftz -align -xAVX -no-prec-div 
FFLAGS_GF= -O3
LDLAGS= 
#LDLAGS= --static --disable-shared --enable-static
FFLAGS= $(FFLAGS_DEB)
#LIBS = -L/usr/lib -llapack  -lblas 
#LIBS =  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread
LIBS =  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
#LIBS =  -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -Wl,--end-group
#LIBSF = -lm -L/usr/lib/i386-linux-gnu -lfftw3
INC = -I/usr/include/ 
#INC = -I/usr/include/ -I/unimore/prod/fftw-3.3.4/include/

OBJ= random.o cav_types.o pedra_friends.o readio.o readio_medium.o spectra.o BEM_medium.o scf.o td_contmed.o QM_coupling.o dissipation.o propagate.o main.o 
OBJ_FREQ= cav_types.o pedra_friends.o readio.o readio_medium.o spectra.o BEM_medium.o scf.o td_contmed.o main_freq.o 
OBJ_TDPLAS= cav_types.o pedra_friends.o readio.o readio_medium.o  BEM_medium.o main_tdplas.o 
OBJ_SPECTRA= cav_types.o pedra_friends.o readio.o readio_medium.o spectra.o main_spectra.o 

all: WaveT.x freq_bem.x make_spectra.x tdplas.x

WaveT.x: $(OBJ) 
	$(FC) $(FFLAGS) -o $@ $(OBJ)  $(LIBS) $(LIBSF)

freq_bem.x: $(OBJ_FREQ) 
	$(FC) $(FFLAGS) -o $@ $(OBJ_FREQ)  $(LIBS) $(LIBSF)

make_spectra.x: $(OBJ_SPECTRA) 
	$(FC) $(FFLAGS) -o $@ $(OBJ_SPECTRA)  $(LIBS) $(LIBSF)

tdplas.x: $(OBJ_TDPLAS) 
	$(FC) $(FFLAGS) -o $@ $(OBJ_TDPLAS)  $(LIBS) $(LIBSF)

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

dissipation.o: dissipation.f90
	$(FC) -c $(FFLAGS)  $<

propagate.o: propagate.f90
	$(FC) -c $(FFLAGS)  $<

spectra.o: spectra.f90
	$(FC) -c $(FFLAGS)  $(INC)  $<

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

main_spectra.o: main_spectra.f90
	$(FC) -c $(FFLAGS) $(INC) $<

main_tdplas.o: main_tdplas.f90
	$(FC) -c $(FFLAGS) $(INC) $<

