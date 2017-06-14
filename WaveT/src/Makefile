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
#LIBS =  -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -Wl,--end-group#
#LIBSF = -lm -L/usr/lib/i386-linux-gnu -lfftw3
LIBSF = -lm -L/usr/lib/x86_64-linux-gnu -lfftw3
INC = -I/usr/include/ 
#INC = -I/usr/include/ -I/unimore/prod/fftw-3.3.4/include/

BINPATH=../bin
TDPLAS=../../TDPlas/src/

OBJ= random.o readio.o spectra.o dissipation.o propagate.o main.o global_tdplas.o pedra_friends.o readio_medium.o BEM_medium.o QM_coupling.o scf.o td_contmed.o  
OBJ_SPECTRA= readio.o spectra.o main_spectra.o global_tdplas.o pedra_friends.o readio_medium.o  

# rules

all: $(BINPATH)/WaveT.x $(BINPATH)/make_spectra.x

$(BINPATH)/make_spectra.x: $(OBJ_SPECTRA) 
	$(FC) $(FFLAGS) -o $@ $(OBJ_SPECTRA) $(LIBS) $(LIBSF)

$(BINPATH)/WaveT.x: $(OBJ) 
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LIBS) $(LIBSF)

clean:
	rm -f *.o *.mod *.il $(BINPATH)/*.x

random.o: random.f90
	$(FC) $(FFLAGS) -c $< -o $@

readio.o: readio.f90
	$(FC) $(FFLAGS) -c $< -o $@

spectra.o: spectra.f90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

main_spectra.o: main_spectra.f90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

dissipation.o: dissipation.f90
	$(FC) $(FFLAGS) -c $< -o $@

propagate.o: propagate.f90
	$(FC) $(FFLAGS) -c $< -o $@

main.o: main.f90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@


global_tdplas.o: $(TDPLAS)/global_tdplas.f90
	$(FC) $(FFLAGS) -c $< -o $@

pedra_friends.o: $(TDPLAS)/pedra_friends.f90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

readio_medium.o: $(TDPLAS)/readio_medium.f90
	$(FC) $(FFLAGS) -c $< -o $@

BEM_medium.o: $(TDPLAS)/BEM_medium.f90
	$(FC) $(FFLAGS) -c $< -o $@

QM_coupling.o: $(TDPLAS)/QM_coupling.f90
	$(FC) $(FFLAGS) -c $< -o $@

scf.o: $(TDPLAS)/scf.f90
	$(FC) $(FFLAGS) -c $< -o $@

td_contmed.o: $(TDPLAS)/td_contmed.f90
	$(FC) $(FFLAGS) -c $< -o $@

# dependencies

random.o: 

readio.o: global_tdplas.o

spectra.o: readio.o pedra_friends.o readio_medium.o

main_spectra.o: readio.o readio_medium.o spectra.o

dissipation.o: readio.o random.o

propagate.o: random.o readio.o spectra.o readio_medium.o pedra_friends.o td_contmed.o

main.o: readio.o spectra.o dissipation.o propagate.o global_tdplas.o readio_medium.o QM_coupling.o  

global_tdplas.o: 

pedra_friends.o: global_tdplas.o

readio_medium.o: global_tdplas.o pedra_friends.o

BEM_medium.o: global_tdplas.o pedra_friends.o readio_medium.o

QM_coupling.o: global_tdplas.o pedra_friends.o readio_medium.o BEM_medium.o

scf.o: global_tdplas.o pedra_friends.o readio_medium.o BEM_medium.o

td_contmed.o: global_tdplas.o pedra_friends.o readio_medium.o BEM_medium.o scf.o
