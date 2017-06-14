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

OBJ_TDPLAS= global_tdplas.o pedra_friends.o readio_medium.o BEM_medium.o main_tdplas.o
OBJ_FREQ= global_tdplas.o pedra_friends.o readio_medium.o BEM_medium.o scf.o td_contmed.o main_freq.o 

# rules

all: $(BINPATH)/freq_bem.x $(BINPATH)/tdplas.x

$(BINPATH)/tdplas.x: $(OBJ_TDPLAS) 
	$(FC) $(FFLAGS) -o $@ $(OBJ_TDPLAS) $(LIBS) $(LIBSF)

$(BINPATH)/freq_bem.x: $(OBJ_FREQ) 
	$(FC) $(FFLAGS) -o $@ $(OBJ_FREQ) $(LIBS) $(LIBSF)

clean:
	rm -f *.o *.mod *.il $(BINPATH)/*.x

global_tdplas.o: global_tdplas.f90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

pedra_friends.o: pedra_friends.f90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

readio_medium.o: readio_medium.f90
	$(FC) $(FFLAGS) -c $< -o $@

BEM_medium.o: BEM_medium.f90
	$(FC) $(FFLAGS) -c $< -o $@

QM_coupling.o: QM_coupling.f90
	$(FC) $(FFLAGS) -c $< -o $@

scf.o: scf.f90
	$(FC) $(FFLAGS) -c $< -o $@

td_contmed.o: td_contmed.f90
	$(FC) $(FFLAGS) -c $< -o $@

main_tdplas.o: main_tdplas.f90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

main_freq.o: main_freq.f90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

# dependencies

global_tdplas.o: 

pedra_friends.o: global_tdplas.o

readio_medium.o: global_tdplas.o pedra_friends.o

BEM_medium.o: global_tdplas.o pedra_friends.o readio_medium.o

QM_coupling.o: global_tdplas.o pedra_friends.o readio_medium.o BEM_medium.o

scf.o: global_tdplas.o pedra_friends.o readio_medium.o BEM_medium.o

td_contmed.o: global_tdplas.o pedra_friends.o readio_medium.o BEM_medium.o scf.o

main_tdplas.o: readio_medium.o BEM_medium.o

main_freq.o: readio_medium.o td_contmed.o


