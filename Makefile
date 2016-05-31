FC = gfortran
FFLAGS_XEON= -O3 -ftz -align -arch pn4 -tune pn4 -xN -tpp7
FFLAGS_PAR= -parallel -par_threshold75 -par_report3
FFLAGS_DEB= -debug -C -g
FFLAGS_P3= -O3 -ftz -align -arch pn3 -tune pn3 -tpp6
FFLAGS_P4= -O3 -ftz -align -arch pn4 -tune pn4 -tpp7 -xN
FFLAGS_MIO= -O3 -ftz -align -arch pn4 -tune pn4 -tpp7 -xB -ipo
FFLAGS_WS= -O3 -ftz -align -xP -no-prec-div -ipo
FFLAGS_GF= -O3
FFLAGS= $(FFLAGS_DEB)
LIBS = -L/usr/lib -llapack  -lblas 
LIBSF = -lm -L/usr/lib/i386-linux-gnu -lfftw3
INC = -I/usr/include/ 

OBJ= cav_types.o pedra_friends.o readio.o readio_medium.o spectra.o BEM_medium.o scf.o td_contmed.o QM_coupling.o propagate.o main.o 

embem.x: $(OBJ) 
	$(FC) $(FFLAGS) -o $@ $(OBJ)  $(LIBS) $(LIBSF)

clean:
	rm -f *.o *.mod *.il

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

