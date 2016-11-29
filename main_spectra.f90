      program make_spectra
      use readio  
      use readio_medium
      use spectra       
      implicit none
!
!     read in the input parameter for the present evolution
      call read_input
      if (mdm.ne."vac") call read_medium
      call init_spectra
!
!     calculate spectra               
      call read_arrays 
      call do_spectra
!         
      stop
      end
