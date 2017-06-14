      program make_spectra
      use readio  
      use spectra       
      use readio_medium
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
      end program make_spectra
