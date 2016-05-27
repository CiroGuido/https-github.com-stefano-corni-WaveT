      program tdcis
      use readio  
      use readio_medium
      use propagate    
      use QM_coupling  
      use spectra       
      implicit none
!
!     read in the input parameter for the present evolution
      call read_input
      if (mdm.ne."vac") call read_medium
      call init_spectra
!
!     create the field 
      call create_field
!
!     propagate or diagonalise matrix
      if(nmodes.eq.0) then
        call prop
        call do_spectra
      else 
        call do_diag    
      endif
!         
      stop
      end
