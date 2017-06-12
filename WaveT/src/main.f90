      program tdcis
      use readio  
      use readio_medium
      use propagate    
      use QM_coupling  
      use spectra
      use dissipation       
      implicit none
      integer :: st,current,rate
!
!     read in the input parameter for the present evolution
      call system_clock(st,rate)
      call read_input
      if (mdm.ne."vac") call read_medium
      call init_spectra
!
!     create the field 
      call create_field
      call system_clock(current)
      write(6,'("Done reading input & setting up the field, took", &
            F10.3,"s")') real(current-st)/real(rate)
!
!     propagate or diagonalise matrix
      if(nmodes.eq.0) then
        call prop
        call do_spectra
      else 
        call do_diag    
      endif
      call system_clock(current)
      write(6,'("Done , total elapsed time", &
            F10.3,"s")') real(current-st)/real(rate)
!         
      stop
      end
