      program tdcis
      use readio  
      use spectra
      use dissipation
      use propagate    
      use global_tdplas, only: set_global_tdplas
      use readio_medium
      use QM_coupling
  
      implicit none
      integer :: st,current,rate
!
!     read in the input parameter for the present evolution
      call system_clock(st,rate)
      call read_input
      call set_global_tdplas(dt,mdm,mol_cc,n_ci,n_ci_read,c_i,e_ci,mut,fmax,omega,Ffld,n_out,n_f)
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
      end program tdcis
