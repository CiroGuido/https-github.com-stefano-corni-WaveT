      program tdplas
      use readio  
      use readio_medium
      use BEM_medium
      implicit none
      integer :: st,current,current1,rate
!
!     read in the input parameter for the present evolution
      call system_clock(st,rate)
      call read_medium_tdplas
      call system_clock(current)
      write(6,'("Done reading input, took", &
            F10.3,"s")') real(current-st)/real(rate)
!
!     create cavity and S and D matrix
      call init_BEM
      call system_clock(current1)
      write(6,'("Done setting S and D, took", &
            F10.3,"s")') real(current1-current)/real(rate)
      if (Fgam.eq.'yes') then
       call allocate_TS_matrix
       call do_TS_matrix
       call deallocate_TS_matrix
      endif
      call system_clock(current)
      write(6,'("Done, total elapsed time", &
            F10.3,"s")') real(current-st)/real(rate)
!         
      stop
      end
