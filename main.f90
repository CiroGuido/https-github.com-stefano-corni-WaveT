      program tdcis
      use readio  
      use readio_medium
      use propagate    
      use QM_coupling  
      use spectra
      use dissipation       
      !use parallel 

      implicit none
!
!#ifdef PARALLEL
!      include "mpif.h"
!#endif

      integer :: st,current,rate

!#ifdef PARALLEL
!     !initialize the MPI environment
!     call init_MPI()
!     call mpi_init(ierr_mpi)
!     call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr_mpi)
!     call mpi_comm_rank(MPI_COMM_WORLD,nproc,ierr_mpi) 
!     call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
!#endif

!      if (myrank.eq.0) then
!
!     read in the input parameter for the present evolution
         call system_clock(st,rate)
         call read_input
         if (mdm.ne."vac") call read_medium
!      endif
      call init_spectra
!
!     create the field 
      call create_field
!      if (myrank.eq.0) then
         call system_clock(current)
         write(6,'("Done reading input & setting up the field, took", &
            F10.3,"s")') real(current-st)/real(rate)
!      endif
!
!     propagate or diagonalise matrix
      if(nmodes.eq.0) then
        call prop
        call do_spectra
      else 
        call do_diag    
      endif

!      if (myrank.eq.0) then
         call system_clock(current)
         write(6,'("Done , total elapsed time", &
            F10.3,"s")') real(current-st)/real(rate)
!      endif
!         
!#ifdef PARALLEL
    !finalize the MPI environment
!    call final_MPI(ierr_mpi)
!#endif
      stop
      end
