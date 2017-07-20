module parallel 

  use readio
  use readio_medium
  use, intrinsic :: iso_c_binding

  integer(4)  :: myrank,nproc,ierr_mpi,init_addr,final_addr,idest,isource


  contains

  subroutine init_MPI()
!------------------------------------------------------------------------
! @brief Initialize the MPI environment 
!
! @date Created   : E. Coccia 18 May 2017
! Modified  :
! @param 
!------------------------------------------------------------------------

    call mpi_init(ierr_mpi)
    call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr_mpi)
    call mpi_comm_rank(MPI_COMM_WORLD,nproc,ierr_mpi)
    call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)

    return 

  end subroutine init_MPI


  subroutine bcast_MPI()
!------------------------------------------------------------------------
! @brief Broadcast variable values 
!
! @date Created   : E. Coccia 18 May 2017
! Modified  :
! @param 
!------------------------------------------------------------------------

    call mpi_bcast(,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
    call mpi_bcast(,3,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr_mpi)

    return

  end subroutine bcast_MPI


  subroutine do_avg_MPI(c)
!------------------------------------------------------------------------
! @brief Master does the sum for population and coherence
!        at each selected step 
!        From slaves to master
!
! @date Created   : E. Coccia 18 May 2017
! Modified  :
! @param 
!------------------------------------------------------------------------

     if (myrank.ne.0) then
        call mpi_send(c,n_ci,mpi_complex,0,myrank,mpi_comm_world,ierr_mpi) 
     else
        call mpi_recv(c,n_ci,mpi_complex)  
     endif

     return

  end subroutine do_avg_MPI


  subroutine final_MPI()
!------------------------------------------------------------------------
! @brief FInalize the MPI environment 
!
! @date Created   : E. Coccia 18 May 2017
! Modified  :
! @param 
!------------------------------------------------------------------------

   call mpi_finalize(ierr_mpi)

   return

  end subroutine final_MPI

end module parallel
