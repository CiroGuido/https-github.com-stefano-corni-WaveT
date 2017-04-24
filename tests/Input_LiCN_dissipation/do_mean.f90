program mean

 implicit none
 integer                 :: nrep,nstates,istate,nsteps,itype
 real(8), allocatable    :: rdum(:,:,:), sum(:), err(:)
 integer, allocatable    :: i(:,:)
 real(8), allocatable    :: t(:,:)
 real(8), allocatable    :: idum(:,:,:), isum(:), ierr(:) 
 
 integer :: j,k,m 
 real(8) :: modavg, phavg, ferr

 read(*,*) nsteps
 read(*,*) nrep
 read(*,*) nstates
 read(*,*) istate
 read(*,*) itype

 open(10,file='tmp.dat')

 allocate(i(nsteps,nrep))
 allocate(t(nsteps,nrep))
 allocate(rdum(nsteps,nstates,nrep))
 allocate(sum(nsteps), err(nsteps))

 if (itype.eq.0) then

   open(11,file='pop.dat')

   do j=1,nsteps
      read(10,*) (i(j,m),t(j,m), (rdum(j,k,m), k=1,nstates),m=1,nrep)
   enddo

 ! Mean
   do j=1,nsteps
      sum(j)=0.d0
      do m=1,nrep
         sum(j) = sum(j) + rdum(j,istate,m) 
      enddo
   enddo
   sum=sum/real(nrep)

 ! Error 
   do j=1,nsteps
      err(j) = 0.d0
      do m=1,nrep
         err(j) = err(j) + (rdum(j,istate,m) - sum(j))**2
      enddo
   enddo
   err=err/real(nrep-1.d0)
   err=sqrt(err)

   do j=1,nsteps
      write(11,*) i(j,1), t(j,1), sum(j), err(j)/sqrt(real(nrep))
   enddo

 elseif (itype.eq.1) then

   open(11,file='dec_m.dat')
   open(12,file='dec_p.dat')

   allocate(idum(nsteps,nstates,nrep))
   allocate(isum(nsteps),ierr(nsteps))


   do j=1,nsteps
      read(10,*) (i(j,m),t(j,m), (rdum(j,k,m), idum(j,k,m), k=1,nstates),m=1,nrep)
   enddo

   ! Mean
   do j=1,nsteps
      sum(j)=0.d0
      isum(j)=0.d0
      do m=1,nrep
         sum(j) = sum(j) + rdum(j,istate,m)
         isum(j) = isum(j) + idum(j,istate,m)
      enddo
   enddo
   sum=sum/real(nrep)
   isum=isum/real(nrep)

 ! Error 
   do j=1,nsteps
      err(j) = 0.d0
      ierr(j) = 0.d0
      do m=1,nrep
         err(j) = err(j) + (rdum(j,istate,m) - sum(j))**2
         ierr(j) = ierr(j) + (idum(j,istate,m) - isum(j))**2 
      enddo
   enddo
   err=err/real(nrep-1.d0)
   ierr=ierr/real(nrep-1.d0)
   err=sqrt(err)
   ierr=sqrt(ierr)

   do j=1,nsteps
      modavg = sqrt(sum(j)**2 + isum(j)**2)
      phavg = atan2(isum(j),sum(j)) 
      ferr = sqrt( (sum(j)/modavg*err(j))**2 + (isum(j)/modavg*ierr(j))**2 )
      ferr = ferr/sqrt(real(nrep)) 
      write(11,*) i(j,1), t(j,1), modavg, ferr !, err(j)/sqrt(real(nrep))
      write(12,*) i(j,1), t(j,1), phavg  !, err(j)/sqrt(real(nrep))
   enddo

   deallocate(idum)
   deallocate(isum)
   deallocate(ierr)

   close(12)

 endif

 deallocate(i) 
 deallocate(t)
 deallocate(rdum)
 deallocate(err)
 deallocate(sum)


 close(10)
 close(11)

end program mean
