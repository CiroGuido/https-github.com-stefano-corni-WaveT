program twodfft

 use constants 
 use, intrinsic :: iso_c_binding

 implicit none

 integer(i4b)               :: dim1,dim2,dim3,i,j,plan,ntot
 integer(i4b)               :: jj,n3,idum,dim3e,m,dscan1
 real(dbl)                  :: rdum,dw1,dw3,dt1,dt3,modd
 real(dbl)                  :: dir(3),mu(3),sigma3
 real(dbl), allocatable     :: dinp(:,:,:),inp(:,:) 
 complex(cmp), allocatable  :: doutp(:,:) 
 character*20               :: filename
 character*100              :: cdum

 include 'fftw3.f03'

 namelist /fft/ dscan1,dim1,dim2,dim3,dt1,dt3,dir,sigma3

 read(*,nml=fft)

 write(*,*) '**********************************************'
 write(*,*) '*                                            *'
 write(*,*) '*               Tool for a                   *'
 write(*,*) '*                 2D FFT                     *'
 write(*,*) '*                                            *'
 write(*,*) '*                   by                       *'
 write(*,*) '*              Stefano Corni                 *'
 write(*,*) '*                  and                       *'
 write(*,*) '*             Emanuele Coccia                *'
 write(*,*) '*                                            *'
 write(*,*) '**********************************************'


 if (dscan1.lt.0) then
    write(*,*) 'ERROR: dscan1 must be positive'
    stop
 endif 

 if (dim1.lt.0) then
    write(*,*) 'ERROR: dim1 must be positive'
    stop
 endif

 if (dim2.lt.0) then
    write(*,*) 'ERROR: dim2 must be positive'
    stop
 endif

 if (dim3.lt.0) then
    write(*,*) 'ERROR: dim3 must be positive'
    stop
 endif

 if (dt1.le.0) then
    write(*,*) 'ERROR: dt1 must be larger than zero'
    stop
 endif

 if (dt3.le.0) then
    write(*,*) 'ERROR: dt3 must be larger than zero'
    stop
 endif

 if (sigma3.lt.0) then
    write(*,*) 'ERROR: sigma3 must be positive'
    stop
 endif 

 if (dir(1).le.0.and.dir(2).le.0.and.dir(3).le.0) then
    write(*,*) 'ERROR: One component of dir_ft must be larger than zero'
    stop
 endif

 write(*,*) ''
 write(*,*) 'Number of points for t1', dscan1
 write(*,*) 'Number of points for t3', dim3
 write(*,*) 'Time step along the first dimension', dt1
 write(*,*) 'Time step along the second dimension', dt3
 write(*,*) 'Direction for FFT', dir(1),dir(2),dir(3)
 write(*,*) 'Sigma of the enevelope of the third pulse', sigma3
 write(*,*) ''
 
 n3=int((2.d0*sigma3)/dt3)
 dim3e=dim3+n3

 if(mod(dim3e,2).gt.0) dim3e=dim3e-1
 if(mod(dscan1,2).gt.0) dscan1=dscan1-1

 allocate (dinp(dscan1,dim3e,3))
 allocate (inp(dscan1,dim3e))
 !allocate (doutp(int(dble(dscan1)/two)+1,int(dble(dim3e)/two)+1))
 allocate (doutp(int(dble(dscan1)/two)+1,dim3e))
 dw1=2*pi/dble(dscan1)/dt1
 dw3=2*pi/dble(dim3e)/dt3

 ntot=dim1+dim2+dim3

 open(8,file='2d_spectrum.dat',status="unknown")

 do m=1,dscan1

    if (m.lt.10) then
        WRITE(filename,'(a,i1.1,a)') "mu_t_",m,".dat"
    elseif (m.ge.10.and.m.lt.100) then
        WRITE(filename,'(a,i2.2,a)') "mu_t_",m,".dat"
    elseif (m.ge.100.and.m.lt.1000) then
        WRITE(filename,'(a,i3.3,a)') "mu_t_",m,".dat"
    elseif (m.ge.1000.and.m.lt.10000) then
        WRITE(filename,'(a,i4.4,a)') "mu_t_",m,".dat"
    elseif (m.ge.10000.and.m.lt.100000) then
        WRITE(filename,'(a,i5.5,a)') "mu_t_",m,".dat"
    endif

    open(7,file=filename)

    dinp(:,:,:)=0.d0
    jj=0
    read(7,*) cdum
    do i=1,ntot
       read(7,*) idum, rdum, mu(:) 
       if (i.ge.(dim1+dim2-n3)) then
          jj=jj+1 
          dinp(m,jj,1) = mu(1)
          dinp(m,jj,2) = mu(2)
          dinp(m,jj,3) = mu(3) 
       endif
    enddo

    close(7)

 enddo

 doutp=cmplx(0.d0,0.d0)

 dir(:)=dir(:)/sqrt(dot_product(dir,dir))

 do j=1,dim3e
    do i=1,dscan1
       inp(i,j)=dot_product(dinp(i,j,:),dir(:))
    enddo
 enddo

 !Homodyne detection
 inp(:,:) = inp(:,:)**2

 call dfftw_plan_dft_r2c_2d(plan,dscan1,dim3e,inp,doutp,FFTW_ESTIMATE)
 call dfftw_execute_dft_r2c(plan,inp,doutp)
 call dfftw_destroy_plan(plan)

 do i=1,int(dscan1/two)
    do j=1,int(dim3e/two)
       !modd=sqrt(real(doutp(i,j))**2+aimag(doutp(i,j))**2)
       modd=aimag(doutp(i,j)) 
       write(8,'(3e20.10)') (i-1)*dw1/ev_to_au, (j-1)*dw3/ev_to_au, modd
    enddo 
 enddo

 close(8)

 deallocate(dinp)
 deallocate(doutp)
 deallocate(inp)

 write(*,*) 'End of the calculation'
 
 stop


end program twodfft
