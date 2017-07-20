      module spectra  
      use cav_types
      use readio
      use readio_medium
      use pedra_friends
      
      use, intrinsic :: iso_c_binding

      implicit none
      real(dbl), allocatable:: Sdip(:,:,:), Sfld(:,:)
      double complex, parameter :: zeroc=(zero,zero)                
      save
      private
      public Sdip, Sfld, do_spectra, init_spectra, read_arrays
!
      contains
!
      subroutine init_spectra
       allocate (Sdip(n_step,3,3),Sfld(n_step,3))
       Sdip(:,:,:)=0.0d0
       Sfld(:,:)=0.0d0

      return
      end subroutine
!
      subroutine read_arrays  
      integer(4) :: file_mol=10,file_fld=8,file_med=9,i,x
      real(8) :: t
       character(20) :: name_f
       write(name_f,'(a5,i0,a4)') "mu_t_",n_f,".dat"
       open (file_mol,file=name_f,status="unknown")
       if (mdm.ne.'vac') then 
         write(name_f,'(a9,i0,a4)') "medium_t_",n_f,".dat"
         open (file_med,file=name_f,status="unknown")
       endif
       write(name_f,'(a5,i0,a4)') "field",n_f,".dat"
       open (file_fld,file=name_f,status="unknown")
       read(file_mol,*)
       !read(file_fld,*)
       if (mdm.ne.'vac') read(file_med,*)
       do i=1,n_step
         read (file_mol,'(i8,f14.4,3e22.10)') x,t,Sdip(i,:,1)       
         if (mdm.ne.'vac') then
           read (file_med,'(i8,f12.2,4e22.10)') x,t,Sdip(i,:,2)       
         endif
         read (file_fld,'(f12.2,3e22.10e3)') t,Sfld(i,:) 
       enddo
       close(file_mol)
       close(file_fld)
       if (mdm.ne.'vac') close(file_med)
      return
      end subroutine
!
      subroutine do_spectra      
      implicit none
      real(dbl), allocatable :: Dinp(:),Finp(:)
      real(dbl) :: dw,fac,absD,refD,phiF,phiD,modF,modD,mdl    
      real(dbl) :: Deq(3),Deq_np(3)    
      integer(i4b) :: i,isp,vdim,istart  
      integer*8 plan
      double complex, allocatable :: Doutp(:),Foutp(:),src       
      character(len=15):: fname
! SC 15/01/2016: changed Makefile from Silvio's version
!      find a better way to include this file
!      than changing source back and forth from f to f03
!      include '/usr/include/fftw3.f03'
!      include '/usr/local/include/fftw3.f03'
      include 'fftw3.f'
    
      ! This is needed to improve the delta_omega and also the quality 
      ! of the DFT with respect to the FT
      istart=int(start/dt)
      vdim=n_step-istart
      allocate (Dinp(vdim),Finp(vdim)) 
      allocate (Doutp(int(dble(vdim)/2)+1))
      allocate (Foutp(int(dble(vdim)/2)+1)) 
      dw=2*pi/dble(vdim)/dt
      ! The minus sign is for the electronic negative charge
      Deq(:)=Sdip(1+istart,:,1)
      Deq_np(:)=Sdip(1+istart,:,2)
      do i=1,vdim
        Sdip(i+istart,:,1)=(Sdip(i+istart,:,1)-Deq(:))*    & 
                                     exp(-i*dt/tau(1))
        Sdip(i+istart,:,2)=(Sdip(i+istart,:,2)-Deq_np(:))* & 
                                     exp(-i*dt/tau(2))
        Sdip(i+istart,:,3)=Sdip(i+istart,:,1)+Sdip(i+istart,:,2)
      enddo
      ! SP 28/10/16: normalize the dir_ft
      mdl=dir_ft(1)*dir_ft(1)+dir_ft(2)*dir_ft(2)+dir_ft(3)*dir_ft(3) 
      dir_ft(:)=dir_ft(:)/sqrt(mdl)
      do isp=1,3
        Doutp=zeroc
        Foutp=zeroc
        do i=1,vdim
          ! SP 28/10/16: FT in the dir_ft direction
          Dinp(i)=dot_product(Sdip(i+istart,:,isp),dir_ft(:)) 
          Finp(i)=dot_product(Sfld(i+istart,:),dir_ft(:))
        enddo
        call dfftw_plan_dft_r2c_1d(plan,vdim,Dinp,Doutp,FFTW_ESTIMATE)
        call dfftw_execute_dft_r2c(plan, Dinp, Doutp)
        call dfftw_destroy_plan(plan)
        call dfftw_plan_dft_r2c_1d(plan,vdim,Finp,Foutp,FFTW_ESTIMATE)
        call dfftw_execute_dft_r2c(plan, Finp, Foutp)
        call dfftw_destroy_plan(plan)
        if (isp.eq.1) &
            write(fname,'(a7,i0,a4)') "sp_mol_",n_f,".dat"
        if (isp.eq.2) &
            write(fname,'(a6,i0,a4)') "sp_np_",n_f,".dat"
        if (isp.eq.3) &
            write(fname,'(a9,i0,a4)') "sp_molnp_",n_f,".dat"
        open(unit=15,file=fname,status="unknown",form="formatted")
        do i=1,int(vdim/2)+1  
          modD=sqrt(real(Doutp(i))**2+aimag(Doutp(i))**2)
          modF=sqrt(real(Foutp(i))**2+aimag(Foutp(i))**2)
          phiD=atan2(aimag(Doutp(i)),real(Doutp(i)))
          phiF=atan2(aimag(Foutp(i)),real(Foutp(i)))
          absD=-(modD/modF)*sin(phiD-phiF)
          refD=(modD/modF)*cos(phiD-phiF)
!          src=1./Foutp(i)
!          absD=aimag(Doutp(i)*src)
!          refD=real(Doutp(i)*src)
          write(15,'(3e20.10)') (i-1)*dw, absD, refD
        enddo 
        close(unit=15)
      enddo
      deallocate (Dinp,Finp) 
      deallocate (Doutp,Foutp) 
      deallocate (Sdip,Sfld)
      return
      end subroutine
!
      end module
