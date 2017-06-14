      Module BEM_medium      
      use readio_medium
      use cav_types
      use pedra_friends
!      use, intrinsic :: iso_c_binding

      implicit none
      real(dbl) :: f_f,f_w    
      real(dbl), allocatable :: cals(:,:),cald(:,:) !Calderon D and S matrices 
      real(dbl), allocatable :: sm1(:,:)
      real(dbl), allocatable :: eigv(:),eigt(:,:)
      real(dbl), allocatable :: sm12(:,:),sp12(:,:)
      real(dbl), allocatable :: matq0(:,:),matqd(:,:)
      real(dbl), allocatable :: matqv(:,:),matqq(:,:)
!      complex(16), allocatable :: matq_omega(:,:)
!SP 22/02/16:  Onsager's factors for reaction and local(x) fields
      real(8) :: f_0,f_d,tau_ons,fx_0,fx_d,taux_ons
      complex(cmp) :: eps,eps_f
      save
      private
      public init_BEM,deallocate_BEM,do_eps,eps,eps_f,eigv,eigt, &
             sm1,f_f,f_w,sm12,f_0,f_d,tau_ons,fx_0,fx_d,taux_ons, &
             matqv,matqq,matqd,matq0,do_TS_matrix, &
             do_charge_freq,allocate_TS_matrix,deallocate_TS_matrix


      contains
!     
      subroutine init_BEM     
      integer(i4b) :: its                    
       if (Fprop.eq.'dip') then
         call do_factors 
       else
         if (Fbem.eq.'wri') then
           !Build the cavity/nanoparticle surface
           if(Fcav.eq.'fil') then 
             call read_cavity_full_file
           elseif(Fcav.eq.'gms') then
             call read_gmsh_file
           else
            if(mdm.eq.'sol') then
             call pedra_int('act')
            elseif (mdm.eq.'nan') then
             call pedra_int('met')
            endif
           endif
           ! write out the cavity/nanoparticle surface
           call output_surf
           write(6,*) "Created output file with surface points"
! SC: 30/09/2016: now explicitly set to zero by Fchr.
           !call write_charges0
           !write(6,*) "Created charges0.inp file with zero charges"
         elseif (Fbem.eq.'rea') then
           ! read in the cavity/nanoparticle surface
           call read_cavity_file
         endif
         if(Fprop.eq."ief") then
           if(Fbem.eq."rea") then
            allocate(matqv(nts_act,nts_act))
            allocate(matqq(nts_act,nts_act))
            allocate(matqd(nts_act,nts_act))
            allocate(matq0(nts_act,nts_act))
           endif
           call do_PCM_propMat
         else
           allocate(sm1(nts_act,nts_act))
           call do_ons_propMat  
         endif
       endif
      return
      end subroutine
!
!
      subroutine output_surf
      integer :: i
! This routine creates the same files also created by Gaussian
      open(unit=7,file="cavity.inp",status="unknown",form="formatted")
       write (7,*) nts_act,nesf_act
       do i=1,nesf_act
         write (7,'(3F22.10)') sfe_act(i)%x,sfe_act(i)%y, &
                               sfe_act(i)%z
       enddo
       do i=1,nts_act
         write (7,'(4F22.10,D14.5)') cts_act(i)%x,cts_act(i)%y,cts_act(i)%z, &
                               cts_act(i)%area,cts_act(i)%rsfe
       enddo
      close(unit=7)
! SC 28/9/2016: write out cavity in GAMESS ineq format
! WHICH UNITS ARE ASSUMED IN GAMESS? CHECK!!
      open(unit=7,file="np_bem.cav",status="unknown",form="formatted")
       write (7,*) nts_act
       do i=1,nts_act
         write (7,'(F22.10)') cts_act(i)%x
       enddo
       do i=1,nts_act
         write (7,'(F22.10)') cts_act(i)%y
       enddo
       do i=1,nts_act
         write (7,'(F22.10)') cts_act(i)%z
       enddo
       do i=1,nts_act
         write (7,'(F22.10)') cts_act(i)%area
       enddo
      close(unit=7)
!      open(unit=7,file="charges0.inp",status="unknown",form="formatted")
!       write (7,*) nts_act
!       do i=1,nts_act
!         write (7,'(2F22.10)') 0.d0, cts_act(i)%area
!       enddo
!      close(unit=7)
      end subroutine
!
      subroutine output_charge_pqr(st)
      integer :: its,i
      character(18) :: fname
      real(8) :: st(:,:)
      real(8) :: area
      area=sum(cts_act(:)%area)
      do i=1,10
       write (fname,'("charge_freq_",I0,".pqr")') i
       open(unit=7,file=fname,status="unknown", &
          form="formatted")
       write (7,*) nts_act
       do its=1,nts_act
         write (7,'("ATOM ",I6," H    H  ",I6,3F11.3,F15.5,"  1.5")') &
              its,its,cts_act(its)%x,cts_act(its)%y,cts_act(its)%z, &
              st(its,i+1)/cts_act(its)%area*area
       enddo
       close(unit=7)
      enddo
      end subroutine
!
      subroutine write_charges0
      integer :: i
      open(unit=7,file="charges0.inp",status="unknown",form="formatted")
      write(7,*) nts_act
      do i=1,nts_act
       write(7,*) 0.d0
      enddo
      close(7)
      end subroutine
!
      subroutine deallocate_BEM     
       if (Fprop.eq.'dip') then
       else
         ! deallocate matrices related to nanoparticle
         if (allocated(matqv)) deallocate(matqv)
         if (allocated(matqq)) deallocate(matqq)
         if (allocated(matqd)) deallocate(matqd)
         if (allocated(matq0)) deallocate(matq0)
         if (allocated(sm1)) deallocate(sm1)
       endif
      return
      end subroutine
!
      subroutine do_factors   
       real(8) :: lambda,r_bc
       ! NB: here assumed prolate cavity with a_cav the largest
       !     fix for generic input!
       r_bc=b_cav/c_cav
       if (r_bc.eq.1.d0) then
         lambda=1./3.
       else
         lambda=-4.d0*pi/(r_bc**2-1)*(r_bc/2.d0/sqrt(r_bc**2-1.)* &
           log((r_bc+sqrt(r_bc**2-1))/(r_bc-sqrt(r_bc**2-1.)))-1.d0)
       endif
       f_0=3.d0*lambda/b_cav/c_cav**2*(eps_0-1.d0)/(eps_0+  &
           lambda/(1.d0-lambda))
       f_d=3.d0*lambda/b_cav/c_cav**2*(eps_d-1.d0)/(eps_d+  &
           lambda/(1.d0-lambda))
       tau_ons=(eps_d+lambda/(1.d0-lambda))/ &
               (eps_0+lambda/(1.d0-lambda))*tau_deb 
       write (6,*) "Onsager"
       write (6,*) "eps_0,eps_d",eps_0,eps_d
       write (6,*) "r_bc",r_bc
       write (6,*) "lambda",lambda
       write (6,*) "f_0",f_0
       write (6,*) "f_d",f_d
       write (6,*) "tau_ons",tau_ons
       if(Floc.eq."loc") then
         fx_0=3.d0/b_cav/c_cav**2*eps_0/(two*eps_0+one)
         fx_d=3.d0/b_cav/c_cav**2*eps_d/(two*eps_d+one)
         taux_ons=(two*eps_d+one)/(two*eps_0+one)
       endif
      return
      end subroutine
!
      subroutine do_PCM_propMat
      integer(i4b) :: i,j
      real(8), allocatable :: scr4(:,:),scr1(:,:)
      real(8), allocatable :: scr2(:,:),scr3(:,:)
      real(dbl), dimension(nts_act) :: fact1,fact2
      real(dbl), dimension(nts_act) :: Kdiag0,Kdiagd
      complex(16), allocatable :: Kdiag_omega(:)
      real(dbl) :: sgn,fac_eps0,fac_epsd
      real(dbl):: temp
! SC 11/08/2016 if this is a writing run, create the matrixes and stop.
! SC TO DO: all the automatic arrays are always created, but some of them
!           are needed for fbem=rea and others for fbem=wri, move to allocatable
      ! Solvent or nanoparticle
      sgn=one                 
      if(mdm.eq."nan") sgn=-one  
!
      If(Fbem.eq.'wri') then
      ! create calderon D, D* and S matrices
       allocate(cals(nts_act,nts_act))
       allocate(cald(nts_act,nts_act))
         do i=1,nts_act
          do j=1,nts_act
            call green_s(i,j,temp)
            cals(i,j)=temp
            call green_d(i,j,temp)
            cald(i,j)=temp
          enddo
         enddo
         open(7,file="mat_SD.inp",status="unknown")
         write(7,*) nts_act
         do j=1,nts_act
          do i=j,nts_act
           write(7,'(2E26.16)')cals(i,j),cald(i,j)
          enddo
         enddo
         close(7)
         deallocate(cals,cald)
         write(6,*) "Matrixes S D have been written out, I'll stop now"
        stop
      elseif(Fbem.eq.'rea') then
!SC 05/11/2016: use the same routine as for frequency calculation
!   to build S^-1/2 (in sm12), S^1/2 (in sp12), T (in eigt) and Lambda (in eigv) 
        call allocate_TS_matrix
        call do_TS_matrix
!
        allocate(scr4(nts_act,nts_act),scr1(nts_act,nts_act))
        allocate(scr2(nts_act,nts_act),scr3(nts_act,nts_act))
!        Test on S^-1/2 D A S^1/2 + S^1/2 A D* S^-1/2 Diagonal passed 
!        
!        Form the Q_w and Q_f and K_d and K_0 for debye propagation
         if (Feps.eq."deb") then
!          debye dielectric function  
           fac_eps0=(eps_0+1)/(eps_0-1)
           Kdiag0(:)=(twp-sgn*eigv(:))/(twp*fac_eps0-sgn*eigv(:))
           fac_epsd=(eps_d+1)/(eps_d-1)
           Kdiagd(:)=(twp-sgn*eigv(:))/(twp*fac_epsd-sgn*eigv(:))
           ! SP: Need to check the signs of the second part for a debye medium localized in space  
           fact1(:)=((twp-sgn*eigv(:))*eps_0+twp+eigv(:))/ &
                    ((twp-sgn*eigv(:))*eps_d+twp+eigv(:))/tau_deb
           fact2(:)=Kdiag0(:)*fact1(:)
         elseif (Feps.eq."drl") then       
!          Drude-Lorentz dielectric function
           Kdiagd(:)=0.d0
           fact2(:)=(twp-sgn*eigv(:))*eps_A/(two*twp)  
! SC: the first eigenvector should be 0 for the NP
           if (mdm.eq.'nan') fact2(1)=0.d0
! SC: no spurious negative square frequencies
           do i=1,nts_act
            if(fact2(i).lt.0.d0) fact2(i)=0d0
           enddo
           if (eps_w0.eq.0.d0) eps_w0=1.d-8
           fact1(:)=fact2(:)+eps_w0*eps_w0  
           Kdiag0(:)=fact2(:)/fact1(:)
! already printed in do_TS_matrix
!           write (6,*) "Squares of resonance frequencies, in a.u."
!           do i=1,nts_act
!            write(6,*) i,fact1(i)
!           enddo
         else
           write(6,*) "Wrong epsilon choice"
           stop
         endif
!
!        Buil S-1/2T in scr3
         scr3=matmul(sm12,eigt)
! SC 31/10/2016: print out total charges associated with eigenvectors
         write (6,*) "Total charge associated to each eigenvector"
         do i=1,nts_act
           write(6,*) i,sum(scr3(:,i))
         enddo
! SC 05/11/2016 write out the transition charges in pqr format
         call output_charge_pqr(scr3)
!
!        Buil T*S1/2 in scr2
         scr2=matmul(transpose(eigt),sp12)  
!        Buil T*S-1/2 in scr4
         scr4=transpose(scr3)
!        Build matqq (R for debye)
         do i=1,nts_act
           scr1(:,i)=scr3(:,i)*fact1(i) 
         enddo
         matqq=matmul(scr1,scr2)
!        Build matqv (\tilde{Q} for debye)
         do i=1,nts_act
           scr1(:,i)=scr3(:,i)*fact2(i) 
         enddo
         matqv=-matmul(scr1,scr4) 
!        Build Q_0 
         do i=1,nts_act
           scr1(:,i)=scr3(:,i)*Kdiag0(i) 
         enddo
         matq0=-matmul(scr1,scr4) 
!        Build Q_d only for debye
         do i=1,nts_act
           scr1(:,i)=scr3(:,i)*Kdiagd(i) 
         enddo
         matqd=-matmul(scr1,scr4) 
!   Print matrices in output
         open(7,file="BEM_matrices.mat",status="unknown")
         write(7,*) nts_act
         do j=1,nts_act
          do i=j,nts_act
           write(7,'(4E26.16)')matqv(i,j),matqd(i,j),matqq(i,j),matq0(i,j)
          enddo
         enddo
         close(7)
         write(6,*) "Written out the propagation BEM matrixes"
         open(7,file="np_bem.mat",status="unknown")
         do j=1,nts_act
          do i=1,nts_act
           write(7,'(D20.12)') matq0(i,j)/cts_act(i)%area
          enddo
         enddo
         close(7)
         write(6,*) "Written out the static matrix for gamess"
         open(7,file="np_bem.mdy",status="unknown")
         do j=1,nts_act
          do i=1,nts_act
           write(7,'(D20.12)') matqd(i,j)/cts_act(i)%area
          enddo
         enddo
         close(7)
         write(6,*) "Written out the static matrix for gamess"
      endif
      deallocate(scr4,scr1)
      deallocate(scr2,scr3)
      call deallocate_TS_matrix
      return
      end subroutine
!
      subroutine green_d (i,j,value)
      integer(i4b), intent(in):: i,j
      real(dbl), intent(out) :: value
      real(dbl):: dist,diff(3)
      if (i.ne.j) then
         diff(1)=(cts_act(i)%x-cts_act(j)%x)
         diff(2)=(cts_act(i)%y-cts_act(j)%y)
         diff(3)=(cts_act(i)%z-cts_act(j)%z)
         dist=sqrt(dot_product(diff,diff))
         value=dot_product(cts_act(j)%n,diff)/dist**3 
      else
         value=-1.0694*sqrt(4.d0*pi*cts_act(i)%area)/(2.d0* &
                cts_act(i)%rsfe)/cts_act(i)%area
      endif
      return
      end subroutine
!
      subroutine green_s (i,j,value)
      integer(i4b), intent(in):: i,j
      real(dbl), intent(out) :: value
      real(dbl):: dist,diff(3)
      if (i.ne.j) then
        diff(1)=(cts_act(i)%x-cts_act(j)%x)
        diff(2)=(cts_act(i)%y-cts_act(j)%y)
        diff(3)=(cts_act(i)%z-cts_act(j)%z)
        dist=sqrt(dot_product(diff,diff))
        value=one/dist 

      else
        value=1.0694*sqrt(4.d0*pi/cts_act(i)%area)
      endif
      return
      end subroutine
!     
      subroutine do_ons_propMat   
      integer(i4b) :: i,j
      real(dbl) :: temp
! Same read/write phylosophy as for ief
       allocate(cals(nts_act,nts_act))
       if(Fbem.eq."wri") then
      ! create calderon  S matrix
        do i=1,nts_act
         do j=1,nts_act
           call green_s(i,j,temp)
           cals(i,j)=temp
         enddo
        enddo
!        Form S^-1
        sm1=inv(cals)
        open(7,file="S_matrix.inp",status="unknown")
        write(7,*) nts_act
        do j=1,nts_act
         do i=j,nts_act
          write(7,'(4E26.16)')sm1(i,j)
         enddo
        enddo
        close(7)
        write(6,*) "Written out the propagation S matrix, now I'll stop"
        stop
       elseif(Fbem.eq."rea") then
        open(7,file="S_matrix.inp",status="old")
        write(7,*) nts_act
        do j=1,nts_act
         do i=j,nts_act
          read(7,*) sm1(i,j)
         enddo
        enddo
        close(7)
        write(6,*) "Read in the propagation S matrix"
        ! Form the f_w=-1/tau+w_0^2 and f_f=A/3                           
        f_w=-(eps_A/3.+eps_w0*eps_w0) 
        f_f=-eps_A/3. 
       endif
       if(allocated(cals)) deallocate(cals)
       return
      end subroutine
!
      subroutine do_eps      
       !eps_gm=eps_gm+f_vel/sfe_act(1)%r
       eps=dcmplx(eps_A,zero)/dcmplx(eps_w0**2-omega**2,-omega*eps_gm)
       eps=eps+onec
       eps_f=(eps-onec)/(eps+twoc)
      return
      end subroutine
!
      function inv(A) result(Ainv)
        real(dbl), dimension(:,:), intent(in) :: A
        real(dbl), dimension(size(A,1),size(A,2)) :: Ainv
      
        real(dbl), dimension(size(A,1)) :: work  ! work array for LAPACK
        integer(i4b), dimension(size(A,1)) :: ipiv   ! pivot indices
        integer(i4b) :: n, info
      
        ! External procedures defined in LAPACK
        external DGETRF
        external DGETRI
      
        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)
      
        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call DGETRF(n, n, Ainv, n, ipiv, info)
      
        if (info /= 0) then
           stop 'Matrix is numerically singular!'
        end if
      
        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call DGETRI(n, Ainv, n, ipiv, work, n, info)
      
        if (info /= 0) then
           stop 'Matrix inversion failed!'
        end if
      end function inv
!
      subroutine allocate_TS_matrix
      allocate(eigv(nts_act))
      allocate(eigt(nts_act,nts_act))
      allocate(sm12(nts_act,nts_act))
      allocate(sp12(nts_act,nts_act))
      return
      end subroutine
!
      subroutine deallocate_TS_matrix
      deallocate(eigv)
      deallocate(eigt)
      deallocate(sm12)
      deallocate(sp12)
      return
      end subroutine
!
      subroutine do_TS_matrix
      integer(i4b) :: i,j,info,lwork,liwork
      real(8), allocatable :: scr1(:,:),scr2(:,:),eigt_t(:,:)
      real(dbl) :: sgn,fac_eps0,fac_epsd
      real(dbl):: temp,fact1,fact2
      character jobz,uplo
      integer(i4b), allocatable :: iwork(:)
      real(8),allocatable :: work(:)
      allocate(scr1(nts_act,nts_act),scr2(nts_act,nts_act))
      allocate(eigt_t(nts_act,nts_act))
      allocate(work(1+6*nts_act+2*nts_act*nts_act))
      allocate(iwork(3+5*nts_act))
      allocate(cals(nts_act,nts_act))
      allocate(cald(nts_act,nts_act))
      sgn=one                 
      if(mdm.eq."nan") sgn=-one  
!
      open(7,file="mat_SD.inp",status="old")
      read(7,*) nts_act
      do j=1,nts_act
       do i=j,nts_act
        read(7,*) cals(i,j), cald(i,j)
        cals(j,i)=cals(i,j)
        cald(j,i)=cald(i,j)
       enddo
      enddo
      close(7)
      write(6,*) "Done reading in S and D"
      ! Form S^1/2 and S^-1/2
      ! Set parameters for diagonalization routine dsyevd
      jobz = 'V'
      uplo = 'U'
      lwork = 1+6*nts_act+2*nts_act*nts_act
      liwork = 3+5*nts_act
      ! Copy the matrix in the eigenvector matrix
      eigt = cals
      call dsyevd (jobz,uplo,nts_act,eigt,nts_act,eigv,work,lwork, &
        iwork,liwork,info)
      do i=1,nts_act
        if(eigv(i).lt.0.d0) then
         write(6,*) "WARNING:",i," eig of S is negative!"
         write(6,*) "   I put it to 1e-8"
         eigv(i)=1.d-8
        endif
        scr1(:,i)=eigt(:,i)*sqrt(eigv(i))
      enddo
      eigt_t=transpose(eigt)
      sp12=matmul(scr1,eigt_t)                   
      do i=1,nts_act
        scr1(:,i)=eigt(:,i)/sqrt(eigv(i))
      enddo
      sm12=matmul(scr1,eigt_t)                   
!     Test on S Diagonal passed 
!     
!     Form the S^-1/2 D A S^1/2 + S^1/2 A D* S^-1/2 , and diagonalize it
      !S^-1/2 D A S^1/2
      do i=1,nts_act
        scr1(:,i)=cald(:,i)*cts_act(i)%area
      enddo
      scr2=matmul(matmul(sm12,scr1),sp12)                   
      !S^-1/2 D A S^1/2+S^1/2 A D* S^-1/2 and diagonalise
      do j=1,nts_act
       do i=1,nts_act
        eigt(i,j)=0.5*(scr2(i,j)+scr2(j,i))
       enddo
      enddo
      call dsyevd (jobz,uplo,nts_act,eigt,nts_act,eigv,work,lwork, &
        iwork,liwork,info)
      open(7,file="TSSK_matrices.mat",status="unknown")
      write(7,*) nts_act
      do j=1,nts_act
       do i=j,nts_act
        write(7,'(4E26.16)')eigt(i,j),sp12(i,j),sm12(i,j)
       enddo
      enddo
      do j=1,nts_act
        write(7,*)eigv(j)
      enddo
      close(7)
      if (Feps.eq.'drl') then
! Do Q0 and Qd for Gamess
       scr1=matmul(sm12,eigt)
       scr2=transpose(scr1)
       write (6,*) "Squares of resonance frequencies, in a.u."
       do i=1,nts_act
        fact2=(twp-sgn*eigv(i))*eps_A/(two*twp)  
        fact1=fact2+eps_w0*eps_w0  
        if (fact1.lt.0.d0) eigv(i)=-twp
        write(6,*) i,fact1
        scr1(:,i)=scr1(:,i)*fact2/fact1
       enddo
       eigt_t=-matmul(scr1,scr2)
         open(7,file="np_bem.mat",status="unknown")
         do j=1,nts_act
          do i=1,nts_act
           write(7,'(D20.12)') eigt_t(i,j)/cts_act(i)%area
          enddo
         enddo
         close(7)
         write(6,*) "Written out the static matrix for gamess"
         open(7,file="np_bem.mdy",status="unknown")
         do j=1,nts_act
          do i=1,nts_act
           write(7,'(D20.12)') 0.d0
          enddo
         enddo
         close(7)
         write(6,*) "Written out the dynamic matrix for gamess"
      endif
      if(allocated(cals).and.allocated(cald)) deallocate(cals,cald)
      deallocate(work)
      deallocate(scr1,scr2,eigt_t)
      write(6,*) "Done setting up T and S" 
      return
      end subroutine
!
      subroutine do_charge_freq(omega_a,pot)
      integer(4) :: its
      complex(16), allocatable :: q_omega(:)
      complex(16) :: mu_omega(3)
      real(8) :: omega_a
      real(8) :: pot(:)
!      real(8), allocatable :: fact1(:),fact2(:)
      complex(16), allocatable :: Kdiag_omega(:)
      integer(4) :: sgn,i
      allocate(Kdiag_omega(nts_act))
!      allocate(fact1(nts_act))
!      allocate(fact2(nts_act))
!      allocate(scr3(nts_act,nts_act))
!      allocate(scr1(nts_act,nts_act))
      sgn=-1
      do i=1,nts_act
       Kdiag_omega(i)=(twp-sgn*eigv(i))/ &
          ((1.d0-omega_a*(omega_a+ui*eps_gm)*2.d0/eps_A)*twp-sgn*eigv(i))
!       if (i.eq.2) write(6,*) omega_a,Kdiag_omega(i)
      enddo
!       write(6,*) 'i,Kdiag_omega(i)'
!      do i=1,nts_act
!       write(6,*) i,Kdiag_omega(i)
!      enddo
      allocate (q_omega(nts_act))
      q_omega=matmul(sm12,pot)
      q_omega=matmul(transpose(eigt),q_omega)
      q_omega=Kdiag_omega*q_omega
      q_omega=matmul(eigt,q_omega)
      q_omega=-matmul(sm12,q_omega)
!      write(6,*) "q_omega"
      mu_omega=0.d0
      do its=1,nts_act
!       write(6,*) its,q_omega(its)
       mu_omega(1)=mu_omega(1)+q_omega(its)*(cts_act(its)%x)
       mu_omega(2)=mu_omega(2)+q_omega(its)*(cts_act(its)%y)
       mu_omega(3)=mu_omega(3)+q_omega(its)*(cts_act(its)%z)
      enddo
      write (6,'(7d15.6)') omega_a,real(mu_omega(:)),aimag(mu_omega(:))
      deallocate (q_omega)
!      scr3=matmul(sm12,eigt)
!      do i=1,nts_act
!        scr1(:,i)=scr3(:,i)*Kdiag_omega(i) 
!      enddo
!      deallocate(Kdiag_omega)
!      matq_omega=-matmul(scr1,transpose(scr3)) 
!      deallocate(scr1,scr3)
      return
      end subroutine
!
      end module
