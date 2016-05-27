      Module BEM_medium      
      use readio
      use readio_medium
      use cav_types
      use pedra_friends
      use spectra
!      use, intrinsic :: iso_c_binding

      implicit none
      real(dbl), parameter :: TOANGS=0.52917724924D+00
      real(dbl), parameter :: ANTOAU=1.0D+00/TOANGS
!LIGHT SPEED IN AU TO BE REVISED
      real(dbl), parameter :: c=1.37036d2                
      real(dbl) :: f_f,f_w    
      real(dbl), allocatable :: cals(:,:),cald(:,:),calds(:,:)
      real(dbl), allocatable :: sm1(:,:)
      real(dbl), allocatable :: eigv(:),eigt(:,:)
      real(dbl), allocatable :: sm12(:,:)
!SP 22/02/16:  Onsager's factors for reaction and local(x) fields
      real(8) :: f_0,f_d,tau_ons,fx_0,fx_d,taux_ons
      complex(cmp) :: eps,eps_f
      save
      private
      public init_BEM,deallocate_BEM,do_eps,eps,eps_f,eigv,eigt, &
             sm1,f_f,f_w,sm12,f_0,f_d,tau_ons,fx_0,fx_d,taux_ons

      contains
!     
      subroutine init_BEM     
      integer(i4b) :: its                    
       if (Fprop.eq.'dip') then
         call do_factors 
       else
         if (.not.readmat) then
           ! Build Matrices for propagation
           allocate(cals(nts_act,nts_act))
           allocate(cald(nts_act,nts_act))
           allocate(calds(nts_act,nts_act))
           allocate(eigv(nts_act))
           allocate(eigt(nts_act,nts_act))
           call do_calderon
           if(Fprop.eq."ief") then
             if (.not.allocated(matqf)) allocate(matqf(nts_act,nts_act))
             if (.not.allocated(matqw)) allocate(matqw(nts_act,nts_act))
             if (.not.allocated(matqd)) allocate(matqd(nts_act,nts_act))
             if (.not.allocated(matq0)) allocate(matq0(nts_act,nts_act))
             if (.not.allocated(sm12)) allocate(sm12(nts_act,nts_act))
             call do_PCMdiag  
           else
             if (.not.allocated(sm1)) allocate(sm1(nts_act,nts_act))
             call do_ons  
           endif
           deallocate(cals,cald,calds)
           deallocate(eigt,eigv)
         endif
       endif
      return
      end subroutine
!
      subroutine deallocate_BEM     
       if (Fprop.eq.'dip') then
       else
         ! deallocate matrices related to nanoparticle
         if (allocated(matqf)) deallocate(matqf)
         if (allocated(matqw)) deallocate(matqw)
         if (allocated(matqd)) deallocate(matqd)
         if (allocated(matq0)) deallocate(matq0)
         if (allocated(sm12)) deallocate(sm12)
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
       if(localf) then
         fx_0=3.d0/b_cav/c_cav**2*eps_0/(two*eps_0+one)
         fx_d=3.d0/b_cav/c_cav**2*eps_d/(two*eps_d+one)
         taux_ons=(two*eps_d+one)/(two*eps_0+one)
       endif
      return
      end subroutine
!
      subroutine do_PCMdiag   
       integer(i4b) :: i,j,info,lwork,liwork
       real(dbl), dimension(nts_act,nts_act) :: scr4,scr1,sp12
       real(dbl), dimension(nts_act,nts_act) :: scr2,scr3
       real(dbl), dimension(nts_act) :: fact1,fact2
       real(dbl), allocatable :: Kdiag0(:),Kdiagd(:)
       real(dbl) :: sgn,fac_eps0,fac_epsd
       character jobz,uplo
       integer(i4b) :: iwork(3+5*nts_act)
       real(dbl) :: work(1+6*nts_act+2*nts_act*nts_act)
!      
!      solvent or nanoparticle
       sgn=one                 
       if(mdm.eq."nan") sgn=-one  
!      Set parameters for diagonalization routine dsyevd
       jobz = 'V'
       uplo = 'U'
       lwork = 1+6*nts_act+2*nts_act*nts_act
       liwork = 3+5*nts_act
       open(7,file="DSD_matrices.inp",status="unknown")
       write(7,*) nts_act
       do j=1,nts_act
        do i=j,nts_act
         write(7,'(3E26.16)')cals(i,j),cald(i,j),calds(i,j)
        enddo
       enddo
       close(7)
!      
!      Form S^1/2 and S^-1/2
!      Copy the matrix in the eigenvector matrix
       eigt = cals
       call dsyevd (jobz,uplo,nts_act,eigt,nts_act,eigv,work,lwork, &
         iwork,liwork,info)
       do i=1,nts_act
         scr1(:,i)=eigt(:,i)*sqrt(eigv(i))
       enddo
       sp12=matmul(scr1,transpose(eigt))                   
       do i=1,nts_act
         scr1(:,i)=eigt(:,i)/sqrt(eigv(i))
       enddo
       sm12=matmul(scr1,transpose(eigt))                   
!      Test on S Diagonal passed 
!      
!      Form the S^-1/2 D A S^1/2 + S^1/2 A D* S^-1/2 , and diagonalize it
       !S^-1/2 D A S^1/2
       do i=1,nts_act
         scr1(:,i)=cald(:,i)*cts_act(i)%area
       enddo
       scr1=matmul(matmul(sm12,scr1),scr4)                   
       !S^1/2 A D* S^-1/2
       do i=1,nts_act
         scr2(i,:)=calds(i,:)*cts_act(i)%area
       enddo
       scr2=matmul(matmul(sp12,scr2),sm12)                   
       !S^-1/2 D A S^1/2+S^1/2 A D* S^-1/2 and diagonalise
       !eigt(:,:)=0.5*(scr1(:,:)+scr2(:,:))
       eigt(:,:)=scr1(:,:)
       call dsyevd (jobz,uplo,nts_act,eigt,nts_act,eigv,work,lwork, &
         iwork,liwork,info)
!      Test on S^-1/2 D A S^1/2 + S^1/2 A D* S^-1/2 Diagonal passed 
!      
!      Form the Q_w and Q_f and K_d and K_0 for debye propagation
       if (Feps.eq."deb") then
!        debye dielectric function  
         ! SP: Need to check the signs of the second part for a debye medium localized in space  
         allocate(Kdiag0(nts_act),Kdiagd(nts_act))
         fact1(:)=((twp-sgn*eigv(:))*eps_0+twp+eigv(:))/ &
                  ((twp-sgn*eigv(:))*eps_d+twp+eigv(:))/tau_deb
         fac_eps0=(eps_0+1)/(eps_0-1)
         Kdiag0(:)=(twp-sgn*eigv(:))/(twp*fac_eps0-sgn*eigv(:))
         fact2(:)=Kdiag0(:)*fact1(:)
         fac_epsd=(eps_d+1)/(eps_d-1)
         Kdiagd(:)=(twp-sgn*eigv(:))/(twp*fac_epsd-sgn*eigv(:))
       elseif (Feps.eq."drl") then       
!        Drude-Lorentz dielectric function
         fact1(:)=(two*pi-sgn*eigv(:))*eps_A/(two*twp)  
         fact2(:)=((two*pi-sgn*eigv(i))*eps_A/(two*twp)+eps_w0*eps_w0)  
       else
         write(6,*) "Wrong epsilon choice"
         stop
       endif
!      Buil S-1/2T in scr3
       scr3=matmul(sm12,eigt)
!      Buil T*S1/2 in scr2
       scr2=matmul(transpose(eigt),sp12)  
!      Buil T*S-1/2 in scr4
       scr4=matmul(transpose(eigt),sm12)  
!      Build matqf (R for debye)
       do i=1,nts_act
         scr1(:,i)=scr3(:,i)*fact1(i) 
       enddo
       matqf=matmul(scr1,scr2)
!      Build matqw (\tilde{Q} for debye)
       do i=1,nts_act
         scr1(:,i)=scr3(:,i)*fact2(i) 
       enddo
       matqw=-matmul(scr1,scr4) 
!      Build Q_0 
       do i=1,nts_act
         scr1(:,i)=scr3(:,i)*Kdiag0(i) 
       enddo
       matq0=-matmul(scr1,scr4) 
!      Build Q_d only for debye
       if (Feps.eq."deb") then
         do i=1,nts_act
           scr1(:,i)=scr3(:,i)*Kdiagd(i) 
         enddo
         matqd=-matmul(scr1,scr4) 
       endif
! Print matrices in output
       open(7,file="BEM_matrices.inp",status="unknown")
       write(7,*) nts_act
       do j=1,nts_act
        do i=j,nts_act
         write(7,'(4E26.16)')matqw(i,j),matqd(i,j),matqf(i,j),matq0(i,j)
        enddo
       enddo
       close(7)
!      Deallocate arrays
       if (allocated(Kdiag0)) deallocate(Kdiag0)
       if (allocated(Kdiagd)) deallocate(Kdiagd)
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
      subroutine do_calderon    
      integer(i4b):: i,j,info,lwork
      real(dbl):: temp
      ! create calderon D, D* and S matrices
      do i=1,nts_act
       do j=1,nts_act
         call green_s(i,j,temp)
         cals(i,j)=temp
         call green_d(i,j,temp)
         cald(i,j)=temp
         call green_d(j,i,temp)
         calds(i,j)=temp
       enddo
      enddo
      return
      end subroutine
!
      subroutine do_ons   
      integer(i4b) :: i,j,info,lwork,liwork
       real(dbl), dimension(nts_act,nts_act) :: scr1,scr2
      character jobz,uplo
      integer(i4b) :: iwork(3+5*nts_act)
      real(dbl) :: work(1+6*nts_act+2*nts_act*nts_act)
!
!      Set parameters for diagonalization routine dsyevd
       jobz = 'V'
       uplo = 'U'
       lwork = 1+6*nts_act+2*nts_act*nts_act
       liwork = 3+5*nts_act
!      
!      Form S^-1
       eigt = cals
       call dsyevd (jobz,uplo,nts_act,eigt,nts_act,eigv,work,lwork, &
         iwork,liwork,info)
       do i=1,nts_act
         scr1(:,i)=eigt(:,i)/eigv(i)
       enddo
       scr2=transpose(eigt)
       sm1=matmul(scr1,scr2)                   
!      
!      Form the f_w=-1/tau+w_0^2 and f_f=A/3                           
       f_w=-(eps_A/3.+eps_w0*eps_w0) 
       f_f=-eps_A/3. 
       return
      end subroutine
!
      subroutine do_eps      
       eps=dcmplx(eps_A,zero)/dcmplx(eps_w0**2-omega**2,-omega*eps_gm)
       eps=eps+onec
       eps_f=(eps-onec)/(eps+twoc)
      return
      end subroutine
!
      end module
