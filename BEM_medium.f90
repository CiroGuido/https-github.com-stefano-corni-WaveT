      Module BEM_medium      
      use readio
      use cav_types
      use pedra_friends
      use spectra
      use, intrinsic :: iso_c_binding

      implicit none
      real(dbl), parameter :: zero=0.d0
      real(dbl), parameter :: one=1.d0
      real(dbl), parameter :: two=2.d0
      integer(i4b), parameter :: one_i=1
      complex(cmp), parameter :: onec=(one,zero)                
      complex(cmp), parameter :: twoc=(two,zero)                
      real(dbl), parameter :: TOANGS=0.52917724924D+00
      real(dbl), parameter :: ANTOAU=1.0D+00/TOANGS
!LIGHT SPEED IN AU TO BE REVISED
      real(dbl), parameter :: c=1.37036d2                
      real(dbl) :: f_f,f_w    
      real(dbl), allocatable :: cals(:,:),cald(:,:),calds(:,:)
      real(dbl), allocatable :: sm1(:,:)
      real(dbl), allocatable :: matqf(:,:),matqw(:,:)
      real(dbl), allocatable :: eigv(:),eigt(:,:)
      complex(cmp) :: eps,eps_f
      save
      private
      public init_nanop,deallocate_nanop,do_eps,eps,eps_f,eigv,eigt, &
             cals,cald,calds,sm1,matqf,matqw,f_f,f_w

      contains
!     
      subroutine init_nanop   
      integer(i4b) :: its                    
       ! Build Matrices for propagation
       allocate(cals(nts_act,nts_act))
       allocate(cald(nts_act,nts_act))
       allocate(calds(nts_act,nts_act))
       allocate(eigv(nts_act))
       allocate(eigt(nts_act,nts_act))
       call do_calderon
       if(metdyn.eq."PCM") then
        allocate(matqf(nts_act,nts_act))
        allocate(matqw(nts_act,nts_act))
        call do_PCMdiag  
       else
        allocate(sm1(nts_act,nts_act))
        call do_ons  
       endif
      return
      end subroutine
!
      subroutine deallocate_nanop   
       ! Build Matrices for propagation
       deallocate(cals,cald,calds)
       if (allocated(matqf)) deallocate(matqf)
       if (allocated(matqw)) deallocate(matqw)
       if (allocated(sm1))deallocate(sm1)
       deallocate(eigt,eigv)
      return
      end subroutine
!
      subroutine do_PCMdiag   
       integer(i4b) :: i,j,info,lwork,liwork
       real(dbl), allocatable :: sp12(:,:),sm12(:,:)
       real(dbl), allocatable :: scr1(:,:),scr2(:,:)
       real(dbl), allocatable :: scr3(:,:)
       character jobz,uplo
       integer(i4b) :: iwork(3+5*nts_act),cont
       real(dbl) :: work(1+6*nts_act+2*nts_act*nts_act), test
!      
       allocate(sp12(nts_act,nts_act))
       allocate(sm12(nts_act,nts_act))
       allocate(scr1(nts_act,nts_act))
       allocate(scr2(nts_act,nts_act))
       allocate(scr3(nts_act,nts_act))
!      Set parameters for diagonalization routine dsyevd
       jobz = 'V'
       uplo = 'U'
       lwork = 1+6*nts_act+2*nts_act*nts_act
       liwork = 3+5*nts_act
!      
!      Form S^1/2 and S^-1/2
!      Copy the matrix in the eigenvector matrix
       eigt = cals
       call dsyevd (jobz,uplo,nts_act,eigt,nts_act,eigv,work,lwork, &
         iwork,liwork,info)
       do i=1,nts_act
         scr1(:,i)=eigt(:,i)*sqrt(eigv(i))
       enddo
       scr2=transpose(eigt)
       sp12=matmul(scr1,scr2)                   
       do i=1,nts_act
         scr1(:,i)=eigt(:,i)/sqrt(eigv(i))
       enddo
       sm12=matmul(scr1,scr2)                   
!      Test on S Diagonal passed 
!      
!      Form the S^-1/2 D A S^1/2 + S^1/2 A D* S^-1/2 , and diagonalize it
       !S^-1/2 D A S^1/2
       do i=1,nts_act
         scr1(:,i)=cald(:,i)*cts_act(i)%area
       enddo
       eigt=matmul(sm12,scr1)                   
       scr1=matmul(eigt,sp12)                   
       !S^1/2 A D* S^-1/2
       do i=1,nts_act
         scr2(i,:)=calds(i,:)*cts_act(i)%area
       enddo
       eigt=matmul(sp12,scr2)                   
       scr2=matmul(eigt,sm12)                   
       !S^-1/2 D A S^1/2+S^1/2 A D* S^-1/2
       do i=1,nts_act
         do j=1,nts_act
           eigt(i,j)=0.5*(scr1(i,j)+scr2(i,j))
         enddo
       enddo
       call dsyevd (jobz,uplo,nts_act,eigt,nts_act,eigv,work,lwork, &
         iwork,liwork,info)
!      Test on S^-1/2 D A S^1/2 + S^1/2 A D* S^-1/2 Diagonal passed 
!      
!      Form the Q_w and Q_f                           
       scr3=matmul(sm12,eigt)
       do i=1,nts_act
         matqw(:,i)=scr3(:,i)*(2*pi+eigv(i))*eps_A/(4*pi)
       enddo
       scr1=transpose(eigt)
       scr2=matmul(scr1,sm12)  
       matqf=-matmul(matqw,scr2)
       do i=1,nts_act
        sm12(:,i)=scr3(:,i)*((2*pi+eigv(i))*eps_A/(4*pi)+eps_w0*eps_w0) 
       enddo
       scr2=matmul(scr1,sp12)
       matqw=-matmul(sm12,scr2) 
! Print Matrices
       write(*,*) "Matrices"
       do i=1,nts_act
         do j=i,nts_act
           write(8,*) matqw(i,j),matqf(i,j)
         enddo
       enddo
       deallocate(sp12)
       deallocate(sm12)
       deallocate(scr1)
       deallocate(scr2)
       deallocate(scr3)
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
      real(dbl), allocatable :: scr1(:,:),scr2(:,:)
      character jobz,uplo
      integer(i4b) :: iwork(3+5*nts_act)
      real(dbl) :: work(1+6*nts_act+2*nts_act*nts_act), test
!
       allocate(scr1(nts_act,nts_act))
       allocate(scr2(nts_act,nts_act))
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
       deallocate(scr1)
       deallocate(scr2)
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
