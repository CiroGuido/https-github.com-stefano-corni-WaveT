      Module QM_coupling    
      use cav_types
      use readio
      use readio_medium
      use pedra_friends
      use spectra
      use BEM_medium
      use, intrinsic :: iso_c_binding

      implicit none
      real(dbl), parameter :: TOANGS=0.52917724924D+00
      real(dbl), parameter :: ANTOAU=1.0D+00/TOANGS
!LIGHT SPEED IN AU TO BE REVISED
      real(dbl), allocatable :: omegax(:),Htot(:,:)                
      real(dbl), parameter :: c=1.37036d2                
      real(dbl), allocatable :: Meigt(:,:)
      real(dbl), allocatable :: Meigv(:)
      save
      private
      public do_diag  
!
      contains
!     
      subroutine do_modes                     
         !double complex, allocatable :: Kdiag(:)                
       call do_eps
       ! obtinaed by solving Re(2pi*(eps+1)/(eps-1)+lambda_i)=0
         !allocate(Kdiag(nts_act))
         !Kdiag(:)=cmplx(2*pi+eigv(:),zero)/                      &
         !          (2*pi*eps_f+cmplx(eigv(:),zero)) 
       omegax(:)=sqrt(eps_w0**2+(2*pi+eigv(:))*eps_A/(4*pi))     
         !deallocate(Kdiag)
       ! omegax(2) Tested against the spherical plasmon
      return
      end subroutine
!
      subroutine do_matrix                       
       real(dbl):: omegai,gFi,dp(nts_act)                   
       integer(4)::i,j,k,p,s     
       real(dbl), allocatable :: tsm(:,:)                
       ! Build Hamiltonian Super-matrix
       ! Build the diagonal superblocs:
       ! H11
       allocate(tsm(nts_act,nts_act))
       do j=1,n_ci
         do k=1,n_ci
           Htot(k,j)=-dot_product(mut(k,j,:),fmax(:))
         enddo
         Htot(j,j)=Htot(j,j)+e_ci(j)
       enddo
       tsm=matmul(transpose(eigt),sm12)
       do i=1,nts_act
         dp(i)=cts_act(i)%x*fmax(1)+cts_act(i)%y*fmax(2)+       &
                                    cts_act(i)%z*fmax(3) 
       enddo 
       do i=2,nmodes+1
         omegai=eps_w0**2+eps_A*(pt5+pt5*pt5*eigv(xmode(i-1))/pi)
         gFi=-dot_product(tsm(xmode(i-1),:),dp(:))* &
                       sqrt((omegai**2-eps_w0**2)/(two*omegai))
         do j=1,n_ci
           p=(i-1)*nmodes+j
           do k=j,n_ci
             s=(i-1)*nmodes+k
             ! H22 
             Htot(s,p)=Htot(k,j)
             Htot(p,s)=Htot(s,p)
             ! H12 no check the indices
             Htot(k,p)=dot_product(tsm(xmode(i-1),:),vts(k,j,:))* &
                       sqrt((omegai**2-eps_w0**2)/(two*omegai))
             Htot(j,s)=Htot(k,p)
             ! H21
             Htot(s,j)=Htot(k,p)
             Htot(p,k)=Htot(s,j)
           enddo
           Htot(p,p)=Htot(p,p)+omegax(i-1)*(occmd(i-1)+pt5)    
           Htot(p,j)=Htot(p,j)+gFi
           Htot(j,p)=Htot(j,p)+gFi
         enddo
       enddo
       deallocate(tsm)
      return
      end subroutine
       
      subroutine diag_matrix(Mdim)                       
       integer(i4b),intent(in) :: Mdim   
       integer(i4b) :: i,j,info,lwork,liwork
       character jobz,uplo
       integer(i4b) :: iwork(3+5*Mdim)
       real(dbl) :: work(1+6*Mdim+2*Mdim*Mdim)
!      Set parameters for diagonalization routine dsyevd
       jobz = 'V'
       uplo = 'U'
       lwork = 1+6*Mdim+2*Mdim*Mdim
       liwork = 3+5*Mdim
       eigt(:,:) = Htot(:,:)
       call dsyevd (jobz,uplo,Mdim,Meigt,Mdim,Meigv,work,lwork, &
         iwork,liwork,info)
      return
      end subroutine
!
      subroutine do_diag                         
       integer(i4b) :: Mdim,i   
       ! Build the diagonal kernel matrix and find poles                    
       allocate(omegax(nts_act))
       call do_modes 
       ! Build Hamiltonian Super-matrix
       Mdim=n_ci*nmodes
       allocate(Htot(Mdim,Mdim))
       call do_matrix
       ! Diagonalize Super-matrix        
       allocate(Meigt(Mdim,Mdim))
       allocate(Meigv(Mdim))
       call diag_matrix(Mdim)
       do i=1,Mdim
         if(Mdim.le.n_ci) then
           write(*,*) i, Meigv(i), e_ci(i)
         else
           write(*,*) i, Meigv(i)
         endif
       enddo
       deallocate(Meigt)
       deallocate(Meigv)
       deallocate(Htot)
       deallocate(omegax)
      return
      end subroutine
!
      end module
