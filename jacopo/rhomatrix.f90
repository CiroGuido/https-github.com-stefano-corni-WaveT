module rhomatrix
implicit none
   !Public
public :: rhocalc,vecrho
   !
   !
  ! type rhoprop
  !   complex(16),allocatable,dimension(:,:) :: mat
  !   complex(16),allocatable,dimension(:)   :: vec
  ! end type rhoprop
  ! !
  ! type(rhoprop), save :: rho
   !
   !Density matrix type
   !
contains
!  
 subroutine rhocalc(nci,coef,rhomat)
  implicit none
  integer(4) :: i,j
  integer(4), intent (IN) :: nci
  complex(16), intent (IN):: coef(nci)  
  complex(16) :: rhomat(nci,nci)
  !allocate(rho%mat(nci,nci))
   do i=1,nci
   rhomat(i,i)=coef(i)*conjg(coef(i))
     do j=i+1,nci
      rhomat(i,j)=coef(i)*conjg(coef(j))
      rhomat(j,i)=conjg(rhomat(i,j))
     enddo
   enddo
  !deallocate(rho%mat) 
 end subroutine                           

 subroutine vecrho(nci,rhomat,rhovec)
 implicit none
 integer(4) :: i, j, k
 integer(4),intent(IN) ::nci
 complex(16), intent(IN)::rhomat(nci,nci)
 complex(16)::rhovec((nci+1)*nci/2)

 !allocate(rho%mat(nci,nci))
 !allocate(rho%vec((nci+1)*nci/2))
   k=0
 do i=1, nci
   do j=i, nci !stora solo la upper part
   k=k+1
   rhovec(k) = rhomat(j,i)
   enddo
 enddo
 !deallocate(rho%mat,rho%vec)
 end subroutine

end module                 
