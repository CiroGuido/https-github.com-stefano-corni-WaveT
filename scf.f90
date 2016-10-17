      Module scf            
      use cav_types
      use readio
      use readio_medium
      use pedra_friends
      use BEM_medium    
      use, intrinsic :: iso_c_binding

      implicit none
      real(8), parameter :: ev_to_au=0.0367493
      real(dbl), parameter :: TOANGS=0.52917724924D+00
      real(dbl), parameter :: ANTOAU=1.0D+00/TOANGS
!LIGHT SPEED IN AU TO BE REVISED
      real(dbl), allocatable :: omegax(:),Htot(:,:)                
      real(dbl), parameter :: c=1.37036d2                
! SP quantities at c-th cycle 
      real(dbl), allocatable :: eigt_c(:,:),eigv_c(:)
      real(dbl) :: fc(3)
! SP quantities at previous cycle 
      real(dbl), allocatable :: eigt_cp(:,:),eigv_cp(:)
! SP Max values for eigenval/vec differences wrt previous cycle
      real(dbl) :: maxe,maxv                    
      save
      private
      public do_scf
!
      contains
!
      subroutine do_scf(qst,c_prev)                  
       implicit none
       real(dbl), intent(INOUT):: qst(:)     
       complex(16), intent(INOUT) :: c_prev(:)
       integer(i4b) :: ncyc=1 
       integer(i4b) :: its,i,j
       logical :: docycle=.true. 
       real(dbl) :: thre,thrv                    
       real(dbl) :: e_scf, e_ini
       integer(i4b):: max_p(1)    
       write(6,*) "SCF Cycle"
       write(6,*) "Max_cycle ", ncycmax
       thrv=10**(-thrshld+2)
       thre=10**(-thrshld)
       write(6,*) "Threshold ", thrv,thre
       allocate(eigv_c(n_ci),eigt_c(n_ci,n_ci))
       allocate(eigv_cp(n_ci),eigt_cp(n_ci,n_ci))
       allocate(Htot(n_ci,n_ci))
       ! scf cycle
       do while (docycle.and.ncyc.le.ncycmax) 
         call do_matrix(qst)
         call diag_matrix(n_ci,Htot)       
         call do_charges(qst)
         call do_energies(e_scf,e_ini)
         if (ncyc.gt.2) then 
           call check_conv(maxe,maxv,n_ci)       
!           if (maxe.le.thre.and.maxv.le.thrv) docycle=.false.         
! SC 24/4/2016: check convergence only on egeinvalues:
!               in case of degeneracy the variation of eigenvector
!               can be erratic
           if (maxe.le.thre) docycle=.false.         
           write(6,*) "cycle ", ncyc, e_scf, e_ini
         endif
         eigt_cp=eigt_c
         eigv_cp=eigv_c
         ncyc=ncyc+1 
         write(6,*) "Max Diff on Eigenvalue ", maxe
         write(6,*) "Max Diff on Eigenvector ", maxv
       enddo
       write(6,*) "SCF Done"
! write out the charges in the charge0.dat file
       if (Fprop.ne.'dip') then
       open(unit=7,file="charges0_scf.inp",status="unknown", &
            form="formatted")
         write (7,*) nts_act
         do its=1,nts_act
          write (7,'(E22.8,F22.10)') qst(its)
         enddo
       close(unit=7)
       write(6,*) "Written out the SCF charges"
       endif
!       write(6,*) "Max Diff on Eigenvalue ", maxe
!       write(6,*) "Max Diff on Eigenvector ", maxv
! Transform the potentials to the new basis
! SP 30/05/16: added the condition on Fint for Fint=ons and Fprop=ief
       if (Fprop.ne.'dip') then
         do its=1,nts_act
          vts(:,:,its)=matmul(vts(:,:,its),eigt_c)
          vts(:,:,its)=matmul(transpose(eigt_c),vts(:,:,its))
         enddo
         open(unit=7,file="ci_pot_scf.inp",status="unknown", &
            form="formatted")
         write (7,*) nts_act
         write (7,*) "V0  check Vnuc"
         do its=1,nts_act
          write (7,*) vts(1,1,its)-vtsn(its),0.d0,vtsn(its)
         enddo
         do j=2,n_ci
           write(7,*) 0,j-1
           do its=1,nts_act
            write(7,*) vts(1,j,its)
           enddo
         enddo
         !Vij
         do i=2,n_ci
          do j=2,i-1   
           write(7,*) i-1,j-1
           do its=1,nts_act
            write(7,*) vts(i,j,its)             
           enddo
          enddo
           write(7,*) i-1,i-1
           do its=1,nts_act
            write(7,*) vts(i,i,its)-vtsn(its)             
           enddo
         enddo
         close(unit=7)
         write(6,*) "Written out the SCF potentials"
       endif
! Transform the dipoles to the new basis
       do its=1,3
        mut(:,:,its)=matmul(mut(:,:,its),eigt_c)
        mut(:,:,its)=matmul(transpose(eigt_c),mut(:,:,its))
       enddo
! write out the scf dipoles
       open(unit=7,file="ci_mut_scf.inp",status="unknown", &
           form="formatted")
       do i=1,n_ci
        write(7,'(4i4,3f15.6)') 0,0,0,0,mut(1,i,1),mut(1,i,2),mut(1,i,3)
       enddo
       do i=2,n_ci
         do j=2,i   
           write(7,'(4i4,3f15.6)') 0,0,0,0,mut(i,j,1),mut(i,j,2),mut(i,j,3)
         enddo
       enddo
       close(unit=7)
       write(6,*) "Written out the SCF dipoles"
! Put the new diagonal energy of the hamiltonian
       e_ci=eigv_c
       open(unit=7,file="ci_energy_scf.inp",status="unknown", &
           form="formatted")
       do i=2,n_ci
        write(7,'(3i4,f15.8)') 0,0,0,(e_ci(i)-e_ci(1))/ev_to_au
       enddo
       close(unit=7)
       write(6,*) "Written out the SCF energies,", &
               " GS has been given zero energy!"
!  find the new eigenvector that is most similar to the old one
       c_i=abs(matmul(c_i,eigt_c))
       max_p=maxloc(c_i)
       write(6,*) 'maxloc',max_p(1)
       c_i=0.d0
       c_i(max_p(1))=1.d0
       c_prev=c_i
       deallocate(eigv_c,eigt_c)
       deallocate(eigv_cp,eigt_cp)
       deallocate(Htot)
      return
      end subroutine
!     
      subroutine do_charges(qts)                      
      real(dbl), intent(OUT):: qts(:)     
      integer(i4b)::i     
      integer(i4b):: max_p(1)    
      real(dbl),allocatable :: pot(:)
      real(dbl) :: c_c(n_ci)
!  find the new eigenvector that is most similar to the old one
      c_c=abs(matmul(c_i,eigt_c))
      max_p=maxloc(c_c)
      write(6,*) 'maxloc',max_p(1)
      c_c=0.d0
      c_c(max_p(1))=1.d0
      ! This is the new state on the basis of the old states  
      c_c=matmul(eigt_c,c_c)
!     write the state
      write(6,*) "State on the basis of original states"
      do i=1,n_ci
       write(6,*) i, c_c(i)
      enddo
      write(6,*)
      select case (Fprop)
      case('dip')
! SC pot is the dipole and qts is the field
       allocate(pot(3))
       do i=1,3
         pot(i)=dot_product(c_c,matmul(mut(:,:,i),c_c))
       enddo
       qts=(1.-mix_coef)*qts+mix_coef*f_0*pot
       write(6,*) pot(3),qts(3),f_0*pot(3)
      case default
       allocate(pot(nts_act))
       do i=1,nts_act
         pot(i)=dot_product(c_c,matmul(vts(:,:,i),c_c))
       enddo 
       qts=(1.-mix_coef)*qts+mix_coef*matmul(matq0,pot)
! SC 12/8/2016: apparently for NP, charge compensation is needed
       if (mdm.eq.'nan') qts=qts-sum(qts)/nts_act
      end select
      if (allocated(pot)) deallocate(pot)
      return
      end subroutine
!       
      subroutine do_matrix(qts)                      
      real(dbl), intent(IN):: qts(:)     
      integer(4)::i,j,k     
      write(6,*) "Htot(j,j)"
      select case (Fint)
      case('dip')
       do j=1,n_ci
         do k=1,j   
           Htot(k,j)=-dot_product(mut(k,j,:),qts(:)-fr0(:))
           Htot(j,k)=Htot(k,j)
         enddo
         Htot(j,j)=Htot(j,j)+e_ci(j)
         write(6,*) j,Htot(j,j)
       enddo
      case default
       do j=1,n_ci
         do k=1,j   
           Htot(k,j)=dot_product(vts(k,j,:),qts(:)-q0(:))
           Htot(j,k)=Htot(k,j)
         enddo
         Htot(j,j)=Htot(j,j)+e_ci(j)
         write(6,*) j,Htot(j,j)
       enddo
      end select
      return
      end subroutine
!       
      subroutine diag_matrix(Mdim,Mtrx)                       
       integer(i4b),intent(in) :: Mdim   
       real(dbl),intent(in) :: Mtrx(Mdim,Mdim)
       integer(i4b) :: i,j,info,lwork,liwork
       character jobz,uplo
       integer(i4b) :: iwork(3+5*Mdim)
       real(dbl) :: work(1+6*Mdim+2*Mdim*Mdim)
!      Set parameters for diagonalization routine dsyevd
       jobz = 'V'
       uplo = 'U'
       lwork = 1+6*Mdim+2*Mdim*Mdim
       liwork = 3+5*Mdim
       eigt_c(:,:)=Mtrx(:,:)
       call dsyevd (jobz,uplo,Mdim,eigt_c,Mdim,eigv_c,work,lwork, &
         iwork,liwork,info)
       write (6,*) "eigv_c"
       do i=1,mdim
         write (6,*) i,eigv_c(i)
       enddo
!        mtrx(:,:)=matmul(mtrx(:,:),eigt_c)
!        mtrx(:,:)=matmul(transpose(eigt_c),mtrx(:,:))
!       write (6,*) "htot dopo diagonalizzazione"
!       do i=1,mdim
!        do j=1,mdim
!         write (6,*) i,j,mtrx(i,j)
!        enddo
!       enddo
      return
      end subroutine
!
      subroutine check_conv(mxe,mxv,Mdim)                       
       real(dbl),intent(inout) :: mxe,mxv
       integer(i4b),intent(in) :: Mdim 
       integer(i4b) :: i,j
       real(dbl):: diff               
       mxv=zero                
       mxe=zero          
!       write(6,*)
       do i=1,Mdim   
         diff=sqrt((eigv_c(i)-eigv_cp(i))**2)
         if(diff.gt.mxe) mxe=diff
!         write (6,*) i,eigt_c(i,:)
         do j=1,Mdim   
! SC 24/4/2016: changed below, otherwise a change of sign result in non-convergence
           diff=abs(eigt_c(j,i)**2-eigt_cp(j,i)**2)
           if(diff.gt.mxv) mxv=diff
         enddo
       enddo
      return
      end subroutine
!
      subroutine do_energies(e_scf,e_ini)
      implicit none
      integer(4) :: i
      real(8) :: e_scf,e_ini
      e_scf=0.d0
      e_ini=0.d0
      do i=1,n_ci
       e_scf=abs(c_i(i))*abs(c_i(i))*eigv_c(i)+e_scf
       e_ini=abs(c_i(i))*abs(c_i(i))*e_ci(i)+e_ini
      enddo
      return
      end subroutine
!
!      subroutine do_rediag(nts,tdm)                  
!       implicit none
!       integer(i4b), intent(IN):: nts,tdm     
!       integer(i4b) :: its,i
!       allocate(eigv_c(tdm),eigt_c(tdm,tdm))
!       allocate(eigv_cp(tdm),eigt_cp(tdm,tdm))
!       allocate(Htot(tdm,tdm))
!       call do_matrix(nts,q0,tdm)
!       call diag_matrix(tdm,Htot)       
!       do its=1,nts
!        vts(:,:,its)=matmul(vts(:,:,its),eigt_c)
!        vts(:,:,its)=matmul(transpose(eigt_c),vts(:,:,its))
!       enddo
!       write (6,*) "Old vs New eigenenergies after rediagonalization"
!       do i=1,tdm
!        write (6,*) i,e_ci(i),eigv_c(i)
!       enddo
!! Put the new diagonal energy of the hamiltonian
!       e_ci=eigv_c
!       deallocate(eigv_c,eigt_c)
!       deallocate(eigv_cp,eigt_cp)
!       deallocate(Htot)
!      return
!      end subroutine
!
      end module
