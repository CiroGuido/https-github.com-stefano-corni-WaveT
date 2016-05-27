      module propagate
      use readio
      use readio_medium
      use td_ContMed
      use cav_types      
      use pedra_friends  
      use spectra
      implicit none
      real(8),allocatable :: f(:,:)
      real(8) :: mu_a(3)
      integer(4) :: file_c=7,file_e=8,file_mu=9
      save
      private
      public create_field, prop
!
      contains
!
      subroutine prop
       implicit none
       integer(4) :: i,j
       complex(16), allocatable :: c(:),c_prev(:),c_prev2(:)
       real(8), allocatable :: h_int(:,:)
       real(8) :: f_prev(3),f_prev2(3)
! OPEN FILES
       open (file_c,file="c_t.dat",status="unknown")
       open (file_e,file="e_t.dat",status="unknown")
       open (file_mu,file="mu_t.dat",status="unknown")
! ALLOCATING
       allocate (c(n_ci))
       allocate (c_prev(n_ci))
       allocate (c_prev2(n_ci))
       allocate (h_int(n_ci,n_ci))

! STEP ZERO: build interaction matrices to do a first evolution
       c_prev2=c_i
       c_prev=c_i
       c=c_i
       f_prev2=f(:,2)
       f_prev=f(:,1)
       h_int=zero  
       if (mdm.ne."vac") then
           call init_mdm(c_prev,c_prev2,f_prev,f_prev2,h_int)
           call prop_mdm(1,c_prev,c_prev2,f_prev,f_prev2,h_int)
       endif
       call add_int_vac(f_prev,h_int)
!
! INITIAL STEP: dpsi/dt=(psi(2)-psi(1))/dt
       c=c_prev+ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))
       c=c/sqrt(dot_product(c,c))
       c_prev=c
       call out_header
       call output(2,c,f_prev,h_int)
!
! PROPAGATION CYCLE: starts the propagation at timestep 3
       do i=3,n_step
         f_prev2=f(:,i-2)
         f_prev=f(:,i-1)
         h_int=zero   
         if (mdm.ne."vac") call prop_mdm(i,c_prev,c_prev2,f_prev, & 
                                                      f_prev2,h_int)
         call add_int_vac(f_prev,h_int)
         c=c_prev2+2.0*ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))
         c=c/sqrt(dot_product(c,c))
         c_prev2=c_prev
         c_prev=c
         call output(i,c,f_prev,h_int)
       enddo

! DEALLOCATION AND CLOSING
       deallocate (c,c_prev,c_prev2,h_int)
       close (file_c)
       close (file_e)
       close (file_mu)
       if(mdm.ne.'vac') then 
         call  end_mdm
       endif
       return
      end subroutine
!
      subroutine create_field
       implicit none
       integer(4) :: i,i_max
       real(8) :: t_a,ti,tf
       allocate (f(3,n_step))
       open (7,file="field.dat",status="unknown")
       select case (tfield)
        case ("mdg")
         do i=1,n_step 
          t_a=dt*(i-1)
          f(:,i)=fmax(:)*exp(-(t_a-t_mid)**2/(sigma**2))*   &
                         sin(omega*t_a)
          write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
         enddo
        case ("mds")
         i_max=int(t_mid/dt)
         do i=1,2*i_max
          t_a=dt*(i-1)
          f(:,i)=fmax*cos(pi*(t_a-t_mid)/(2*t_mid))**2/2.d0* &
                 sin(omega*t_a)
          write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
         enddo
         do i=2*i_max+1,n_step
          t_a=dt*(i-1)
          f(:,i)=0.
          write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
         enddo
        case ("pip")
         i_max=int(t_mid/dt)
         do i=1,n_step
          t_a=dt*(i-1)
          f(:,i)=0.
          if (abs(t_a-t_mid).lt.sigma) then
            f(:,i)=fmax(:)*(cos(pi*(t_a-t_mid)/(2*sigma)))**2* &
               cos(omega*(t_a-t_mid))
          endif
          write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
         enddo
        case ("sin")
         do i=1,n_step
          t_a=dt*(i-1)
          f(:,i)=fmax(:)*sin(omega*t_a)
          write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
         enddo
        case ("gau")
         do i=1,n_step
          t_a=dt*(i-1)
          f(:,i)=fmax(:)*exp(-(t_a-t_mid)**2/(sigma**2))
          write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
         enddo
        case ("css")
         ti=t_mid-sigma/two
         tf=t_mid+sigma/two
         do i=1,n_step
          t_a=dt*(i-1)
          f(:,i)=zero
          if (t_a.gt.ti.and.t_a.le.tf) then
            f(:,i)=fmax(:)*(cos(pi*(t_a-t_mid)/(sigma)))**2
          endif
          write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
         enddo
        case default
         write(*,*)  "Error: wrong field type !"
         stop
        end select
        close(7)
       return
      end subroutine
!
      subroutine output(i,c,f_prev,h_int)     
       implicit none
       integer(4), intent(IN) :: i
       complex(16), intent(IN) :: c(n_ci)
       real(8), intent(IN) :: h_int(n_ci,n_ci)
       real(8), intent(IN) :: f_prev(3)
       real(8) :: e_a,e_vac,t,g_neq_t,g_neq2_t,g_eq_t,f_med(3)
        t=(i-1)*dt
        mu_a(1)=dot_product(c,matmul(mut(:,:,1),c))
        mu_a(2)=dot_product(c,matmul(mut(:,:,2),c))
        mu_a(3)=dot_product(c,matmul(mut(:,:,3),c))
        e_a=dot_product(c,e_ci*c+matmul(h_int,c))
! SC 07/02/16: added printing of g_neq, g_eq 
        if(mdm.ne.'vac') then 
         g_eq_t=e_a
         g_neq_t=e_a
         g_neq2_t=e_a
         e_vac=e_a
         call get_gneq(e_vac,g_eq_t,g_neq_t,g_neq2_t)
         write (file_e,'(i8,f14.4,5e20.8)') i,t,e_a,e_vac, &
                   g_eq_t,g_neq2_t,g_neq_t
        else
         write (file_e,'(i8,f14.4,e22.10)') i,t,e_a
        endif
        write (file_c,*) i,t,real(c(:)*conjg(c(:)))
        write (file_mu,'(i8,f14.4,3e22.10)') i,t,mu_a(:)
        Sdip(i,:,1)=mu_a(:)
        Sfld(i,:)=f(:,i-1)
       return
      end subroutine
!
      subroutine add_int_vac(f_prev,h_int)
       implicit none
       real(8), intent(IN) :: f_prev(3)
       real(8), intent(INOUT) :: h_int(n_ci,n_ci)
       ! create the field term of the hamiltonian
! SC 16/02/2016: changed to - sign, 
       h_int(:,:)=h_int(:,:)-mut(:,:,1)*f_prev(1)-             &
                   mut(:,:,2)*f_prev(2)-mut(:,:,3)*f_prev(3)
      return
      end subroutine
!
      subroutine out_header
! SC write headers to output files, to be completed!
      implicit none
      write(file_e,*) '#   istep time   <H(t)>-E_gs(0)  DE_vac(t) DG_eq(t)    DG_neq(t)    Const'   
      return
      end subroutine

!
      end module
