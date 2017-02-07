      module propagate
      use readio
      use readio_medium
      use td_ContMed
      use cav_types      
      use pedra_friends  
      use spectra
      use random
      use dissipation

      implicit none
      real(8),allocatable :: f(:,:)
! SC mu_a is the dipole moment at current step,
!    int_rad is the classical radiated power at current step
!    int_rad_int is the integral of the classical radiated power at current step
      real(8) :: mu_a(3),int_rad,int_rad_int
      integer(4) :: file_c=10,file_e=8,file_mu=9,file_p=11 
      integer(4) :: file_dm=12, file_dp=13
      save
      private
      public create_field, prop
!
      contains
!
      subroutine prop
       implicit none
       integer(4)                :: i,j,k,istop
       complex(16), allocatable  :: c(:),c_prev(:),c_prev2(:)
       real(8),     allocatable  :: h_int(:,:), h_dis(:,:), h_rnd(:,:)
       real(8),     allocatable  :: pjump(:) 
       real(8)                   :: f_prev(3),f_prev2(3)
       real(8)                   :: mu_prev(3),mu_prev2(3),mu_prev3(3),&
                                    mu_prev4(3), mu_prev5(3)
       real(8)                   :: eps  
       real(8),     allocatable  :: w(:), w_prev(:)
       character(20)             :: name_f
       logical                   :: first=.true. 
! OPEN FILES
       write(name_f,'(a4,i0,a4)') "c_t_",n_f,".dat"
       open (file_c,file=name_f,status="unknown")
       write(name_f,'(a4,i0,a4)') "e_t_",n_f,".dat"
       open (file_e,file=name_f,status="unknown")
       write(name_f,'(a5,i0,a4)') "mu_t_",n_f,".dat"
       open (file_mu,file=name_f,status="unknown")
       write(name_f,'(a4,i0,a4)') "p_t_",n_f,".dat"
       open (file_p,file=name_f,status="unknown")
       write(name_f,'(a5,i0,a4)') "dm_t_",n_f,".dat"
       open (file_dm,file=name_f,status="unknown")
       write(name_f,'(a5,i0,a4)') "dp_t_",n_f,".dat"
       open (file_dp,file=name_f,status="unknown")
! ALLOCATING
       allocate (c(n_ci))
       allocate (c_prev(n_ci))
       allocate (c_prev2(n_ci))
       allocate (h_int(n_ci,n_ci))
       if (dis) then
          allocate (h_dis(n_ci,n_ci))
          if (qjump) then
             allocate (pjump(3*nexc+1))
          elseif (dis.and..not.qjump) then 
             allocate (h_rnd(n_ci,n_ci))
             allocate (w(n_ci), w_prev(n_ci))
          endif
       endif

! STEP ZERO: build interaction matrices to do a first evolution
       c_prev2=c_i
       c_prev=c_i
       c=c_i
       f_prev2=f(:,2)
       f_prev=f(:,1)
       h_int=zero  
       mu_prev=0.d0
       mu_prev2=0.d0
       mu_prev3=0.d0
       mu_prev4=0.d0
       mu_prev5=0.d0
       int_rad_int=0.d0
       if(rad.eq."arl".or.dis) call seed_random_number_sc(iseed)
       call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
       if (mdm.ne."vac") then
           call init_mdm(c_prev,c_prev2,f_prev,f_prev2,h_int)
           call prop_mdm(1,c_prev,c_prev2,f_prev,f_prev2,h_int)
       endif
       call add_int_vac(f_prev,h_int)

! EC 20/12/16
! Dissipation according to the Markovian SSE (eq 25 J. Phys: Condens.
! Matter vol. 24 (2012) 273201)
! OR
! Dissipation according to the non-Markovian SSE (eq 24 J. Phys:
! Condens. Matter vol. 24 (2012) 273201) 
! Add a random fluctuation for the stochastic propagation (.not.qjump)
       if (dis) then
          call define_h_dis(h_dis,n_ci,tdis)
          if (.not.qjump) call add_h_rnd(h_rnd,n_ci) 
       endif
!
! INITIAL STEP: dpsi/dt=(psi(2)-psi(1))/dt
       c=c_prev-ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))
       if (dis) then
          c=c-dt*matmul(h_dis,c_prev)
          if (.not.qjump) then
             call rnd_noise(w,w_prev,n_ci,first,tdis)
             first=.false.
             if (tdis.eq.0) then
             ! Euler-Maruyama
                c=c-ui*sqrt(dt)*matmul(h_rnd,c_prev)*w*dt
             elseif (tdis.eq.1) then
             ! Leimkuhler-Matthews
                c=c-ui*0.5d0*sqrt(dt)*matmul(h_rnd,c_prev)*(w+w_prev)
             endif
          endif
          !c_prev=c
       endif
       c=c/sqrt(dot_product(c,c))
       call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
       call out_header
       if (mod(2,n_out).eq.0) call output(2,c,f_prev,h_int)
!
! PROPAGATION CYCLE: starts the propagation at timestep 3
       do i=3,n_step
         f_prev2=f(:,i-2)
         f_prev=f(:,i-1)
         h_int=zero 
         if (mdm.ne."vac") call prop_mdm(i,c_prev,c_prev2,f_prev, & 
                                                      f_prev2,h_int)
         call add_int_vac(f_prev,h_int)
! SC field
         if (rad.eq."arl".and.i.gt.5) call add_int_rad(mu_prev,mu_prev2,mu_prev3, &
                                                mu_prev4,mu_prev5,h_int)

! No dissipation in the propagation -> .not.dis
! Markovian dissipation (quantum jump) -> dis.and.qjump
! Markovian dissipation (Euler-Maruyama) -> dis.and.not.qjump
         if (.not.dis) then
            c=c_prev2-2.d0*ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))
            c=c/sqrt(dot_product(c,c))
            c_prev2=c_prev
         elseif (dis.and.qjump) then
! Quantum jump (spontaneous or nonradiative relaxation, pure dephasing)
! Algorithm from J. Opt. Soc. Am. B. vol. 10 (1993) 524
            c=c_prev-ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))-dt*matmul(h_dis,c_prev)
! loss_norm computes: 
! norm = 1 - dtot
! dtot = dsp + dnr + dde
! Loss of the norm, dissipative events simulated
! eps -> uniform random number in [0,1]
            call loss_norm(c,n_ci,pjump)

            call random_number(eps)  
            if (dtot.gt.eps)  then
               call quan_jump(c,n_ci,pjump)
            else
               c=c/sqrt(dot_product(c,c))  
            endif
         elseif (dis.and..not.qjump) then
! Dissipation by a continuous stochastic propagation
            call rnd_noise(w,w_prev,n_ci,first,tdis)
            if (tdis.eq.0) then
            ! Euler-Maruyama 
                c=c_prev-ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))-dt*matmul(h_dis,c_prev)-ui*sqrt(dt)*matmul(h_rnd,c_prev)*w
            elseif (tdis.eq.1) then
                c=c_prev-ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))-dt*matmul(h_dis,c_prev)-ui*0.5d0*sqrt(dt)*matmul(h_rnd,c_prev)*(w+w_prev)
            endif
            c=c/sqrt(dot_product(c,c))
         endif

         !c_prev2=c_prev
         c_prev=c
         call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
         if (mod(i,n_out).eq.0)  call output(i,c,f_prev,h_int)
         
       enddo

! DEALLOCATION AND CLOSING
       deallocate (c,c_prev,c_prev2,h_int)
       if (dis) then
          call deallocate_dis()
          deallocate(h_dis)
          if (.not.qjump) then
             deallocate(h_rnd)
             deallocate(w)
             deallocate(w_prev)
          endif 
          if (qjump) then
             write(*,*)
             write(*,*) 'Total number of quantum jumps', i_sp+i_nr+i_de, '(',  real(i_sp+i_nr+i_de)/real(n_step),')'
             write(*,*) 'Spontaneous emission quantum jumps', i_sp,'(',  real(i_sp)/real(n_step),')'
             write(*,*) 'Nonradiarive relaxation quantum jumps', i_nr,'(', real(i_nr)/real(n_step),')'
             write(*,*) 'Dephasing relaxation quantum jumps', i_de, '(', real(i_de)/real(n_step),')'
             write(*,*)
          endif
       endif
       close (file_c)
       close (file_e)
       close (file_mu)
       close (file_p)
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
       character(15) :: name_f
       allocate (f(3,n_step))
       write(name_f,'(a5,i0,a4)') "field",n_f,".dat"
       open (7,file=name_f,status="unknown")
        f(:,:)=0.d0
        select case (Ffld)
        case ("mdg")
         do i=1,n_step 
          t_a=dt*(i-1)
          f(:,i)=fmax(:)*exp(-(t_a-t_mid)**2/(sigma**2))*   &
                         sin(omega*t_a)
          if (mod(i,n_out).eq.0) &
           write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
         enddo
        case ("mds")
         i_max=int(t_mid/dt)
         do i=1,2*i_max
          t_a=dt*(i-1)
          f(:,i)=fmax*cos(pi*(t_a-t_mid)/(2*t_mid))**2/2.d0* &
                 sin(omega*t_a)
          if (mod(i,n_out).eq.0) &
            write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
         enddo
         do i=2*i_max+1,n_step
          t_a=dt*(i-1)
          f(:,i)=0.
          if (mod(i,n_out).eq.0) &
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
          if (mod(i,n_out).eq.0) &
            write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
         enddo
        case ("sin")
         do i=1,n_step
          t_a=dt*(i-1)
          f(:,i)=fmax(:)*sin(omega*t_a)
          if (mod(i,n_out).eq.0) &
            write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
         enddo
        case ("snd")
         do i=1,n_step
          t_a=dt*(i-1)
          if (t_a.gt.t_mid) then
            f(:,i)=fmax(:)*sin(omega*t_a)
          else
            if (t_a.gt.0.d0) f(:,i)=fmax(:)*t_a/t_mid*sin(omega*t_a)
          endif            
          if (mod(i,n_out).eq.0) &
            write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
         enddo
        case ("gau")
         do i=1,n_step
          t_a=dt*(i-1)
          f(:,i)=fmax(:)*exp(-(t_a-t_mid)**2/(sigma**2))
          if (mod(i,n_out).eq.0) &
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
          if (mod(i,n_out).eq.0) &
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
      subroutine do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
       implicit none
       complex(16), intent(IN) :: c(n_ci)
       real(8)::mu_prev(3),mu_prev2(3),mu_prev3(3),mu_prev4(3), &
                mu_prev5(3)
       mu_a(1)=dot_product(c,matmul(mut(:,:,1),c))
       mu_a(2)=dot_product(c,matmul(mut(:,:,2),c))
       mu_a(3)=dot_product(c,matmul(mut(:,:,3),c))
! SC save previous mu for radiative damping
       mu_prev5=mu_prev4
       mu_prev4=mu_prev3
       mu_prev3=mu_prev2
       mu_prev2=mu_prev
       mu_prev=mu_a
      return
      end subroutine

      subroutine output(i,c,f_prev,h_int)     
       implicit none
       integer(4), intent(IN) :: i
       complex(16), intent(IN) :: c(n_ci)
       real(8), intent(IN) :: h_int(n_ci,n_ci)
       real(8), intent(IN) :: f_prev(3)
       real(8) :: e_a,e_vac,t,g_neq_t,g_neq2_t,g_eq_t,f_med(3)
       character(20) :: fmt_ci,fmt_ci1
        t=(i-1)*dt
        e_a=dot_product(c,e_ci*c+matmul(h_int,c))
! SC 07/02/16: added printing of g_neq, g_eq 
        if(mdm.ne.'vac') then 
         g_eq_t=e_a
         g_neq_t=e_a
         g_neq2_t=e_a
         e_vac=e_a
         call get_gneq(e_vac,g_eq_t,g_neq_t,g_neq2_t)
         write (file_e,'(i8,f14.4,7e20.8)') i,t,e_a,e_vac, &
                   g_eq_t,g_neq2_t,g_neq_t,int_rad,int_rad_int
        else
         write (file_e,'(i8,f14.4,3e22.10)') i,t,e_a,int_rad,int_rad_int
        endif
        write (fmt_ci,'("(i8,f14.4,",I0,"e13.5)")') n_ci
        write (fmt_ci1,'("(i8,f14.4,",I0,"e13.5)")') n_ci*(n_ci-1)/2
        write (file_c,fmt_ci) i,t,real(c(:)*conjg(c(:)))
        write (file_p,fmt_ci) i,t,atan2(aimag(c),real(c))
        call wrt_decoherence(i,t,file_dm,file_dp,fmt_ci1,c,n_ci)
        write (file_mu,'(i8,f14.4,3e22.10)') i,t,mu_a(:)
        Sdip(i,:,1)=mu_a(:)
        !Sfld(i,:)=f(:,i-1)
        ! SP 25/10/16: changed for FT
        Sfld(i,:)=f(:,i)
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
      subroutine add_int_rad(mu_prev,mu_prev2,mu_prev3,mu_prev4, &
                                                   mu_prev5,h_int)
! SC calculate the Aharonov Lorentz radiative damping
       implicit none
       real(8),intent(in) :: mu_prev(3),mu_prev2(3),mu_prev3(3), &
                 mu_prev4(3), mu_prev5(3)
       real(8), intent(INOUT) :: h_int(n_ci,n_ci)
       real(8) :: d3_mu(3),d2_mu(3),d_mu(3),d2_mod_mu,coeff,scoeff, &
                  force(3),de
!       d3_mu=(mu_prev-3.*mu_prev2+3.*mu_prev3-mu_prev4)/(dt*dt*dt)
       d3_mu=(2.5*mu_prev-9.*mu_prev2+12.*mu_prev3-7.*mu_prev4+ &
              1.5*mu_prev5)/(dt*dt*dt)
! SC 08/06/2016: i'm confused on the right formula for istantaneous radiated intensity:
!  Novotny-Hecht is pro. to (d^2/dt^2 |mu|)^2 (eq. 8.70), however in Jackson
! for Larmor one has pro. to (d^2/dt^2 mu \cdot d^2/dt^2 mu). Since Novotny-Hech
! in eq. 8.82 seems to use jacson definition, I also use it here
!       d2_mod_mu=(2.*sqrt(dot_product(mu_prev,mu_prev)) &
!                  -5.*sqrt(dot_product(mu_prev2,mu_prev2))+ &
!                  4.*sqrt(dot_product(mu_prev3,mu_prev3)) &
!               -sqrt(dot_product(mu_prev4,mu_prev4)))/(dt*dt)
       d2_mu=(2.*mu_prev-5.*mu_prev2+4.*mu_prev3-mu_prev4)/(dt*dt)
       d_mu=(1.5*mu_prev-2.*mu_prev2+0.5*mu_prev3)/dt
!SC coefficient 1/(6 pi eps0 c^3) in atomic units
       coeff=2.d0/3.d0/137.036**3.
       scoeff=coeff
!SC: added random force             
!       d3_mu=d3_mu*coeff
       force(1)=random_normal()*scoeff
       force(2)=random_normal()*scoeff
       force(3)=random_normal()*scoeff
       d3_mu=d3_mu*coeff
       de=dot_product(d_mu,force)
       if(de.lt.0) d3_mu=d3_mu+force
!       write (6,*) d3_mu
!SC Instantenous emitted intensity (from Novotny Hech eq. 8.70)
!       int_rad=coeff*d2_mod_mu*d2_mod_mu
       int_rad=coeff*dot_product(d2_mu,d2_mu)+de/dt
       int_rad_int=int_rad_int+int_rad*dt
       h_int(:,:)=h_int(:,:)-mut(:,:,1)*d3_mu(1)-             &
                   mut(:,:,2)*d3_mu(2)-mut(:,:,3)*d3_mu(3)
       return
      end subroutine
!
      subroutine out_header
! SC write headers to output files, to be completed!
      implicit none
      write(file_e,'(8a)') '#   istep time',' <H(t)>-E_gs(0)', &
              ' DE_vac(t)',' DG_eq(t)',' DG_neq(t)',  '  Const', &
              '  Rad. Int', '  Rad. Ene'   
      return
      end subroutine
!
      subroutine wrt_decoherence(i,t,int1,int2,char20,c,nci)
!------------------------------------------------------------------------
! Print the tridiagional C*_iC_j (i.ne.j) matrix 
! corresponding to the decoherence 
! 
! Created   : E. Coccia 3 Feb 2017
! Modified  :
!------------------------------------------------------------------------
  
        implicit none

        integer(4),    intent(in)    :: i, nci
        integer(4),    intent(in)    :: int1, int2
        character(20), intent(in)    :: char20 
        real(8),       intent(in)    :: t
        complex(16),   intent(in)    :: c(nci)
        integer(4)                   :: j, k
        complex(16)                  :: tmp
        real(8)                      :: r(nci,nci*(nci-1)/2)
        real(8)                      :: p(nci,nci*(nci-1)/2)

        do j=2,nci
           tmp = conjg(c(j))*c(1)
           r(1,j) = abs(tmp)
           p(1,j) = atan2(aimag(tmp),real(tmp))
        enddo  

        do j=2,nci
           do k=j+1,nci
              tmp = conjg(c(k))*c(j)
              r(j,k) = abs(tmp)
              p(j,k) = atan2(aimag(tmp),real(tmp)) 
           enddo
         enddo

         write(int1,char20) i, t, ((r(j,k), k=j+1,nci), j=1,nci)
         write(int2,char20) i, t, ((p(j,k), k=j+1,nci), j=1,nci)    
         !((FORM(K,L), L=1,10), K=1,10,2)
    
         return

      end subroutine wrt_decoherence

      end module
