      module propagate
      use constants   
      use global_wavet
      use readio
      !use pedra_friends  
      use spectra
      use random
      use dissipation

      implicit none
      real(dbl),allocatable :: f(:,:)
! SC mu_a is the dipole moment at current step,
!    int_rad is the classical radiated power at current step
!    int_rad_int is the integral of the classical radiated power at current step
      real(dbl) :: mu_a(3),int_rad,int_rad_int
      integer(i4b) :: file_c=10,file_e=8,file_mu=9,file_p=15 
      integer(i4b) :: file_dm=12, file_dp=13, file_d=14
      save
      private
      public create_field, prop
!
      contains
!
      subroutine prop
!------------------------------------------------------------------------
! @brief Propogate C(t) using a second
! order Euler algorithm 
! 
! 
! @date Created   : 
! Modified  : E. Coccia Dec-Apr 2017
!------------------------------------------------------------------------

       implicit none
       integer(i4b)                :: i,j,k,istop,ijump=0
       complex(cmp), allocatable  :: c(:),c_prev(:),c_prev2(:), h_rnd(:,:), h_rnd2(:,:)
       real(dbl),     allocatable  :: h_int(:,:), h_dis(:,:)
       real(dbl),     allocatable  :: pjump(:) 
       real(dbl)                   :: f_prev(3),f_prev2(3)
       real(dbl)                   :: mu_prev(3),mu_prev2(3),mu_prev3(3),&
                                    mu_prev4(3), mu_prev5(3)
       real(dbl)                   :: eps  
       real(dbl),     allocatable  :: w(:), w_prev(:)
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
       write(name_f,'(a4,i0,a4)') "d_t_",n_f,".dat"
       open (file_d,file=name_f,status="unknown")
! ALLOCATING
       allocate (c(n_ci))
       allocate (c_prev(n_ci))
       allocate (c_prev2(n_ci))
       allocate (h_int(n_ci,n_ci))
! SP 17/07/17: new flags
       !if (dis) then
       if (Fdis(1:3).eq."mar".or.Fdis(1:3).eq."nma") then
          allocate (h_dis(n_ci,n_ci))
          !if (qjump) then
          if (Fdis(5:9).eq."qjump") then
             allocate (pjump(2*nf+nexc+1))
          !elseif (dis.and..not.qjump) then 
          else
             allocate (h_rnd(n_ci,n_ci))
             allocate (h_rnd2(n_ci,n_ci))
             allocate (w(3*n_ci), w_prev(3*n_ci))
          endif
       endif

! STEP ZERO: build interaction matrices to do a first evolution
       c_prev2=c_i
       c_prev=c_i
       c=c_i
! SP: 17/07/17: changed the following f_prev2 <= f_prev
       !f_prev2=f(:,2)
       f_prev2=f(:,1)
       f_prev=f(:,1)
       h_int=zero  
       mu_prev=0.d0
       mu_prev2=0.d0
       mu_prev3=0.d0
       mu_prev4=0.d0
       mu_prev5=0.d0
       int_rad_int=0.d0
       !if(rad.eq."arl".or.dis.or.ernd) call seed_random_number_sc(iseed)
       if(Frad.eq."arl".or.Fdis.ne."nodis") & 
                               call seed_random_number_sc(iseed)
       call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
       i=1
       if (Fmdm(1:3).ne."vac") then
           call init_medium(c_prev,f_prev,h_int)
           call prop_medium(i,c_prev,f_prev,h_int)
       endif
       call add_int_vac(f_prev,h_int)

! SP 16/07/17: added call to output at step 0 to have full output in outfiles
       call out_header
       call output(1,c,f_prev,h_int)

! EC 20/12/16
! Dissipation according to the Markovian SSE (eq 25 J. Phys: Condens.
! Matter vol. 24 (2012) 273201)
! OR
! Dissipation according to the non-Markovian SSE (eq 24 J. Phys:
! Condens. Matter vol. 24 (2012) 273201) 
! Add a random fluctuation for the stochastic propagation (.not.qjump)
       !if (dis) then
       if (Fdis(1:3).eq."mar".or.Fdis(1:3).eq."nma") then
          call define_h_dis(h_dis,n_ci)
          !if (.not.qjump) then
          if (Fdis(5:9).ne."qjump") then 
             call rnd_noise(w,w_prev,n_ci,first)
             first=.false.
             call add_h_rnd(h_rnd,n_ci,w,w_prev) 
             call add_h_rnd2(h_rnd2,n_ci)
          endif
       endif
!
! INITIAL STEP: dpsi/dt=(psi(2)-psi(1))/dt
       c=c_prev-ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))
!       if (ernd) then
       if (Fdis.eq."ernd") then
          do j=1,n_ci
             c(j) = c(j) - ui*dt*krnd*random_normal()*c_prev(j)
          enddo
       endif
       !if (dis) then
       if (Fdis(1:3).eq."mar".or.Fdis(1:3).eq."nma") then
          c=c-dt*matmul(h_dis,c_prev)
          !if (.not.qjump) then
          !if (tdis.eq.0) then
          if (Fdis(5:9).eq."EuMar") then
          ! Euler-Maruyama
            c=c-ui*sqrt(dt)*matmul(h_rnd,c_prev)-dt*matmul(h_rnd2,c_prev)               
          !elseif (tdis.eq.1) then
          elseif (Fdis(5:9).eq."LeiMa") then
          ! Leimkuhler-Matthews
            c=c-ui*0.5d0*sqrt(dt)*matmul(h_rnd,c_prev)-dt*matmul(h_rnd2,c_prev)
          endif
       endif
       c=c/sqrt(dot_product(c,c))
       c_prev=c

! SP 16/07/17: added call to medium propagation at step 2 to have full output
       i=2
       if (Fmdm(1:3).ne."vac") call prop_medium(i,c_prev,f_prev,h_int)
       call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
! SP 16/07/17: heder called at step 1                                        
       !call out_header
       if (mod(2,n_out).eq.0) call output(2,c,f_prev,h_int)
!
!
! PROPAGATION CYCLE: starts the propagation at timestep 3
! Markovian dissipation (quantum jump) -> dis.and.qjump
!       if (dis.and.qjump) then
       if (Fdis(5:9).eq."qjump") then
          do i=3,n_step
            f_prev2=f(:,i-2)
            f_prev=f(:,i-1)
            h_int=zero 
            if (Fmdm(1:3).ne."vac") call prop_medium(i,c_prev,f_prev,h_int)
            call add_int_vac(f_prev,h_int)
! SC field
            if (Frad.eq."arl".and.i.gt.5) call add_int_rad(mu_prev,mu_prev2,mu_prev3, &
                                                mu_prev4,mu_prev5,h_int)
! Quantum jump (spontaneous or nonradiative relaxation, pure dephasing)
! Algorithm from J. Opt. Soc. Am. B. vol. 10 (1993) 524
            if (i.eq.ijump+1) then
               c=c_prev-ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))-dt*matmul(h_dis,c_prev)
            else 
               c=c_prev2-2.d0*ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))-2.d0*dt*matmul(h_dis,c_prev)
            endif 
! loss_norm computes: 
! norm = 1 - dtot
! dtot = dsp + dnr + dde
! Loss of the norm, dissipative events simulated
! eps -> uniform random number in [0,1]
            call loss_norm(c_prev,n_ci,pjump)
            call random_number(eps)  
            if (dtot.gt.eps)  then
               call quan_jump(c,c_prev,n_ci,pjump)
               ijump=i
               write(*,*) 'Quantum jump at step:', i, (i-1)*dt 
               c_prev=c
            else
                c=c/sqrt(dot_product(c,c))
                c_prev2=c_prev
                c_prev=c
            endif
            call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
            if (mod(i,n_out).eq.0)  call output(i,c,f_prev,h_int)
          enddo
! Markovian dissipation (Euler-Maruyama) -> dis.and.not.qjump
!       elseif (dis.and..not.qjump) then
       elseif (Fdis(1:3).eq."mar") then
          do i=3,n_step
            f_prev2=f(:,i-2)
            f_prev=f(:,i-1)
            h_int=zero 
            if (Fmdm(1:3).ne."vac") call prop_medium(i,c_prev,f_prev,h_int)
            call add_int_vac(f_prev,h_int)
! SC field
            if (Frad.eq."arl".and.i.gt.5) call add_int_rad(mu_prev,mu_prev2,mu_prev3, &
                                                mu_prev4,mu_prev5,h_int)
! Dissipation by a continuous stochastic propagation
            call rnd_noise(w,w_prev,n_ci,first)
            call add_h_rnd(h_rnd,n_ci,w,w_prev)
            if (Fdis(5:9).eq."EuMar") then
            ! Euler-Maruyama 
              c=c_prev-ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))-dt* &
                matmul(h_dis,c_prev)-ui*sqrt(dt)*matmul(h_rnd,c_prev)- &
                dt*matmul(h_rnd2,c_prev)
            elseif (Fdis(5:9).eq."LeiMa") then
              c=c_prev-ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))-dt* &
                matmul(h_dis,c_prev)-ui*0.5d0*sqrt(dt)* &
                matmul(h_rnd,c_prev)-dt*matmul(h_rnd2,c_prev)
            endif
            !c=c/sqrt(dot_product(c,c))
            c_prev=c
            call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
            if (mod(i,n_out).eq.0)  call output(i,c,f_prev,h_int)
          enddo
       elseif (Fdis.eq."nodis".or.Fdis.eq."ernd") then
! No dissipation in the propagation -> .not.dis
          do i=3,n_step
            f_prev2=f(:,i-2)
            f_prev=f(:,i-1)
            h_int=zero 
            if (Fmdm(1:3).ne."vac") call prop_medium(i,c_prev,f_prev,h_int)
            call add_int_vac(f_prev,h_int)
! SC field
            if (Frad.eq."arl".and.i.gt.5) call add_int_rad(mu_prev,mu_prev2,mu_prev3, &
                                                mu_prev4,mu_prev5,h_int)
            c=c_prev2-2.d0*ui*dt*(e_ci*c_prev+matmul(h_int,c_prev))
            if (Fdis.eq."ernd") then
               do j=1,n_ci
                  c(j) = c(j) - 2.d0*ui*dt*krnd*random_normal()*c_prev(j)
               enddo
            endif
            c=c/sqrt(dot_product(c,c))
            c_prev2=c_prev
            c_prev=c
            call do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
            if (mod(i,n_out).eq.0)  call output(i,c,f_prev,h_int)
          enddo
       endif 

! DEALLOCATION AND CLOSING
       deallocate (c,c_prev,c_prev2,h_int)
!       if (dis) then
       if (Fdis(1:3).eq."mar".or.Fdis(1:3).eq."nma") then
          call deallocate_dis()
          deallocate(h_dis)
          !if (.not.qjump) then
          if (Fdis(5:9).ne."qjump") then
             deallocate(h_rnd)
             deallocate(h_rnd2)
             deallocate(w)
             deallocate(w_prev)
          endif 
!          if (qjump) then
          if (Fdis(5:9).eq."qjump") then
             write(*,*)
             write(*,*) 'Total number of quantum jumps',i_sp+i_nr+i_de,&
                        '(',  real(i_sp+i_nr+i_de)/real(n_step),')'
             write(*,*) 'Spontaneous emission quantum jumps', i_sp,    &
                        '(',  real(i_sp)/real(n_step),')'
             write(*,*) 'Nonradiarive relaxation quantum jumps', i_nr, &
                        '(', real(i_nr)/real(n_step),')'
             write(*,*) 'Dephasing relaxation quantum jumps', i_de,    &
                        '(', real(i_de)/real(n_step),')'
             write(*,*)
          endif
       endif
       close (file_c)
       close (file_e)
       close (file_mu)
       close (file_p)
       if(Fmdm(1:3).ne.'vac') call finalize_medium
       return
      end subroutine prop
!
      subroutine create_field
       implicit none
       integer(i4b) :: i,i_max
       real(dbl) :: t_a,ti,tf,arg
       character(15) :: name_f
       allocate (f(3,n_step))
       write(name_f,'(a5,i0,a4)') "field",n_f,".dat"
       open (7,file=name_f,status="unknown")
        f(:,:)=0.d0
        select case (Ffld)
        case ("mdg")
        ! Gaussian modulated sinusoid: exp(-(t-t0)^2/s^2) * sin(wt) 
         select case (npulse)
          case (1)
           do i=1,n_step 
              t_a=dt*(i-1)
              f(:,i)=fmax(:)*exp(-(t_a-t_mid)**2/(sigma**2))*   &
                         sin(omega*t_a)
              if (mod(i,n_out).eq.0) &
              write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
           enddo
          case (2)
           do i=1,n_step
              t_a=dt*(i-1)
              f(:,i)=fmax(:)*exp(-(t_a-t_mid)**2/(sigma**2))*   &
                         sin(omega*t_a)
              f(:,i)=f(:,i) + fmax(:)*exp(-(t_a-(t_mid+tdelay))**2 &
                     /(sigma1**2))*sin(omega1*t_a+pshift) 
              if (mod(i,n_out).eq.0) &
              write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
           enddo
         end select 
        case ("mds")
        ! Cosine^2 modulated sinusoid: 1/2* cos^2(pi(t-t0)/(2t0)) * sin(wt) 
        !          f=0 for t>t0
         i_max=int(t_mid/dt)
         do i=1,2*i_max
          t_a=dt*(dble(i)-1)
          f(:,i)=fmax(:)*cos(pi*(t_a-t_mid)/(2*t_mid))**2/2.d0* &
                 sin(omega*t_a)
         enddo
         do i=2*i_max+1,n_step
          t_a=dt*(i-1)
          f(:,i)=0.
         enddo
        case ("pip")
        ! Pi pulse: cos^2(pi(t-t0)/(2s)) * cos(w(t-t0)) 
         i_max=int(t_mid/dt)
         do i=1,n_step
          t_a=dt*(dble(i)-1)
          f(:,i)=0.
          if (abs(t_a-t_mid).lt.sigma) then
            f(:,i)=fmax(:)*(cos(pi*(t_a-t_mid)/(2*sigma)))**2* &
               cos(omega*(t_a-t_mid))
          endif
         enddo
        case ("sin")
        ! Sinusoid:  sin(wt) 
         do i=1,n_step
          t_a=dt*(dble(i)-1)
          f(:,i)=fmax(:)*sin(omega*t_a)
         enddo
        case ("snd")
        ! Linearly modulated (up to t0) Sinusoid:
        !         0 < t < t0 : t/to* sin(wt) 
        !             t > t0 :       sin(wt) 
         do i=1,n_step
          t_a=dt*(dble(i)-1)
          if (t_a.gt.t_mid) then
            f(:,i)=fmax(:)*sin(omega*t_a)
          else
            if (t_a.gt.zero) f(:,i)=fmax(:)*t_a/t_mid*sin(omega*t_a)
          endif            
         enddo
        case ("gau")
        ! Gaussian pulse: exp(-(t-t0)^2/s^2) 
          select case (npulse)
          case (1)
           do i=1,n_step
              t_a=dt*(i-1)
              f(:,i)=fmax(:)*exp(-pt5*(t_a-t_mid)**2/(sigma**2))
              if (mod(i,n_out).eq.0) &
                write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
           enddo
          case(2)
           do i=1,n_step
              t_a=dt*(i-1)
              f(:,i)=fmax(:)*exp(-pt5*(t_a-t_mid)**2/(sigma**2))
              f(:,i)=f(:,i)+fmax(:)*exp(-pt5*(t_a-(t_mid+tdelay))**2/(sigma**2))
              if (mod(i,n_out).eq.0) &
                write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
           enddo
          end select
! SP 270817: the following (commented) is probably needed for spectra  
         !do i=1,n_step
         ! t_a=dt*(dble(i)-1)
         ! arg=-pt5*(t_a-t_mid)**2/(sigma**2)
         ! if(arg.lt.-50.d0) then 
         !   f(:,i)=zero
         ! else
         !   f(:,i)=fmax(:)*exp(arg)
         ! endif
         !enddo
        case ("css")
        ! Cos^2 pulse (only half a period): cos^2(pi*(t-t0)/(s)) 
         ti=t_mid-sigma/two
         tf=t_mid+sigma/two
         do i=1,n_step
          t_a=dt*(dble(i)-1)
          f(:,i)=zero
          if (t_a.gt.ti.and.t_a.le.tf) then
            f(:,i)=fmax(:)*(cos(pi*(t_a-t_mid)/(sigma)))**2
          endif
         enddo
        case default
         write(*,*)  "Error: wrong field type !"
         stop
        end select
        ! write out field 
        do i=1,n_step
         t_a=dt*(i-1)
         if (mod(i,n_out).eq.0) &
           write (7,'(f12.2,3e22.10e3)') t_a,f(:,i)
        enddo
        close(7)
       return
      end subroutine create_field
!
      subroutine do_mu(c,mu_prev,mu_prev2,mu_prev3,mu_prev4,mu_prev5)
       implicit none
       complex(cmp), intent(IN) :: c(n_ci)
       real(dbl)::mu_prev(3),mu_prev2(3),mu_prev3(3),mu_prev4(3), &
                mu_prev5(3)
       mu_a(1)=dot_product(c,matmul(mut(1,:,:),c))
       mu_a(2)=dot_product(c,matmul(mut(2,:,:),c))
       mu_a(3)=dot_product(c,matmul(mut(3,:,:),c))
! SC save previous mu for radiative damping
       mu_prev5=mu_prev4
       mu_prev4=mu_prev3
       mu_prev3=mu_prev2
       mu_prev2=mu_prev
       mu_prev=mu_a
      return
      end subroutine do_mu

      subroutine output(i,c,f_prev,h_int)     
       implicit none
       integer(i4b), intent(IN) :: i
       complex(cmp), intent(IN) :: c(n_ci)
       real(dbl), intent(IN) :: h_int(n_ci,n_ci)
       real(dbl), intent(IN) :: f_prev(3)
       real(dbl) :: e_a,e_vac,t,g_neq_t,g_neq2_t,g_eq_t,f_med(3)
       real(dbl) :: ctmp(n_ci)
       character(20) :: fmt_ci,fmt_ci1,fmt_ci2
       integer(i4b)    :: itmp,j
        t=(i-1)*dt 
        e_a=dot_product(c,e_ci*c+matmul(h_int,c))
! SC 07/02/16: added printing of g_neq, g_eq 
        if(Fmdm(1:3).ne.'vac') then 
         g_eq_t=e_a
         g_neq_t=e_a
         g_neq2_t=e_a
         e_vac=e_a
         call get_energies(e_vac,g_eq_t,g_neq_t,g_neq2_t)
         write (file_e,'(i8,f14.4,7e20.8)') i,t,e_a,e_vac, &
                   g_eq_t,g_neq2_t,g_neq_t,int_rad,int_rad_int
        else
         write (file_e,'(i8,f14.4,3e22.10)') i,t,e_a,int_rad,int_rad_int
        endif
        itmp = n_ci*(n_ci-1)/2
        write (fmt_ci,'("(i8,f14.4,",I0,"e13.5)")') n_ci
        write (fmt_ci1,'("(i8,f14.4,",I0,"e13.5)")') itmp 
        write (fmt_ci2,'("(i8,f14.4,",I0,"2e13.5)")') 2*itmp
        ctmp(:) = real(c(:)*conjg(c(:)))
        do j=1,n_ci
           if (ctmp(j).lt.10.d-50) ctmp(j)=0.d0
        enddo
        write (file_c,fmt_ci) i,t,ctmp(:)
        write (file_p,fmt_ci) i,t,atan2(aimag(c),real(c))
! SP 10/07/17: commented the following, strange error 
        call wrt_decoherence(i,t,file_dm,file_dp,file_d,fmt_ci1,fmt_ci2,c,n_ci)
        write (file_mu,'(i8,f14.4,3e22.10)') i,t,mu_a(:)
        j=int(dble(i)/dble(n_out))
        if(j.lt.1) j=1
        Sdip(:,1,j)=mu_a(:)
! SP 270817: using get_* functions to communicate with TDPlas
        if(Fmdm(1:3).ne."vac") call get_medium_dip(Sdip(:,2,j))
        Sfld(:,j)=f(:,i)
       return
      end subroutine output
!
      subroutine add_int_vac(f_prev,h_int)
       implicit none
       real(dbl), intent(IN) :: f_prev(3)
       real(dbl), intent(INOUT) :: h_int(n_ci,n_ci)
       ! create the field term of the hamiltonian
! SC 16/02/2016: changed to - sign, 
       h_int(:,:)=h_int(:,:)-mut(1,:,:)*f_prev(1)-             &
                   mut(2,:,:)*f_prev(2)-mut(3,:,:)*f_prev(3)
      return
      end subroutine add_int_vac
!
      subroutine add_int_rad(mu_prev,mu_prev2,mu_prev3,mu_prev4, &
                                                   mu_prev5,h_int)
! SC calculate the Aharonov Lorentz radiative damping
       implicit none
       real(dbl),intent(in) :: mu_prev(3),mu_prev2(3),mu_prev3(3), &
                 mu_prev4(3), mu_prev5(3)
       real(dbl), intent(INOUT) :: h_int(n_ci,n_ci)
       real(dbl) :: d3_mu(3),d2_mu(3),d_mu(3),d2_mod_mu,coeff,scoeff, &
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
       h_int(:,:)=h_int(:,:)-mut(1,:,:)*d3_mu(1)-             &
                   mut(2,:,:)*d3_mu(2)-mut(3,:,:)*d3_mu(3)
       return
      end subroutine add_int_rad
!
      subroutine out_header
! SC write headers to output files, to be completed!
      implicit none
      write(file_e,'(8a)') '#   istep time',' <H(t)>-E_gs(0)', &
              ' DE_vac(t)',' DG_eq(t)',' DG_neq(t)',  '  Const', &
              '  Rad. Int', '  Rad. Ene'   
      write(file_mu,'(5a)') '#   istep time',' dipole-x ', &
              ' dipole-y ',' dipole-z '
      return
      end subroutine out_header
!
      subroutine wrt_decoherence(i,t,int1,int2,int3,char20,char1_20,c,nci)
!------------------------------------------------------------------------
! @brief Print the tridiagional C*_iC_j (i.ne.j) matrix 
! corresponding to the decoherence 
! 
! @date Created   : E. Coccia 3 Feb 2017
! Modified  :
!------------------------------------------------------------------------
  
        implicit none

        integer(i4b),    intent(in)    :: i, nci
        integer(i4b),    intent(in)    :: int1, int2, int3
        character(20), intent(in)    :: char20, char1_20 
        real(dbl),       intent(in)    :: t
        complex(cmp),   intent(in)    :: c(nci)
        complex(cmp)                  :: cc(nci,nci)
        integer(i4b)                   :: j, k
        complex(cmp)                  :: tmp
        real(dbl)                      :: r(nci,nci)
        real(dbl)                      :: p(nci,nci)
        real(dbl)                      :: rtmp, itmp

        tmp=cmplx(0.d0,0.d0)
        r=0.d0
        p=0.d0

        do k=2,nci
           tmp = conjg(c(k))*c(1)
           rtmp=real(tmp)
           itmp=aimag(tmp)
           if (abs(rtmp).lt.10.d-50) tmp=cmplx(0.d0,itmp)
           if (abs(itmp).lt.10.d-50) tmp=cmplx(rtmp,0.d0)
           r(1,k) = abs(tmp)
           p(1,k) = atan2(aimag(tmp),real(tmp))
           cc(1,k) = tmp
        enddo  

        if (nci.gt.2) then
           do j=2,nci
              do k=j+1,nci
                 tmp = conjg(c(k))*c(j)
                 rtmp=real(tmp)
                 itmp=aimag(tmp)
                 if (abs(rtmp).lt.10.d-50) tmp=cmplx(0.d0,itmp)
                 if (abs(itmp).lt.10.d-50) tmp=cmplx(rtmp,0.d0)
                 r(j,k) = abs(tmp)
                 p(j,k) = atan2(aimag(tmp),real(tmp)) 
                 cc(j,k) = tmp
              enddo
           enddo
        endif

        if (nci.gt.2) then
           write(int1,char20,advance='no') i, t, ((r(j,k), k=j+1,nci), j=1,nci-1)
           write(int2,char20,advance='no') i, t, ((p(j,k), k=j+1,nci), j=1,nci-1)    
           write(int3,char1_20) i, t, ((cc(j,k), k=j+1,nci),j=1,nci-1)
        else 
           write(int1,char20,advance='no') i, t, r(1,2)
           write(int2,char20,advance='no') i, t, p(1,2)
           write(int3,char1_20) i, t, cc(1,2)
        endif
         !((FORM(K,L), L=1,10), K=1,10,2)
    
        return

      end subroutine wrt_decoherence

      end module
