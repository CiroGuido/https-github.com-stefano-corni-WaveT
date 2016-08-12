      Module td_ContMed       
      use readio
      use readio_medium
      use cav_types
      use pedra_friends
      use spectra
      use BEM_medium
      use scf         
      use, intrinsic :: iso_c_binding

      implicit none
      real(dbl), parameter :: TOANGS=0.52917724924D+00
      real(dbl), parameter :: ANTOAU=1.0D+00/TOANGS
      complex(cmp), parameter :: zeroc=(zero,zero)                
!LIGHT SPEED IN AU TO BE REVISED
      real(dbl), parameter :: c=1.37036d2                
      real(dbl) :: qtot,ref    
      real(8), allocatable :: vtsd(:,:,:),vtsv(:,:,:),vtsq(:,:,:)
! SC this is the medium contribution to the hamiltonian
! SP 26/02/16: h_int0 contains the GS PCM contribution
      real(8), allocatable :: h_mdm(:,:),h_mdm0(:,:)
! SP 160316: qst starting charges may be optimised through scf
      real(8), allocatable :: qst(:)
! SC this is the field acting on the molecule due to the medium: fl_t
      real(8) :: fl_t(3),fl_tp(3),fl_tp2(3)
!SP 22/02/16:  mu_mdm is the medium dipole when propagating charges
      real(8) :: dip(3),t,mu_mdm(3)
      real(8) :: fr_t(3), fr_tp(3), fr_tp2(3),dfr_t(3)
      real(8) :: mu_a(3), mu_prev(3),mu_prev2(3), mu0(3)
!SC 06/02/16: eq and neq free energies; h_0 is the energy of the 
      real(8) :: e_vac,g_eq,g_eq_gs,g_neq_0,g_neq1,g_neq1_part,g_neq2, &
                 g_neq2_0
!SP 22/02/16:  k_r constant for dipole propagation                    
      real(8) :: k_r         
!SP 22/02/16:  External charges (qext)/ potential on tesserae (potx)
!SP 22/02/16:  Reaction charges (q)/ potential on tesserae (pot)
!SP 29/05/16:  pot_gs contains vts(1,1,:), but works also with Fint='ons'
      real(dbl), allocatable :: pot_t(:),pot_tp(:),pot_gs(:)
      real(dbl), allocatable :: potx_t(:),potx_tp(:)
      real(dbl), allocatable :: q_t(:),dq_t(:)
      real(dbl), allocatable :: q_tp(:),dq_tp(:)
      real(dbl), allocatable :: qext_t(:),dqext_t(:)
      real(dbl), allocatable :: qext_tp(:),dqext_tp(:)
      integer(4) :: file_med=11 !file numbers
!SP 22/02/16: debug                                                   
      real(8) :: q1,q2,q3,q4 
      save
      private
!SC 07/02/16: added output_gneq
      public init_mdm,prop_mdm,end_mdm,qtot,ref,get_gneq, &
             get_ons 
      contains
!     
! INTERFACE ROUTINES:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SP 23/02/16 : added f_prev2 for onsager local field propagation
      subroutine init_mdm(c_prev,c_prev2,f_prev,f_prev2,h_int)   
      real(dbl), intent(INOUT) :: f_prev(3),f_prev2(3)
      complex(16), intent(INOUT) :: c_prev(n_ci),c_prev2(n_ci)
      real(dbl), intent(INOUT):: h_int(n_ci,n_ci)
      integer(i4b) :: its,i,j                   
      character(20) :: name_f
! OPEN FILES
      write(name_f,'(a9,i0,a4)') "medium_t_",n_f,".dat"
      open (file_med,file=name_f,status="unknown")
      allocate(h_mdm(n_ci,n_ci),h_mdm0(n_ci,n_ci))
      h_mdm=zero
      h_mdm0=zero
      ! Buil matrices/factors for propagation
      call init_BEM 
! SP 18/02/16: changed to try to base everything on Fprop 
! SC: First dipole propagation...
      if(Fprop.eq."dip") then
        call init_dip(c_prev,f_prev)
! SC: ...then charges propagation        
      else 
! SC 03/05/2016: create a new matq0=matqq^-1*matqv that should avoid
!                spurious charge dynamics for stationary states
        !call init_matq0
        call init_potential(c_prev)
        call init_charges(c_prev)
        ! SP initialize some matrices for propagation  coeff
        if (Fint.eq.'pcm') call init_ief      
      endif
      c_prev2=c_prev
      call correct_hamiltonian
      ! SP 25/02/16 Initial gebug routine:
      if(debug) call test_dbg
! SC set the initial values of the solvent component of the 
! neq free energies
      g_neq1=zero
      g_neq2=zero
      return
      end subroutine
!     
      subroutine prop_mdm(i,c_prev,c_prev2,f_prev,f_prev2,h_int)   
      real(dbl), intent(INOUT) :: f_prev(3),f_prev2(3)
      complex(16), intent(IN) :: c_prev(n_ci),c_prev2(n_ci)
      real(dbl), intent(INOUT):: h_int(n_ci,n_ci)
      real(dbl) :: f_d_mat(3,3)
      integer(i4b), intent(IN) :: i                    
      integer(i4b) :: its,k,j                    
      ! Propagate medium only every n_q timesteps  
      !  otherwise, just resum the interaction hamiltonian and exit
       if(mod(i,n_q).ne.0) then
         h_int(:,:)=h_int(:,:)+h_mdm(:,:)
         return
       endif
       t=(i-1)*dt
!SP 18/02/16: changed to try to base everything on Fint and Fprop 
!SP 29/05/16: changed to allow the calculation of g_neq for Fprop=ief and Fint=ons
       if (Fprop.eq.'dip'.or.Fint.eq.'ons') then
         f_d_mat=0.d0
         f_d_mat(1,1)=f_d
         f_d_mat(2,2)=f_d
         f_d_mat(3,3)=f_d
       endif
       if (Fprop.eq.'dip') then
       ! Dipole propagation: 
         call prop_dip(c_prev,f_prev,f_prev2)
         call do_gneq(c_prev,mut,dfr_t,fr_t,f_d_mat,3,-1)
       else
       ! Charges propagation: 
         ! Calculate external potential on tesserae for local field       
         if(localf) call do_ext_potential(f_prev)
         ! Calculate the molecule potential on tesserae
         call do_potential_ts(c_prev,pot_t)
         call prop_chr(c_prev,c_prev2)
         ! Charges normalization 
         ! call norm_ch
         ! Calculate medium's dipole from charges 
         call do_dipole        
         ! Calculate Reaction Field from charges
!SP 29/05/16: changed to allow the calculation of g_neq for Fprop=ief and Fint=ons
         if (Fint.eq.'ons') dfr_t=-fr_t
         call do_field(q_t,fr_t)
         if (Fint.eq.'ons') dfr_t=dfr_t+fr_t
         ! Calculate Local Field from external charges
         if(localf) call do_field(qext_t,fl_t)
         !if (time.eq.endtime) call print_ch
       ! SC calculate free energy:
!SP 29/05/16: changed to allow the calculation of g_neq for Fprop=ief and Fint=ons
         if (Fint.eq.'ons') then 
           call do_gneq(c_prev,mut,dfr_t,fr_t,f_d_mat,3,-1)
         else
           call do_gneq(c_prev,vts,dq_t,q_t,matqd,nts_act,1)
         endif
       endif
       ! Calculate a reference value (testing purposes)
       if(pref.gt.0)  call do_ref
       ! Build the interaction Hamiltonian Reaction/Local
       call do_interaction_h
       if (i.eq.1) then
        write(6,*) "h_mdm at the first propagation step"
        do j=1,n_ci
         do k=1,n_ci
          write (6,*) j,k,h_mdm(j,k)
         enddo
        enddo
       endif
       ! Update the interaction Hamiltonian 
       h_int(:,:)=h_int(:,:)+h_mdm(:,:)
       ! SP 24/02/16  Write output
      if (mod(i,n_out).eq.0) call out_mdm(i)
      return
      end subroutine
!     
      subroutine end_mdm   
       deallocate(h_mdm)
       call end_prop
       call deallocate_BEM    
       close (file_med)
      return
      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!! Initialization/deallocation !!!!!!!!!!!!!!!
!
!
      ! Initialize potential
      subroutine init_potential(c)
       complex(16), intent(IN) :: c(:)
       complex(16) :: c_gs(n_ci)
       integer(i4b) :: its  
       allocate(pot_t(nts_act))
       allocate(pot_gs(nts_act))
! SC 08/04/2016: a routine to test by calculating the potentials from the dipoles
!       call do_vts_from_dip
       call do_potential_ts(c,pot_t)
       if (Fint.eq.'ons') call do_field(q0,fr_tp)
       c_gs(:)=zeroc
       c_gs(1)=onec
       call do_potential_ts(c_gs,pot_gs)
       if(localf) then
         allocate(potx_t(nts_act))
         fl_t(:)=zero
         potx_t(:)=zero     
       endif
       if(Feps.eq."drl") then
         allocate(pot_tp(nts_act))
         pot_tp(:)=pot_t(:)
         if(localf) then
           allocate(potx_tp(nts_act))      
           potx_tp(:)=zero    
         endif
       elseif(Feps.eq."deb") then
         allocate(pot_tp(nts_act))
         pot_tp(:)=pot_t(:)
         if(localf) then
           allocate(potx_tp(nts_act))      
           potx_tp(:)=zero    
         endif
       endif
       return
      end subroutine
!
      ! Initialize charges
      subroutine init_charges(c_prev)
       implicit none
       complex(16), intent(INOUT) :: c_prev(n_ci)
       integer(i4b):: its
       allocate(q_t(nts_act))
       allocate(q_tp(nts_act))
       allocate(dq_t(nts_act))
       if (.not.allocated(q0)) allocate (q0(nts_act))
       ! init the state and the RF before propagation
!SP 29/05/16: pot_gs replaces vts(1,1,:) to allow treating Fprop=ief and Fint=ons
       q0(:)=matmul(matq0,pot_gs)
       g_eq_gs=0.5d0*dot_product(q0,pot_gs)
       write(6,*) 'Medium contribution to ground state free energy:', &
                   g_eq_gs
       ! see readio_medium for mdm_init_prop
       select case (mdm_init_prop)
        case ('nsc')

!         g_neq_0=0.5*f_d*dot_product(mu_prev,mu_prev) &
!                -f_d*dot_product(mu_prev,mu0) &
!                +0.5*f_d*dot_product(mu0,mu0)
!         g_neq_0=-g_neq_0

! SC in principle a non equilibrium self consistency if eps_d=1 is needed
!    here we use the non self-consitent density but use the correct RF
         q_tp=matmul(matqd,pot_gs)
         g_neq_0=0.5d0*dot_product(q_tp,pot_gs)
         g_neq2_0=g_neq_0
         g_neq_0=g_neq_0-dot_product(q_tp,pot_t)
         q_tp=matmul(matqd,pot_t)
         g_neq_0=g_neq_0+0.5d0*dot_product(q_tp,pot_t)
         q_tp=q0+matmul(matqd,(pot_t-pot_gs))
        case ('sce')
         q_tp=mix_coef*matmul(matq0,pot_t)+(1.-mix_coef)*q0
         call do_scf(q_tp,c_prev)
! SP 30/05/16: changed for Fprop=ief and Fint=ons
! update the potential
!         do its=1,nts_act
!          pot_t(its)=dot_product(c_prev,matmul(vts(:,:,its),c_prev))
!         enddo
         call do_potential_ts(c_prev,pot_t)
         g_neq_0=0.5*dot_product(q_tp,pot_t)
         q_tp=matmul(matqd,pot_t)
         g_neq2_0=0.5*dot_product(q_tp,pot_t)
         q_tp=matmul(matq0,pot_t)
       end select
       write(6,*) 'G_neq at t=0:',g_neq_0
       q_t(:)=q_tp(:)
       dq_t(:)=zero  
       if(Fint.eq."ons") call do_field(q_t,fr_0)
       if(localf) then
         allocate(qext_t(nts_act))
         allocate(qext_tp(nts_act))
         allocate(dqext_t(nts_act))
         qext_t(:)=zero 
         qext_tp(:)=zero
         dqext_t(:)=zero 
       endif
       if(Feps.eq."drl") then
         allocate(dq_tp(nts_act))
         dq_tp(:)=zero
         if(localf) then
           allocate(dqext_tp(nts_act))
           dqext_tp(:)=zero
         endif
       endif
       return
      end subroutine
!
      subroutine init_dip(c_prev,f_prev)
      ! inizialize onsager 
       implicit none
       complex(16), intent(INOUT) :: c_prev(:)
       real(dbl), intent(IN) :: f_prev(3)
       mu0(:)=mut(1,1,:)
       if (debug) mu0(:)=zero
       fr_0(:)=f_0*mu0(:)
       g_eq_gs=-0.5d0*dot_product(fr_0,mu0)
       write(6,*) 'Medium contribution to ground state free energy:', &
                   g_eq_gs
       mu_prev(1)=dot_product(c_i,matmul(mut(:,:,1),c_i))
       mu_prev(2)=dot_product(c_i,matmul(mut(:,:,2),c_i))
       mu_prev(3)=dot_product(c_i,matmul(mut(:,:,3),c_i))
       mu_prev2=mu_prev 
       select case (mdm_init_prop)
        case ('nsc')
! SC in principle a non equilibrium self consistency if eps_d=1 is needed
!    here we use the non self-consitent dipole but use the correct RF
         fr_t=fr_0+f_d*(mu_prev-mu0)
         g_neq_0=0.5*f_d*dot_product(mu_prev,mu_prev) &
                -f_d*dot_product(mu_prev,mu0) &
                +0.5*f_d*dot_product(mu0,mu0)
         g_neq_0=-g_neq_0
         g_neq2_0=-0.5*f_d*dot_product(mu0,mu0)
        case ('sce')
         fr_t=mix_coef*f_0*mu_prev+(1.-mix_coef)*fr_0
         call do_scf(fr_t,c_prev)
         mu_prev(1)=dot_product(c_prev,matmul(mut(:,:,1),c_prev))
         mu_prev(2)=dot_product(c_prev,matmul(mut(:,:,2),c_prev))
         mu_prev(3)=dot_product(c_prev,matmul(mut(:,:,3),c_prev))
         mu_prev2=mu_prev
         g_neq_0=-0.5*f_0*dot_product(mu_prev,mu_prev)
         g_neq2_0=-0.5*f_d*dot_product(mu_prev,mu_prev)
       end select
       write(6,*) 'G_neq at t=0:',g_neq_0
       ! assume that fr_tp2=fr_tp is the field from the dipole 
       fr_tp=fr_t
       fr_tp2=fr_t
      ! Local field only for spherical cavity at the moment
      if(localf) then
        fl_t=fx_0*f_prev 
        fl_tp=fl_t
        fl_tp2=fl_t
      endif
      return
      end subroutine
!
      subroutine init_ief
      integer(i4b) :: i,j  
!     Build the contribution in \tilde{Q} Vij and Q_d dVij/dt
      if (Feps.eq.'drl') allocate(vtsq(n_ci,n_ci,nts_act))
      if (Feps.eq.'deb') allocate(vtsv(n_ci,n_ci,nts_act))
      if (Feps.eq.'deb') allocate(vtsd(n_ci,n_ci,nts_act))
      do i=1,n_ci    
        do j=1,n_ci    
          if (Feps.eq.'drl') vtsq(i,j,:)=matmul(matqq(:,:),vts(i,j,:)) 
          if (Feps.eq.'deb') vtsv(i,j,:)=matmul(matqv(:,:),vts(i,j,:))
          if (Feps.eq.'deb') vtsd(i,j,:)=matmul(matqd(:,:),vts(i,j,:))
        enddo
      enddo
      return
      end subroutine
!
      subroutine end_prop
      integer(i4b) :: its  
      if(allocated(vtsq)) deallocate (vtsq)
      if(allocated(vtsd)) deallocate (vtsd)
      if(allocated(vtsv)) deallocate (vtsv)
      if(allocated(pot_t)) deallocate (pot_t)
      if(allocated(potx_t)) deallocate (potx_t)
      if(allocated(pot_tp)) deallocate (pot_tp)
      if(allocated(potx_tp)) deallocate (potx_tp)
      if(allocated(q_t)) deallocate(q_t)
      if(allocated(q_tp)) deallocate(q_tp)
      if(allocated(dq_t)) deallocate(dq_t)
      if(allocated(dq_tp)) deallocate(dq_tp)
      return
      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
!!!!!!!!!!!!!!!!!!!!!!!! routines for field/dipole/potential/
!                            hamiltonian calculation !!!!!!!!!!!!!!!!!! 
!
      subroutine do_ext_potential(fld)
      ! Computes the potential on tesserae given a position-independent 
      !  external field  
       real(dbl), intent(in):: fld(3) 
       integer(i4b) :: its  
        ! Field
        potx_t(:)=zero
        do its=1,nts_act
          potx_t(its)=potx_t(its)-fld(1)*cts_act(its)%x           
          potx_t(its)=potx_t(its)-fld(2)*cts_act(its)%y          
          potx_t(its)=potx_t(its)-fld(3)*cts_act(its)%z         
        enddo
       return
      end subroutine
!
      !
      subroutine do_potential_ts(c,pot)
      ! Builds dipole from CIS coefficients and calculates the potential
      !  on tesserae                 
       complex(16), intent(IN) :: c(n_ci)
       real(dbl), intent(OUT) :: pot(nts_act)
       real(dbl):: diff(3)  
       real(dbl):: dist
       integer(i4b) :: its  
        if(Fint.eq.'ons') then
          dip(1)=dot_product(c,matmul(mut(:,:,1),c))
          dip(2)=dot_product(c,matmul(mut(:,:,2),c))
          dip(3)=dot_product(c,matmul(mut(:,:,3),c))
          pot(:)=zero
          do its=1,nts_act
            diff(1)=-(mol_cc(1)-cts_act(its)%x)
            diff(2)=-(mol_cc(2)-cts_act(its)%y)
            diff(3)=-(mol_cc(3)-cts_act(its)%z)
            dist=sqrt(dot_product(diff,diff))
            pot(its)=pot(its)+dot_product(diff,dip)/(dist**3)
          enddo
        elseif(Fint.eq.'pcm') then
          do its=1,nts_act
            pot(its)=dot_product(c,matmul(vts(:,:,its),c))
          enddo 
        endif
       return
      end subroutine
!
      subroutine do_field(q,f)
      ! Builds field on molecule's center of charge from BEM charges
       real(dbl),intent(out):: f(3)  
       real(dbl),intent(in):: q(nts_act)  
       real(dbl):: diff(3)  
       real(dbl):: dist
       integer(i4b) :: its  
        f(:)=zero
        do its=1,nts_act
          diff(1)=(mol_cc(1)-cts_act(its)%x)
          diff(2)=(mol_cc(2)-cts_act(its)%y)
          diff(3)=(mol_cc(3)-cts_act(its)%z)
          dist=sqrt(dot_product(diff,diff))
          f(:)=f(:)+q(its)*diff(:)/(dist**3)
        enddo
       return
      end subroutine
!
      subroutine do_dipole
      ! Builds dipole from BEM charges
      integer(i4b) :: its  
      mu_mdm(:)=zero
      do its=1,nts_act
        mu_mdm(1)=mu_mdm(1)+q_t(its)*(cts_act(its)%x)
        mu_mdm(2)=mu_mdm(2)+q_t(its)*(cts_act(its)%y)
        mu_mdm(3)=mu_mdm(3)+q_t(its)*(cts_act(its)%z)
      enddo
      return
      end subroutine

      subroutine do_interaction_h
       ! Builds interaction terms 
       integer(i4b):: i,j
! SC 19/03/2016: added a temporarty array to avoid summing the local field
! to ther reaction field for onsager
       real(dbl)::ft_t(3)
! SC 19/03/2016: moved from prop_mdm (so that h_mdm has values assined only here 
       h_mdm(:,:)=zero
       if (Fint.eq.'ons') then
         ft_t(:)=fr_t(:)
!         write(6,*) "ft_t,fr_0",ft_t(3),fr_0(3)
         if(localf) ft_t(:)=ft_t(:)+fl_t(:) 
         h_mdm(:,:)=-h_mdm0(:,:)-mut(:,:,1)*ft_t(1)-mut(:,:,2)*ft_t(2) &
                                                   -mut(:,:,3)*ft_t(3)
       elseif(Fint.eq.'pcm') then
         if(localf) q_t(:)=q_t(:)+qext_t(:)
         do j=1,n_ci   
           do i=1,j       
             h_mdm(i,j)=-h_mdm0(i,j)+dot_product(q_t(:),vts(i,j,:))
             h_mdm(j,i)=h_mdm(i,j)
           enddo
         enddo
       else
         write(*,*) "wrong interaction type "
         stop
       endif
      return
      end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!! Propagation routines !!!!!!!!!!!!!!!!! 
! SP 22/02/16: A general remark: matmul and DGEMV used randomly, 
!             DGEMV shoud be efficient for big matrices              
!
      subroutine prop_chr(c,c_prev)
      complex(16), intent(IN) :: c(n_ci),c_prev(n_ci)
       ! Propagate
       if(Feps.eq."deb") then
         if(Fprop.eq."ief") then
           !call prop_dbg(c,c_prev) 
           call prop_ief_deb
           !call prop_ief_deb_c(c,c_prev) 
         elseif(Fprop.eq."ied") then
           call prop_ied_deb
         else
           write(*,*) "not implemented yet"
           stop
         endif
       elseif(Feps.eq."drl") then
         if(Fprop.eq."ief") then
           call prop_ief_drl 
           !call prop_ief_drl_c(c_prev) 
         elseif(Fprop.eq."csm") then
           call prop_csm_drl  
         else
           write(*,*) "Wrong propagation method"
         endif
         dq_tp=dq_t
         if(localf) dqext_tp=dqext_t
       endif
       if(eq_rf0) q_t=matmul(matq0,pot_t)
       ! SP 22/02/16 dq_t calculated for g_neq
       dq_t=q_t-q_tp
       q_tp=q_t
       pot_tp=pot_t
       if(localf) qext_tp=qext_t
       if(localf) potx_tp=potx_t
      return
      end subroutine
!
!     Onsager dipole propagation
      subroutine prop_dip(c_prev,f_prev,f_prev2)
      ! evolve onsager field/charges 
       implicit none
       complex(16), intent(IN) :: c_prev(n_ci)
       real(dbl), intent(IN):: f_prev(3),f_prev2(3)  
       mu_a(1)=dot_product(c_prev,matmul(mut(:,:,1),c_prev))
       mu_a(2)=dot_product(c_prev,matmul(mut(:,:,2),c_prev))
       mu_a(3)=dot_product(c_prev,matmul(mut(:,:,3),c_prev))
       if (debug) mu_a(:)=dip(:)
       ! propagate onsager factor
! SP 29/02/16 new alghoritm consistent with charge propagation
       !fr_t=f_d*(mu_a-mu_prev2)+(f_0*mu_prev &
       !      -fr_tp)/tau_ons*2.d0*dt*n_q+fr_tp2
       fr_t=f_d*(mu_a-mu_prev)+(f_0*mu_prev &
             -fr_tp)/tau_ons*dt*n_q+fr_tp
       ! update the data      
       if(eq_rf0) fr_t=f_0*mu_a
       dfr_t=fr_t-fr_tp
       mu_prev2=mu_prev
       mu_prev=mu_a
       fr_tp2=fr_tp
       fr_tp=fr_t
       ! Local Field
       if(localf) then
         fl_t=fx_d*(f_prev-f_prev2)+(fx_0*f_prev2 &
               -fl_tp)/taux_ons*dt*n_q+fl_tp
         fl_tp2=fl_tp
         fl_tp=fl_t
         ! Dipole nanoparticle p=(eps-1/(eps+2)/a^3*EF
       endif
      return
      end subroutine
!
! DANGER: initialize factors for solvent/nanoparticle and check signs !!!
      subroutine prop_csm_drl
      ! charge propagation with drude/lorentz and cosmo/onsager equations
      integer(i4b) :: its  
      ! Reaction Field
      q_t=q_tp+dt*dq_tp
      call DGEMV('N',nts_act,nts_act,dt*f_f,sm1,nts_act,pot_t,one_i,&
                     zero,dq_t,one_i)
      dq_t=dq_t+dt*f_w*q_tp+(1-dt*eps_gm)*dq_tp
      ! Local Field
      if(localf) then
        qext_t=qext_tp+dt*dqext_tp
        call DGEMV('N',nts_act,nts_act,dt*f_f,sm1,nts_act,potx_t,one_i,&
                       zero,dqext_t,one_i)
        dqext_t=dqext_t+dt*f_w*qext_tp+(1-dt*eps_gm)*dqext_tp
      endif
      return
      end subroutine
!
      subroutine prop_ief_drl
      ! Charge propagation with drude/lorentz and IEF equations 
      integer(i4b) :: its  
      ! Reaction Field
      q_t=q_tp+dt*dq_tp
      call DGEMV('N',nts_act,nts_act,dt,matqv,nts_act,pot_t,one_i, &
                     zero,dq_t,one_i)
      call DGEMV('N',nts_act,nts_act,-dt,matqq,nts_act,q_tp,one_i,  &
                     one,dq_t,one_i)
      dq_t=dq_t+(1-dt*eps_gm)*dq_tp 
      ! Local Field
      if(localf) then
        qext_t=qext_tp+dt*dqext_tp
        call DGEMV('N',nts_act,nts_act,dt,matqv,nts_act,potx_t,one_i, &
                       zero,dqext_t,one_i)
        call DGEMV('N',nts_act,nts_act,-dt,matqq,nts_act,qext_tp,one_i,&
                       one,dqext_t,one_i)
        dqext_t=dqext_t+(1.d0-dt*eps_gm)*dqext_tp 
      endif
      return
      end subroutine
!
      subroutine prop_ief_drl_c(c)
      ! PCM interaction with drude/lorentz and IEF equations (c) 
      complex(16), intent(IN) :: c(n_ci)
      integer(i4b) :: its  
       ! Reaction Field
       q_t=q_tp+dt*dq_tp
       do its=1,nts_act
         dq_t(its)=dt*dot_product(c,matmul(vtsv(:,:,its),c)) 
       enddo
       call DGEMV('N',nts_act,nts_act,-dt,matqq,nts_act,q_tp,one_i,  &
                      one,dq_t,one_i)
       dq_t=dq_t+(1-dt*eps_gm)*dq_tp 
       ! Local Field
       if(localf) then
         qext_t=qext_tp+dt*dqext_tp
         dqext_t=dt*matmul(matqv,potx_t)
         call DGEMV('N',nts_act,nts_act,-dt,matqq,nts_act,qext_tp,one_i,&
                        one,dqext_t,one_i)
         dqext_t=dqext_t+(1-dt*eps_gm)*dqext_tp 
       endif
      return
      end subroutine
!
      subroutine prop_ief_deb
      ! Charge propagation with debye and IEF equations 
      integer(i4b) :: its  
      ! Reaction Field
      q_t(:)=q_tp(:)-dt*matmul(matqq,q_tp)+dt*matmul(matqv,pot_tp) &
                    +matmul(matqd,pot_t-pot_tp)
      ! Local Field eq.47 JPCA 2015
      if(localf) then
        qext_t(:)=qext_tp(:)-dt*matmul(matqq,qext_tp)+  &
                  dt*matmul(matqv,potx_tp)+matmul(matqd,potx_t-potx_tp)
      endif
      return
      end subroutine
!
      subroutine prop_ied_deb
      ! Charge propagation with dedye and IEF equations one taud 
      integer(i4b) :: its  
      ! Reaction Field
      tau_ons=(eps_d+1.d0)/(eps_0+1.d0)*tau_deb 
      q_t(:)=(1.d0-dt/tau_ons)*q_tp(:)+dt/tau_ons*matmul(matq0,pot_tp) &
                                   +matmul(matqd,pot_t-pot_tp)
      ! Local Field eq.47 JPCA 2015
      if(localf) then
        qext_t(:)=(1.d0-dt/tau_ons)*qext_tp(:)+                          &
          dt/tau_ons*matmul(matq0,potx_tp)+matmul(matqd,potx_t-potx_tp)
      endif
      return
      end subroutine
!
      subroutine prop_ief_deb_c(c,c_prev)
      ! PCM interaction with debye and IEF equations 
      complex(16), intent(IN) :: c(n_ci),c_prev(n_ci)
      integer(i4b) :: its  
      ! Reaction Field
      q_t(:)=q_tp(:)
      call DGEMV('N',nts_act,nts_act,-dt,matqq,nts_act,q_tp,one_i, &
                     one,q_t,one_i)
      ! Add the contribution in \tilde{Q} Vij and Q_d dVij/dt
      do its=1,nts_act
        q_t(its)=q_t(its)                            & 
          +dt*dot_product(c_prev,matmul(vtsv(:,:,its),c_prev))     
      enddo
      do its=1,nts_act
        q_t(its)=q_t(its)                             & 
         +dot_product(c_prev,matmul(vtsd(:,:,its),c)) &
         +dot_product(c,matmul(vtsd(:,:,its),c_prev)) &
         -two*dot_product(c_prev,matmul(vtsd(:,:,its),c_prev))
      enddo
      ! Local Field eq.47 JPCA 2015
      if(localf) then
        qext_t(:)=qext_tp(:)
        call DGEMV('N',nts_act,nts_act,-dt,matqq,nts_act,qext_tp,one_i,&
                       one,qext_t,one_i)
        qext_t(:)=qext_t(:)+dt*matmul(matqv,potx_tp)+ &
                         matmul(matqd,potx_t-potx_tp)  
      endif
      return
      end subroutine
!
      subroutine prop_dbg(c,c_prev)
      ! PCM interaction with debye and IEF equations 
      complex(16), intent(IN) :: c(n_ci),c_prev(n_ci)
      integer(i4b) :: its  
      real(dbl):: dpot, tmp(nts_act)              
      ! Reaction Field
      !q_t(:)=q_tp(:)
      q_t(:)=matmul(matq0,pot_t)
      !
      !tau_ons=(eps_d+one)/(eps_0+one)*tau_deb 
      !q_t(:)=q_t(:)-dt/tau_ons*q_tp(:)
      !q_t(:)=q_t(:)-dt*matmul(matqq,q_tp)
      !call DGEMV('N',nts_act,nts_act,-dt,matqq,nts_act,q_tp,one_i, &
      !               one,q_t,one_i)

      !q_t(:)=q_t(:)+dt/tau_ons*matmul(matq0,pot_tp)
      !q_t(:)=q_t(:)+dt*matmul(matqv,pot_tp)
      !do its=1,nts_act
      !  tmp(its)=dot_product(c,matmul(vts(:,:,its),c))     
      !enddo
      !q_t(:)=q_t(:)+dt/tau_deb*matmul(matq0,tmp) 
      !q_t(:)=q_t(:)+dt*matmul(matqv,tmp) 

      !q_t(:)=q_t(:)+matmul(matqd,pot_t-pot_tp)
      !do its=1,nts_act
      !   tmp(its)=dot_product(c_prev,matmul(vts(:,:,its),c)) &
      !           +dot_product(c,matmul(vts(:,:,its),c_prev)) &
      !   -two*dot_product(c_prev,matmul(vts(:,:,its),c_prev))
      !enddo
      !q_t(:)=q_t(:)+matmul(matqd,tmp)
 
      ! Local Field eq.47 JPCA 2015
      if(localf) then
        qext_t(:)=qext_tp(:)
        call DGEMV('N',nts_act,nts_act,-dt,matqq,nts_act,qext_tp,one_i,&
                       one,qext_t,one_i)
        qext_t(:)=qext_t(:)+dt*matmul(matqv,potx_t)+ &
                         matmul(matqd,potx_t-potx_tp)  
      endif
      return
      end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!! free-energy   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
! SC 06/02/2016 added routine to update non equilibrium
!    free energy
      subroutine do_gneq(c,mu_or_v,df_or_dq,f_or_q,fact_d, &
                           n_coor_or_ts,sig)
       implicit none
       complex(16), intent(in) :: c(n_ci)
       integer(i4b),intent(in) :: n_coor_or_ts,sig
       real(dbl),intent(in) :: mu_or_v(n_ci,n_ci,n_coor_or_ts)
       real(dbl),intent(in) :: df_or_dq(n_coor_or_ts), &
                               f_or_q(n_coor_or_ts),&
                               fact_d(n_coor_or_ts,n_coor_or_ts)
       real(dbl),allocatable :: v_avg(:)
       integer(i4b) :: its
       g_eq=0.d0
       allocate(v_avg(n_coor_or_ts))
       do its=1,n_coor_or_ts
         v_avg(its)=dot_product(c,matmul(mu_or_v(:,:,its),c))
         g_neq1_part=g_neq1_part+sig*v_avg(its)*df_or_dq(its)
         g_eq=g_eq+sig*f_or_q(its)*v_avg(its)
       enddo
       g_eq=0.5*g_eq
       e_vac=2*(g_eq_gs-g_eq)
!       g_neq1=g_neq_0+g_neq1_part-g_eq
       g_neq1=g_neq_0-g_neq1_part
       g_neq2=-sig*(dot_product((mu_or_v(1,1,:)-v_avg), &
               (f_or_q-matmul(fact_d,v_avg)))- &
               0.5*dot_product(v_avg,matmul(fact_d,v_avg)))+ &
               g_neq2_0+e_vac
!       write (6,*) 'g_eq,g_neq1_part',g_eq,g_neq1_part
       g_eq=g_eq_gs-g_eq
! SC to be completed with other means to calculate gneq
! SC: Caricato et al. JCP 2006, currently only for Onsager
       deallocate(v_avg)
       return
      end subroutine
!
! SC 07/02/2016: a small routine to retrive the solvent 
!                contribution to free energies
      subroutine get_gneq(e_vac_t,g_eq_t,g_neq_t,g_neq2_t)
       real(dbl),intent(inout):: e_vac_t,g_eq_t,g_neq_t,g_neq2_t
       e_vac_t=e_vac_t+e_vac
       g_eq_t=g_eq+g_eq_t
       g_neq_t=g_neq1+g_neq_t
       g_neq2_t=g_neq2+g_neq2_t
       return
      end subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!! miscellaneous !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
! SP 24/02/16 WARNING: the normalization must be checked carefully
      subroutine norm_ch       
      integer(i4b) :: its,isph,dlt,md,istart,iend,nsph  
      !dlt=int(nts_act/nesf_act)
      !md=mod(nts_act,nesf_act)
      !nsph=nesf_act 
      !if (md.ne.0) then
      !  nsph=1
      !  dlt=nts_act
      !endif   
      !do isph=1,nsph
      !  istart=1+(isph-1)*dlt
      !  iend=istart+dlt-1
        qtot=zero 
      !  do its=istart,iend
        do its=1,nts_act
          qtot=qtot+q_t(its)
        enddo
      !  do its=istart,iend
        do its=1,nts_act
          !q_t(its)=q_t(its)-qtot/nts_act      
          q_t(its)=q_t(its)-(qtot-qtot0)/nts_act      
        enddo
      !enddo
      return
      end subroutine
!
      subroutine print_ch       
      integer(i4b) :: its  
      open (7,file="fin_ch.dat",status="unknown")
        do its=1,nts_act
          write(*,*) q_t(its)
        enddo  
      return
      end subroutine
!
      subroutine correct_energies
      ! A routine to coorect ES energies with the PCM GS energy 
       real(dbl) :: e0 
       integer(4) :: i
       do i=1,nts_act
         pot_t(i)=dot_product(c_i,matmul(vts(:,:,i),c_i))
       enddo 
       e0=dot_product(matmul(matq0,pot_t),pot_t)
       do i=1,n_ci     
         e_ci(i)=e_ci(i)-e0
       enddo
       return
      end subroutine
!
      subroutine correct_hamiltonian
      ! A routine to coorect the hamiltonian matrx with the PCM initial state.
       implicit none
       integer(4) :: its,i,j 
       if (Fint.eq.'ons') then
         h_mdm0(:,:)=-mut(:,:,1)*fr_tp(1)-mut(:,:,2)*fr_tp(2) &
                                       -mut(:,:,3)*fr_tp(3)
       elseif (Fint.eq.'pcm') then
!SC 04/05/2016: updated to be coerent with both GS and SCF initialization
         do its=1,nts_act     
           h_mdm0(:,:)=h_mdm0(:,:)+q_tp(its)*vts(:,:,its)
         enddo
       endif
       write(6,*) "in correct hamiltonian"
       do i=1,n_ci
        do j=1,n_ci
         write(6,*) i,j,h_mdm0(i,j)
        enddo
       enddo
       return
      end subroutine
!
      subroutine correct_potentials      
      ! A routine to correct the potential matrix with the GS PCM potential.
       integer(4) :: i
       do i=2,n_ci     
         vts(i,i,:)=vts(i,i,:)+pot_gs(:)
       enddo
       return
      end subroutine
!
      subroutine do_ref      
      complex(cmp) :: refc, eps
      integer(i4b) :: its  
      real(dbl):: cc(3),dist,pos(3),dp,emol(3),fac  
       select case (pref)
         case (1)
           pos(1)=sfe_act(1)%x-mol_cc(1)   
           pos(2)=sfe_act(1)%y-mol_cc(2) 
           pos(3)=sfe_act(1)%z-mol_cc(3)
           dist=sqrt(dot_product(pos,pos))
           dp=dot_product(dip,pos)
           emol(:)=(3*dp*pos(:)/dist**2-dip(:))/dist**3
           ref=emol(3)*cts_act(1)%rsfe**3
         case (2)
           call do_eps
           refc=eps_f*ui*exp(-ui*omega*t)
           ref=fmax(3)*real(refc)*cts_act(1)%rsfe**3
         case default
           write(*,*)  "Error: wrong ref type !"
           stop
       end select
      return
      end subroutine
!
! SC a small subroutine to fetch the current value of the onsager field
! from elsewhere
      subroutine get_ons(f_med)
      real(dbl):: f_med(3)
      f_med=fr_t
      return
      end subroutine
!
      subroutine out_mdm(i)
      integer(i4b),intent(IN) :: i  
      real(dbl):: fm(3)
      Sdip(i,:,2)=mu_mdm(:)
      if(Fprop.eq.'dip') then
        write (file_med,'(i8,f12.2,4e22.10)') i,t,fr_t(:),ref
      else
        call do_field(q_t,fm)
        write (file_med,'(i8,f12.2,9e22.10)') i,t,mu_mdm(:),fm(:),  &
                                                       qtot,qtot0
      endif
      return
      end subroutine
!
      subroutine test_dbg 
      integer(i4b):: its,i,j      
      real(dbl) ::q(nts_act),diff(3),dist
       dip(1)=zero
       dip(2)=zero
       dip(3)=one 
       if(Fprop.eq.'dip') then
         write(*,*) "Onsager limits "
         write(*,*) "mu_d: ", f_d
         write(*,*) "mu_0: ", f_0
       else
         do i=1,n_ci
           do j=i,n_ci
             write(*,*) "Potential ", i, j
             dip(:)=mut(i,j,:)
             pot_t(:)=zero
             do its=1,nts_act
               diff(1)=-(mol_cc(1)-cts_act(its)%x)
               diff(2)=-(mol_cc(2)-cts_act(its)%y)
               diff(3)=-(mol_cc(3)-cts_act(its)%z)
               dist=sqrt(dot_product(diff,diff))
               pot_t(its)=dot_product(diff,dip)/(dist**3)
               if (i.eq.j) write(*,*) pot_t(its), vts(i,j,its)+vtsn(its)
               if (i.ne.j) write(*,*) pot_t(its), vts(i,j,its)
             enddo
           enddo
         enddo 
         write(*,*) "q0 read    q0 calculated "
         q=matmul(matq0,vtsn)
         do its=1,nts_act
          write(*,*) q0(its),q_tp(its)+q(its)
         enddo
         q_tp(:)=zero 
         pot_tp(:)=zero
         pot_t(:)=zero
         do its=1,nts_act
           diff(1)=cts_act(its)%x
           diff(2)=cts_act(its)%y
           diff(3)=cts_act(its)%z
           dist=sqrt(dot_product(diff,diff))
           pot_t(its)=pot_t(its)+dot_product(diff,dip)/(dist**3)
         enddo
       endif
      return
      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine init_matq0
      implicit none
      real(dbl),allocatable :: mat_t(:,:),mat_t2(:,:)
      integer(i4b) :: info,its
      integer(i4b),allocatable :: ipiv(:)
      allocate(mat_t(nts_act,nts_act))
      allocate(mat_t2(nts_act,nts_act))
      allocate(ipiv(nts_act))
      mat_t=matqq
      mat_t2=matqv
      call dgetrf(nts_act,nts_act,mat_t,nts_act,ipiv,info)
      call dgetrs('N',nts_act,nts_act,mat_t,nts_act,ipiv,mat_t2,nts_act,info)
      matq0=mat_t2
      deallocate (mat_t)
      deallocate (mat_t2)
      deallocate (ipiv)
      return
      end subroutine
!
      subroutine do_vts_from_dip
       integer(4) :: i,j,its
       real(dbl) :: diff(3),dist,vts_dip
       do its=1,nts_act
        diff(1)=(mol_cc(1)-cts_act(its)%x)
        diff(2)=(mol_cc(2)-cts_act(its)%y)
        diff(3)=(mol_cc(3)-cts_act(its)%z)
        dist=sqrt(dot_product(diff,diff))
        do i=1,n_ci
         do j=i,n_ci
          vts_dip=-dot_product(mut(j,i,:),diff)/dist**3
          if(its.eq.nts_act) write (6,'(2i6,3f8.3,2e13.5)') i,j, &
                          cts_act(its)%x,cts_act(its)%y, &
                          cts_act(its)%z,vts_dip,vts(i,j,its)
          vts(j,i,its)=vts_dip
          vts(i,j,its)=vts_dip
         enddo
        enddo
       enddo
!       do its=1,nts_act
!        write(6,'(i6,3f8.3,1e13.5)') its,&
!          cts_act(its)%x,cts_act(its)%y, &
!          cts_act(its)%z, vts(1,1,its)
!       enddo
       return
       end subroutine
!
      end module
