      module readio_medium
      use cav_types       
      use readio          
      use pedra_friends
      implicit none
      save
!
      real(8), parameter    :: ev_to_au=0.0367493
      real(8), parameter    :: debye_to_au=0.393456
      real(dbl), parameter  :: TOANGS=0.52917724924D+00
      real(dbl), parameter  :: ANTOAU=1.0D+00/TOANGS
      integer(4), parameter :: mod_max=10
      integer(4) :: nts
      real(8), allocatable :: vts(:,:,:) !transition potentials on tesserae from cis
      real(8) :: a_cav,b_cav,c_cav,eps_0,eps_d,mdm_dip(3),tau_deb
! SP 150316: ncycmax max number of scf cycles
      integer(4) :: n_q,nmodes,ncycmax,ns
      integer(4), parameter :: nsmax=10
      real(8)      :: xr(nsmax),yr(nsmax),zr(nsmax),rr(nsmax)
! EC 3/5/17
      !integer(4), allocatable :: xmode(:),occmd(:)
      integer   :: xmode(mod_max), occmd(mod_max)
      character(3) :: Fint,Feps,Fprop,Fbem,Fcav,Fchr,Fdeb,Floc
! SP 220216: added equilibrium reaction field eq_rf0
! SP 220216: added local field keyword localf
! SP 180516: added pot_file to see if a file with potentials is present
! SC 11/08/2016: instead of automatic detection of pot_file, now it is part 
!                of the Fbem mechanism (reading or writing)
      logical :: molint,iBEM
! SP 220216: total charge as read from input qtot0
! SP 150316: thrshld threshold for scf convergence
      real(8) :: eps_A,eps_gm,eps_w0,f_vel,qtot0,thrshld,mix_coef
! SP 220216: vtsn: nuclear potential on tesserae
! SC: q0 are the charges in equilibrium with ground state 
! SP: q0 are the initial charges if read_chr is true 
      real(dbl), allocatable :: q0(:),vtsn(:)
!     fr0 is the Onsager reaction field in equilibrium with the GS
      real(dbl) :: fr0(3)
      character(3) :: mdm_init_prop
      character(3) :: which_int, which_eps, which_prop, &
                       which_bound_mat,which_cavity, which_charges, &
                       which_debug,which_localF
      character(6) :: which_mdm_init_prop

      
      private
      public read_medium,deallocate_medium,Fint,Feps,Fprop,a_cav,  &
             b_cav,c_cav,eps_0,eps_d,tau_deb,n_q,eps_A,molint,     &
             nmodes,xmode,occmd,nts,vts,eps_gm,eps_w0,f_vel,  &
             iBEM,q0,fr0,Floc,mdm_dip,qtot0, &
             Fdeb,vtsn,mdm_init_prop, &
             ncycmax,thrshld,mix_coef,Fbem,Fcav,Fchr,read_medium_freq,&
             read_medium_tdplas     
!
      contains

      subroutine read_medium_freq
      mdm='nan'
      Fint='pcm'
      Feps='drl'
      Fprop='ief'
      Fbem='rea'
      Fcav='fil'             

      namelist /freq/ eps_0,eps_d,eps_A,eps_gm,eps_w0,f_vel,fmax
      read(*,nml=freq) 
      return
      end subroutine
!
      subroutine read_medium
       implicit none
       integer(4)   :: i,lc,db,rf,sc!,ns
       !integer(4), parameter :: nsmax=10
       !real(8)      :: xr(nsmax),yr(nsmax),zr(nsmax),rr(nsmax)
       !character(3) :: which_int, which_eps, which_prop, &
       !                which_bound_mat,which_cavity, which_charges, &
       !                which_debug,which_localF
       !character(6) :: which_mdm_init_prop 
       nts_act=0

       if (mdm.eq.'nan') then
          namelist /nanoparticle/ n_q,which_eps,eps_0,eps_d,eps_A, &
                                  eps_gm,eps_w0,f_vel,which_int,   &
                                  which_prop,nmodes,xmode,occmd,   &
                                  which_bound_mat,which_cavity,    &
                                  which_mdm_init_prop,thrshld,     &
                                  mix_coef,ncycmax,which_localF,   &
                                  which_charges,which_debug,tau_deb 
       elseif (mdm.eq.'sol') then
          namelist /solvent/ n_q,which_eps,eps_0,eps_d,eps_A, &
                             eps_gm,eps_w0,f_vel,which_int,   &
                             which_prop,a_cav,b_cav,c_cav,    &
                             which_bound_mat,which_cavity,    &
                             which_mdm_init_prop,thrshld,     &
                             mix_coef,ncycmax,which_localF,   &
                             which_charges,which_debug,tau_deb 
       endif

       if (mdm.eq.'nan') then
          call init_nml_nanoparticle() 
          read(*,nml=nanoparticle) 
          call write_nml_nanoparticle()
       elseif (mdm.eq.'sol') then
          call init_nml_solvent()
          read(*,nml=solvent)
          call write_nml_solvent()
       endif

       return

      end subroutine

!     Used in main_tdplas.f90 
      subroutine read_medium_tdplas
       implicit none
       integer(4)   :: i,lc,db,rf,sc,ns
       integer(4),parameter :: nsmax=10
       character(3) :: which_eps,which_cavity
       real(8)      :: xr(nsmax),yr(nsmax),zr(nsmax),rr(nsmax)
       nts_act=0
       mdm='nan'

       namelist /tdplas/ which_eps,eps_0,eps_d,tau_deb,eps_A,eps_gm, &
                         eps_w0,f_vel,which_cavity,xr,yr,zr,rr,ns
       read(*,nml=tdplas)
! read dielectric function type and parameters
       select case (which_eps)
         case ('deb','Deb','DEB')
           Feps='deb'
           write(6,*) "Debye dielectric constant"
         case ('drl','Drl','DRL')
           Feps='drl'
! SP 30/05/16: changed to have eps_0 and eps_d in input for solids
           write(6,*) "Drude-Lorentz dielectric constant"
         case default
           write(*,*) "Error, specify eps(omega) type DEB or DRL"
           stop
       end select
! Interaction is PCM
       Fint='pcm'
! This run prepare matrices
       Fbem='wri'
       select case(which_cavity)
       case ('fil','FIL','Fil')
        Fcav='fil'             
       case ('gms','GMS','Gms')
        Fcav='gms'             
       case ('bui','Bui','BUI')
        Fcav='bui'             
        if (ns.gt.nsmax) then
           write(*,*) 'ERROR: ns is larger than', nsmax
           stop
        endif
        call read_act(xr,yr,zr,rr,ns,nsmax)
       case default
        write(6,*) "Please choose: build or read cavity?"
        stop
       end select
! This is to follow the correct route in init_BEM
       Fprop='ief'
       return
      end subroutine
!

      ! read transition potentials on tesserae
      subroutine read_gau_out_medium
       integer(4) :: i,j,its,nts
       real(dbl)  :: scr       
       open(7,file="ci_pot.inp",status="old")
       read(7,*) nts
       if(nts_act.eq.0.or.nts.eq.nts_act) then
         nts_act=nts
       else
         write(*,*) "Tesserae number conflict"
         stop
       endif
       allocate (vts(n_ci,n_ci,nts_act))
       allocate (vtsn(nts_act))
       ! V00
       read(7,*) 
       do its=1,nts_act
        read(7,*) vts(1,1,its),scr,vtsn(its)
! SC 08/04/2016: moved following lines inside the do to avoid error
       ! add nuclear potential
        vts(1,1,its)=vts(1,1,its)+vtsn(its)
       enddo
       !V0j
       do j=2,n_ci_read
         read(7,*) 
         if (j.le.n_ci) then
          do its=1,nts_act
          read(7,*) vts(1,j,its)
          enddo
          vts(j,1,:)=vts(1,j,:)
         else
          do its=1,nts_act
           read(7,*)
          enddo
         endif
       enddo
       !Vij
       do i=2,n_ci_read
        do j=2,i   
         read(7,*) 
         if (i.le.n_ci.and.j.le.n_ci) then
          do its=1,nts_act
           read(7,*) vts(i,j,its)             
          enddo
          vts(j,i,:)=vts(i,j,:)
         else
          do its=1,nts_act
           read(7,*) 
          enddo
         endif
        enddo
        ! add nuclear potential
        if (i.le.n_ci) then
         do its=1,nts_act
          vts(i,i,its)=vts(i,i,its)+vtsn(its)
         enddo
        endif
       enddo
       close(7)
       return
      end subroutine
!
!
      subroutine output_surf
       integer :: i
       open(unit=7,file="surface.xyz",status="unknown",form="formatted")
        write (7,*) nts_act
        do i=1,nts_act
          write (7,'(3F22.10)') cts_act(i)%x,cts_act(i)%y,cts_act(i)%z
        enddo
       close(unit=7)
      end subroutine
!
      subroutine deallocate_medium
       if(allocated(q0)) deallocate(q0)
       if(allocated(vts)) deallocate(vts)
       !if(allocated(xmode)) deallocate(xmode)
       !if(allocated(occmd)) deallocate(occmd)
       return
      end subroutine

      subroutine init_nml_nanoparticle()
!------------------------------------------------------------------------
! @brief Initialize variables in the namelist nanoparticle 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param n_q,which_eps,eps_0,eps_d,eps_A,eps_gm,eps_w0,f_vel,which_int,
!        which_prop,nmodes,xmode,occmd,which_bound_mat,which_cavity,
!        which_mdm_init_prop,thrshld,mix_coef,ncycmax,which_localF,
!        which_charges,which_debug,tau_deb         
!------------------------------------------------------------------------

       n_q=1
       which_eps='deb'
       eps_0=35.688
       eps_d=1.806874
       eps_A=0.110224
       eps_gm=0.000757576
       eps_w0=0.0
       f_vel=0.0
       which_int='pcm'
       which_prop='ief'
       nmodes=0
       xmode=0.0
       occmd=0.0
       which_bound_mat='rea'
       which_cavity='bui'
       which_mdm_init_prop='non'
       thrshld=10
       mix_coef=0.2
       ncycmax=600
       which_localF='non'
       which_charges='fro'
       which_debug='non'
       tau_deb=1000.

       return

      end subroutine init_nml_nanoparticle

      subroutine init_nml_solvent()
!------------------------------------------------------------------------
! @brief Initialize variables in the namelist solvent 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param n_q,which_eps,eps_0,eps_d,eps_A,eps_gm,eps_w0,f_vel,which_int, 
!        which_prop,a_cav,b_cav,c_cav,which_bound_mat,which_cavity,
!        which_mdm_init_prop,thrshld,mix_coef,ncycmax,which_localF,
!        which_charges,which_debug,tau_deb
!------------------------------------------------------------------------

       n_q=1
       which_eps='deb'
       eps_0=35.688
       eps_d=1.806874
       eps_A=0.110224
       eps_gm=0.000757576
       eps_w0=0.0
       f_vel=0.0
       which_int='pcm' 
       which_prop='ief'
       a_cav=10.0
       b_cav=10.0
       c_cav=10.0
       which_bound_mat='rea'
       which_cavity='bui'
       which_mdm_init_prop='non'
       thrshld=10
       mix_coef=0.2
       ncycmax=600
       which_localF='non'
       which_charges='fro'
       which_debug='non'
       tau_deb=1000.

       return

      end subroutine init_nml_solvent

      subroutine write_nml_nanoparticle()
!------------------------------------------------------------------------
! @brief Write variables in the namelist nanoparticle and put conditions 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param n_q,which_eps,eps_0,eps_d,eps_A,eps_gm,eps_w0,f_vel,which_int,
!        which_prop,nmodes,xmode,occmd,which_bound_mat,which_cavity,
!        which_mdm_init_prop,thrshld,mix_coef,ncycmax,which_localF,
!        which_charges,which_debug,tau_deb 
!------------------------------------------------------------------------

       write(*,*) 'Frequency of updating the interaction potential', n_q 
       select case (which_eps)
         case ('deb','Deb','DEB')
           Feps='deb'
         case ('drl','Drl','DRL')
           Feps='drl'
! SP 30/05/16: changed to have eps_0 and eps_d in input for solids
!SC 16/02/2016: perhaps this correction should not stay in the input
!routine
           ! correct for finite size effects of nanoparticle with
           ! respect to bulk
!           eps_gm=eps_gm+f_vel/sfe_act(1)%r
         case default
           write(*,*) "Error, specify eps(omega) type DEB or DRL"
           stop
       end select        
! read interaction type: PCM or Onsager
! which_int refers to the coupling between the molecule and the reaction field
!   ons: the reaction field is considered constant in space and the coupling is
!   -mu*F_RF
!   pcm: the reaction potential is used, the coupling is \int dr rho(r) V_RF(r)
       select case (which_int)
        case ('ons','Ons','ONS')
         Fint='ons'
         write(*,*) 'The coupling is of Onsager type'
        case ('pcm','Pcm','PCM')
         Fint='pcm'
         write(*,*) 'The coupling is of PCM type'
        case default
         write(*,*) "Error, specify interaction type ONS or PCM"
         stop
       end select
! which_prop refers to which quantity is propagated by equations of motions
!   dip: only the dipolar (i.e., Onsager) reaction field is propagated
!   ief,ied,csm: the apparent charges are propagated, within different schemes
       select case (which_prop)
        case ('ief','Ief','IEF')
         write(*,*) 'Apparent charges are propagated'
         Fprop='ief'
        case ('ied','Ied','IED')
         write(*,*) 'Apparent charges are propagated'
         Fprop='ied'
        case ('csm','Csm','CSM')
       ! COSMO propagation PCM, Onsager Nanoparticle
         write(*,*) 'Apparent charges are propagated'
         Fprop='csm'
        case ('dip','Dip','DIP')
       ! Dipole propagation no charges involved 
         write(*,*) 'Only the dipolar reaction field is propagated'
         Fprop='dip'
        case default
         write(*,*) "Error, specify the propagation type "
         stop
       end select
       if (Fprop.eq.'dip') then
          write(*,*)"Error: propagation of dipolar reaction field", &
            "for the nanoparticle not implemented yet"
          stop
       else
          if (nmodes.eq.0) then
            write (6,*) "This is a standard PCM nanoparticle run"
          else
            write (6,*) "This is a quantum PCM nanoparticle run"
            write(*,*) 'Number of modes', nmodes
            write(*,*) 'Modes', xmode(1:nmodes)
            write(*,*) 'Occupation numbers', occmd(1:nmodes)
          endif
       endif
! if this is a run with an explicit boundary (fprop.neq.dip),
! then
! we need to know if the boundary and the matrix are read from outside
! or are produced here and just written out, with no real propagation
       if (Fprop.ne.'dip') then
         select case(which_bound_mat)
         case ('rea','Rea','REA')
          Fbem='rea'
          write(6,*) "This is full run reading matrix and boundary"
         case ('wri','Wri', 'WRI')
          Fbem='wri'
          write(6,*) "This run just writes matrices and boundary"
         case default
          write(*,*) "Error, specify if boundary data & matrices", &
             "are read (read) or made and written out (write) "
          stop
         end select
       endif
! Read in potentials and charges from outside
! matrices are read in BEM_medium
       if(Fbem.eq.'rea') then
         call read_gau_out_medium
!         call read_charges_gau
       elseif (Fbem.eq.'wri') then
! read in the data needed for building the cavity/nanoparticle
         select case(which_cavity)
         case ('fil','FIL','Fil')
          Fcav='fil'
         case ('gms','GMS','Gms')
          Fcav='gms'
         case ('bui','Bui','BUI')
          Fcav='bui'
           if (ns.gt.nsmax) then
               write(*,*) 'ERROR: ns is larger than', nsmax
               stop
          endif
          call read_act(xr,yr,zr,rr,ns,nsmax)
         case default
          write(6,*) "Please choose: build or read cavity?"
          stop
         end select
       endif
       select case(which_mdm_init_prop)
! SC & SP: determine the starting state & charges of the simulation:
!      SCF_ES: self consistent optimization of the states
!              in the RF of the state defined by ci_ini.inp, charges for
!              such state
!      Default: no self-consistent optimization, initial charges
!               are calculated for the ground state
        case('SCF_ES')
         mdm_init_prop='sce'
         write(6,*) "Do SCF calculation for the state in ci_ini.inp"
         write(*,*) 'SCF threshold', thrshld
         write(*,*) 'Mixing coefficient', mix_coef
         write(*,*) 'Max number of SCF cycles', ncycmax
        case default
         mdm_init_prop='nsc'
         write(6,*) "Use GS RF as coded in the ci_energy.inp values"
       end select
       select case(which_localF)
       case ('loc','Loc','LOC')
        Floc='loc'
        write(6,*) "Local field effects are included"
       case default
        write(6,*) "Local field effects are NOT included"
       end select
       select case(which_charges)
       case ('vac','VAC','Vac')
        Fchr='vac'
       case ('fro','Fro','FRO')
        Fchr='fro'
       case ('rea','Rea','REA')
        Fchr='rea'
       case default
        write(6,*) "Please choose how to initialize charges"
        stop
       end select
       select case(which_debug)
       case ('deb','Deb','DEB')
        Fdeb='deb'
       case ('equ','Equ','EQU')
        Fdeb='equ'
        write(6,*) "DEBUG: Equilibrium reaction field calculation"
       case ('n-r','N-r','n-R','N-R')
        Fdeb='n-r'
        Floc='non'
        write(6,*) "DEBUG: Nanoparticle reaction field"
       case ('n-l','N-l','n-L','N-L')
        Fdeb='n-l'
        Ffld='snd'
        Floc='loc'
        write(6,*) "DEBUG: Nanoparticle local field"
       case ('vmu','Vmu','VMU','VMu')
        Fdeb='vmu'
        write(6,*) "DEBUG: Potentials calculated from Dipoles "
       case ('off','Off','OFF')
        Fdeb='off'
        write(6,*) "DEBUG: Molecule - Medium interaction turned off"
       case default
        Fdeb='non'
       end select

       return

      end subroutine write_nml_nanoparticle

      subroutine write_nml_solvent()
!------------------------------------------------------------------------
! @brief Write variables in the namelist solvent and put conditions 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param n_q,which_eps,eps_0,eps_d,eps_A,eps_gm,eps_w0,f_vel,which_int, 
!        which_prop,a_cav,b_cav,c_cav,which_bound_mat,which_cavity,
!        which_mdm_init_prop,thrshld,mix_coef,ncycmax,which_localF,
!        which_charges,which_debug,tau_deb
!------------------------------------------------------------------------

       write(*,*) 'Frequency of updating the interaction potential', n_q 
       select case (which_eps)
         case ('deb','Deb','DEB')
           Feps='deb'
         case ('drl','Drl','DRL')
           Feps='drl'
! SP 30/05/16: changed to have eps_0 and eps_d in input for solids
!SC 16/02/2016: perhaps this correction should not stay in the input
!routine
           ! correct for finite size effects of nanoparticle with
           ! respect to bulk
!           eps_gm=eps_gm+f_vel/sfe_act(1)%r
         case default
           write(*,*) "Error, specify eps(omega) type DEB or DRL"
           stop
       end select
! read interaction type: PCM or Onsager
! which_int refers to the coupling between the molecule and the reaction
! field
!   ons: the reaction field is considered constant in space and the
!   coupling is
!   -mu*F_RF
!   pcm: the reaction potential is used, the coupling is \int dr rho(r)
!   V_RF(r)
       select case (which_int)
        case ('ons','Ons','ONS')
         Fint='ons'
         write(*,*) 'The coupling is of Onsager type'
        case ('pcm','Pcm','PCM')
         Fint='pcm'
         write(*,*) 'The coupling is of PCM type'
        case default
         write(*,*) "Error, specify interaction type ONS or PCM"
         stop
       end select
! which_prop refers to which quantity is propagated by equations of
! motions
!   dip: only the dipolar (i.e., Onsager) reaction field is propagated
!   ief,ied,csm: the apparent charges are propagated, within different
!   schemes
       select case (which_prop)
        case ('ief','Ief','IEF')
         write(*,*) 'Apparent charges are propagated'
         Fprop='ief'
        case ('ied','Ied','IED')
         write(*,*) 'Apparent charges are propagated'
         Fprop='ied'
        case ('csm','Csm','CSM')
       ! COSMO propagation PCM, Onsager Nanoparticle
         write(*,*) 'Apparent charges are propagated'
         Fprop='csm'
        case ('dip','Dip','DIP')
       ! Dipole propagation no charges involved 
         write(*,*) 'Only the dipolar reaction field is propagated'
         Fprop='dip'
        case default
         write(*,*) "Error, specify the propagation type "
         stop
       end select
       if (Fint.eq.'ons') then
        if (Fprop.ne.'dip') then
          write(*,*) &
      "Error: Ons coupling with solvent PCM charges not implemented yet"
          stop
        else
          write(6,*) 'This is a standard Onsager run in solution'
          a_cav=a_cav*antoau
          b_cav=b_cav*antoau
          c_cav=c_cav*antoau
          write(*,*) 'Size of cavity:'
          write(*,*) 'a (a.u)', a_cav
          write(*,*) 'b (a.u)', b_cav
          write(*,*) 'c (a.u)', c_cav
        endif
       elseif (Fint.eq.'pcm') then
        if (Fprop.eq.'dip') then
          write(*,*) &
      "Error: PCM incompatible with dipolar reaction field propagation"
          stop
        else
          write(6,*) 'This is a standard PCM run in solution'
        endif
       endif
! if this is a run with an explicit boundary (fprop.neq.dip), then
! we need to know if the boundary and the matrix are read from outside
! or are produced here and just written out, with no real propagation
       if (Fprop.ne.'dip') then
         select case(which_bound_mat)
         case ('rea','Rea','REA')
          Fbem='rea'
          write(6,*) "This is full run reading matrix and boundary"
         case ('wri','Wri', 'WRI')
          Fbem='wri'
          write(6,*) "This run just writes matrices and boundary"
         case default
          write(*,*) "Error, specify if boundary data & matrices", &
             "are read (read) or made and written out (write) "
          stop
         end select
       endif
! Read in potentials and charges from outside
! matrices are read in BEM_medium
       if(Fbem.eq.'rea') then
         call read_gau_out_medium
!         call read_charges_gau
       elseif (Fbem.eq.'wri') then
! read in the data needed for building the cavity/nanoparticle
         select case(which_cavity)
         case ('fil','FIL','Fil')
          Fcav='fil'
         case ('gms','GMS','Gms')
          Fcav='gms'
         case ('bui','Bui','BUI')
          Fcav='bui'
           if (ns.gt.nsmax) then
               write(*,*) 'ERROR: ns is larger than', nsmax
               stop
          endif
          call read_act(xr,yr,zr,rr,ns,nsmax)
         case default
          write(6,*) "Please choose: build or read cavity?"
          stop
         end select
       endif
       select case(which_mdm_init_prop)
! SC & SP: determine the starting state & charges of the simulation:
!      SCF_ES: self consistent optimization of the states
!              in the RF of the state defined by ci_ini.inp, charges for
!              such state
!      Default: no self-consistent optimization, initial charges
!               are calculated for the ground state
        case('SCF_ES')
         mdm_init_prop='sce'
         write(6,*) "Do SCF calculation for the state in ci_ini.inp"
         write(*,*) 'SCF threshold', thrshld
         write(*,*) 'Mixing coefficient', mix_coef
         write(*,*) 'Max number of SCF cycles', ncycmax
        case default
         mdm_init_prop='nsc'
         write(6,*) "Use GS RF as coded in the ci_energy.inp values"
       end select
       select case(which_localF)
       case ('loc','Loc','LOC')
        Floc='loc'
        write(6,*) "Local field effects are included"
       case default
        write(6,*) "Local field effects are NOT included"
       end select
       select case(which_charges)
       case ('vac','VAC','Vac')
        Fchr='vac'
       case ('fro','Fro','FRO')
        Fchr='fro'
       case ('rea','Rea','REA')
        Fchr='rea'
       case default
        write(6,*) "Please choose how to initialize charges"
        stop
       end select
       select case(which_debug)
       case ('deb','Deb','DEB')
        Fdeb='deb'
       case ('equ','Equ','EQU')
        Fdeb='equ'
        write(6,*) "DEBUG: Equilibrium reaction field calculation"
       case ('n-r','N-r','n-R','N-R')
        Fdeb='n-r'
        Floc='non'
        write(6,*) "DEBUG: Nanoparticle reaction field"
       case ('n-l','N-l','n-L','N-L')
        Fdeb='n-l'
        Ffld='snd'
        Floc='loc'
        write(6,*) "DEBUG: Nanoparticle local field"
       case ('vmu','Vmu','VMU','VMu')
        Fdeb='vmu'
        write(6,*) "DEBUG: Potentials calculated from Dipoles "
       case ('off','Off','OFF')
        Fdeb='off'
        write(6,*) "DEBUG: Molecule - Medium interaction turned off"
       case default
        Fdeb='non'
       end select

       return

      end subroutine write_nml_solvent



 end module
