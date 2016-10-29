      module readio_medium
      use cav_types       
      use readio          
      use pedra_friends
      implicit none
      save
!
      real(8), parameter :: ev_to_au=0.0367493
      real(8), parameter :: debye_to_au=0.393456
      real(dbl), parameter :: TOANGS=0.52917724924D+00
      real(dbl), parameter :: ANTOAU=1.0D+00/TOANGS
      integer(4) :: nts
      real(8), allocatable :: vts(:,:,:) !transition potentials on tesserae from cis
      real(8), allocatable :: cals(:,:),cald(:,:) !Calderon D and S matrices 
      real(8) :: a_cav,b_cav,c_cav,eps_0,eps_d,mdm_dip(3),tau_deb
! SP 150316: ncycmax max number of scf cycles
      integer(4) :: n_q,nmodes,ncycmax
      integer(4), allocatable :: xmode(:),occmd(:)
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

      
      private
      public read_medium,deallocate_medium,Fint,Feps,Fprop,a_cav,  &
             b_cav,c_cav,eps_0,eps_d,tau_deb,n_q,eps_A,molint,     &
             nmodes,xmode,occmd,nts,vts,eps_gm,eps_w0,f_vel,  &
             iBEM, cals,cald,q0,fr0,Floc,mdm_dip,qtot0, &
             Fdeb,vtsn,mdm_init_prop, &
             ncycmax,thrshld,mix_coef,Fbem,Fcav,Fchr,read_medium_freq     
!
      contains
!
      subroutine read_medium_freq
      Fint='pcm'
      Feps='drl'
      Fprop='ief'
      Fbem='rea'
      Fcav='fil'             
      read(5,*) eps_0,eps_d,eps_A,eps_gm,eps_w0,f_vel
      read(5,*) fmax
      
      return
      end subroutine
!
      subroutine read_medium
       implicit none
       integer(4):: i,lc,db,rf,sc
       character(3) :: which_int, which_eps, which_prop, &
                       which_bound_mat,which_cavity, which_charges, &
                       which_debug,which_localF
       character(6) :: which_mdm_init_prop 
       nts_act=0
! read frequency of updating the interaction potential
       read(5,*) n_q
       read(5,*) 
! read dielectric function type and parameters
       read(5,*) which_eps
       select case (which_eps)
         case ('deb','Deb','DEB')
           Feps='deb'
           read(5,*) eps_0,eps_d,tau_deb
         case ('drl','Drl','DRL')
           Feps='drl'
! SP 30/05/16: changed to have eps_0 and eps_d in input for solids
           read(5,*) eps_0,eps_d,eps_A,eps_gm,eps_w0,f_vel
!SC 16/02/2016: perhaps this correction should not stay in the input routine
           ! correct for finite size effects of nanoparticle with respect to bulk
!           eps_gm=eps_gm+f_vel/sfe_act(1)%r
         case default
           write(*,*) "Error, specify eps(omega) type DEB or DRL"
           stop
       end select
       read(5,*) 
! read interaction type: PCM or Onsager
! which_int refers to the coupling between the molecule and the reaction field
!   ons: the reaction field is considered constant in space and the coupling is -mu*F_RF
!   pcm: the reaction potential is used, the coupling is \int dr rho(r) V_RF(r)
       read(5,*) which_int
       select case (which_int)
        case ('ons','Ons','ONS')
         Fint='ons'
        case ('pcm','Pcm','PCM')
         Fint='pcm'
        case default
         write(*,*) "Error, specify interaction type ONS or PCM"
         stop
       end select
! which_prop refers to which quantity is propagated by equations of motions
!   dip: only the dipolar (i.e., Onsager) reaction field is propagated
!   ief,ied,csm: the apparent charges are propagated, within different schemes
       read(5,*) which_prop
       select case (which_prop)
        case ('ief','Ief','IEF')
         Fprop='ief'
        case ('ied','Ied','IED')
         Fprop='ied'
        case ('csm','Csm','CSM')
       ! COSMO propagation PCM, Onsager Nanoparticle
         Fprop='csm'
        case ('dip','Dip','DIP')
       ! Dipole propagation no charges involved 
         Fprop='dip'
        case default
         write(*,*) "Error, specify the propagation type "
         stop
       end select
! Only some of the combintations mdm/which_int/which_prop make sense/are implemented
! here we sort out what can be done and what cannot, and read in specific input for each option
!
! Here the solvent...
       if (mdm.eq.'sol') then
        if (Fint.eq.'ons') then
         if (Fprop.ne.'dip') then    
          write(*,*) &
      "Error: Ons coupling with solvent PCM charges not implemented yet"
          stop
         else 
          write(6,*) 'This is a standard Onsager run in solution'
          read(5,*) a_cav,b_cav,c_cav
          a_cav=a_cav*antoau
          b_cav=b_cav*antoau
          c_cav=c_cav*antoau
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
! ...and here the NP
       elseif (mdm.eq.'nan') then
         if (Fprop.eq.'dip') then    
          write(*,*)"Error: propagation of dipolar reaction field", & 
            "for the nanoparticle not implemented yet"
          stop
         else 
          read (5,*)
          read (5,*) nmodes
          allocate(xmode(nmodes))
          allocate(occmd(nmodes))
          read (5,*) (xmode(i),i=1,nmodes) 
          read (5,*) (occmd(i),i=1,nmodes) 
          if (nmodes.eq.0) then
            write (6,*) "This is a standard PCM nanoparticle run"
          else
            write (6,*) "This is a quantum PCM nanoparticle run"
          endif
         endif
       endif
! 
! if this is a run with an explicit boundary (fprop.neq.dip), then
! we need to know if the boundary and the matrix are read from outside
! or are produced here and just written out, with no real propagation
       if (Fprop.ne.'dip') then
         read(5,*) 
         read(5,*) which_bound_mat
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
!
! Read in potentials and charges from outside
! matrices are read in BEM_medium
       if(Fbem.eq.'rea') then
         call read_gau_out_medium
!         call read_charges_gau
       elseif (Fbem.eq.'wri') then
! read in the data needed for building the cavity/nanoparticle
         read(5,*) which_cavity
         select case(which_cavity)
         case ('fil','FIL','Fil')
          Fcav='fil'             
         case ('bui','Bui','BUI')
          Fcav='bui'             
          call read_act(5)
         case default
          write(6,*) "Please choose: build or read cavity?"
          stop
         end select
       endif
       read(5,*)
       read(5,*) which_mdm_init_prop
       select case(which_mdm_init_prop)
! SC & SP: determine the starting state & charges of the simulation:
!      SCF_ES: self consistent optimization of the states
!              in the RF of the state defined by ci_ini.inp, charges for such state
!      Default: no self-consistent optimization, initial charges
!               are calculated for the ground state
        case('SCF_ES')
         mdm_init_prop='sce'
         read(5,*) thrshld, mix_coef,ncycmax
         write(6,*) "Do SCF calculation for the state in ci_ini.inp" 
        case default
         mdm_init_prop='nsc'
         write(6,*) "Use GS RF as coded in the ci_energy.inp values" 
       end select
       read(5,*) which_localF
       select case(which_localF)
       case ('loc','Loc','LOC')
        Floc='loc'             
        write(6,*) "Local field effects are included"
       case default
        write(6,*) "Local field effects are NOT included"
       end select
       read(5,*) which_charges
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
       read(5,*) which_debug  
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
       case default
        Fdeb='non'             
       end select
! NEED TO CHECK COMBINATIONS OF OPTIONS THAT ARE NOT IMPLEMENTED YET
       return
      end subroutine
!
!
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
       do j=2,n_ci
         read(7,*) 
         do its=1,nts_act
          read(7,*) vts(1,j,its)
         enddo
         vts(j,1,:)=vts(1,j,:)
       enddo
       !Vij
       do i=2,n_ci
        do j=2,i   
         read(7,*) 
         do its=1,nts_act
          read(7,*) vts(i,j,its)             
         enddo
         vts(j,i,:)=vts(i,j,:)
        enddo
        ! add nuclear potential
        do its=1,nts_act
         vts(i,i,its)=vts(i,i,its)+vtsn(its)
        enddo
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
       if(allocated(cals)) deallocate(cals)
       if(allocated(cald)) deallocate(cald)
       if(allocated(q0)) deallocate(q0)
       if(allocated(vts)) deallocate(vts)
       if(allocated(xmode)) deallocate(xmode)
       if(allocated(occmd)) deallocate(occmd)
       return
      end subroutine
!
!
      end module
