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
      integer(4) :: n_q,pref,nmodes,ncycmax
      integer(4), allocatable :: xmode(:),occmd(:)
      character(3) :: Fint,Feps,Fprop
! SP 220216: added equilibrium reaction field eq_rf0
! SP 220216: added local field keyword localf
! SP 180516: added pot_file to see if a file with potentials is present
      logical :: molint,iBEM,localf,readmat,pot_file,debug,eq_rf0
! SP 220216: total charge as read from input qtot0
! SP 150316: thrshld threshold for scf convergence
      real(8) :: eps_A,eps_gm,eps_w0,f_vel,qtot0,thrshld,mix_coef
! SP 220216: vtsn: nuclear potential on tesserae
! SC: q0 are the charges in equilibrium with ground state
      real(dbl), allocatable :: q0(:),vtsn(:)
!     fr_0 is the Onsager reaction field in equilibrium with the GS
      real(dbl) :: fr_0(3)
      character(3) :: mdm_init_prop

      
      private
      public read_medium,deallocate_medium,Fint,Feps,Fprop,a_cav,  &
             b_cav,c_cav,eps_0,eps_d,tau_deb,n_q,eps_A,molint,     &
             pref,nmodes,xmode,occmd,nts,vts,eps_gm,eps_w0,iBEM,   &
             cals,cald,q0,fr_0,localf,mdm_dip,readmat,qtot0, &
             debug,vtsn,eq_rf0,mdm_init_prop, &
             ncycmax,thrshld,mix_coef
!
      contains
!
      subroutine read_medium
       implicit none
       integer(4):: i,lc,db,rf,sc
       character(3) :: which_int, which_eps, which_prop
       character(6) :: which_mdm_init_prop
       nts_act=0
       read(5,*) 
! read frequency of updating the interaction potential
       read(5,*) n_q
! read interaction type: PCM or dipolar (onsager)
       read(5,*) which_int
       select case (which_int)
         case ('ons','Ons','ONS')
           Fint='ons'
           if (mdm.eq."sol") then
             read(5,*) a_cav,b_cav,c_cav
             a_cav=a_cav*antoau
             b_cav=b_cav*antoau
             c_cav=c_cav*antoau
           endif
         case ('pcm','Pcm','PCM')
           Fint='pcm'
         case default
           write(*,*) "Error, specify interaction type ONS or PCM"
           stop
       end select
       if (mdm.eq."nan") call read_nanop
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
           eps_gm=eps_gm+f_vel/sfe_act(1)%r
         case default
           write(*,*) "Error, specify eps(omega) type DEB or DRL"
           stop
       end select
       read(5,*) 
! read charges propagation type: ief, ied or dpc or dip 
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
       if(Fprop.ne.'dip'.and.Fint.eq.'pcm') then
         call read_gau_out_medium
         if (mdm.eq."sol") then
           !we need to read CIS0 charges and positions at time 0 
           call read_mat 
           call read_charges_gau
!           call read_interface_gau 
           readmat=.true.
         endif
! SC 08/04/2016: a routine to test by calculating the potentials from the dipoles
         !call do_vts_from_dip
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
       read(5,*) 
       read(5,*) lc
       if (lc.gt.0) localf=.true.   
       read(5,*) db
       if (db.gt.0) debug=.true.   
       read(5,*) rf
       if (rf.gt.0) eq_rf0=.true.
! NEED TO CHECK COMBINATIONS OF OPTIONS THAT ARE NOT IMPLEMENTED YET
       return
      end subroutine
!
      subroutine read_charges_gau
       integer(4) :: i,nts
       open(7,file="charges0.inp",status="old")
         read(7,*) nts
         if(nts_act.eq.0.or.nts.eq.nts_act) then
           nts_act=nts
         else
           write(*,*) "Tesserae number conflict"
           stop
         endif
         if(.not.allocated(q0)) allocate (q0(nts_act))
         if(.not.allocated(cts_act)) allocate (cts_act(nts_act))
         qtot0=zero
         do i=1,nts_act 
           read(7,*) q0(i),cts_act(i)%area
           qtot0=qtot0+q0(i)
         enddo
! SP 16/05/2016: we don't have matq0 here now, correct in the propagation routine 
         !q0=q0-matmul(matq0,vtsn)
       close(7)
       return
      end subroutine
!
      subroutine read_interface_gau
       integer(4) :: i,nts,nsphe
       real(dbl)  :: x,y,z,s,r      
       open(7,file="cavity.inp",status="old")
         !read(7,*)  
         read(7,*) nts,nsphe
         if(nts_act.eq.0.or.nts.eq.nts_act) then
           nts_act=nts
         else
           write(*,*) "Tesserae number conflict"
           stop
         endif
         if(.not.allocated(cts_act)) allocate (cts_act(nts_act))
         do i=1,nsphe
          read(7,*)  x,y,z
         enddo
         do i=1,nts_act 
           read(7,*) x,y,z,s
           cts_act(i)%x=x!*antoau 
           cts_act(i)%y=y!*antoau
           cts_act(i)%z=z!*antoau 
           cts_act(i)%area=s!*antoau*antoau 
           ! SP: this is only for a sphere: test purposes
           cts_act(i)%rsfe=sqrt(x*x+y*y+z*z)
           cts_act(i)%n(1)=cts_act(i)%x/cts_act(i)%rsfe 
           cts_act(i)%n(2)=cts_act(i)%y/cts_act(i)%rsfe
           cts_act(i)%n(3)=cts_act(i)%z/cts_act(i)%rsfe 
         enddo
       close(7)
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
! SP 25/02/16 changed the reading loop N*(N+1)/2
      ! read matrices for propagation
      subroutine read_mat  
       integer :: i,j,nts
       open(7,file="mat_SD.inp",status="old")
       read(7,*) nts
       if(nts_act.eq.0.or.nts.eq.nts_act) then
         nts_act=nts
       else
         write(*,*) "Tesserae number conflict"
         stop
       endif
       allocate (cald(nts_act,nts_act),cals(nts_act,nts_act))
       do j=1,nts_act
        do i=j,nts_act
         read(7,*) cals(i,j), cald(i,j)
         cals(j,i)=cals(i,j)
         cald(j,i)=cald(i,j)
        enddo
       enddo
       close(7)
       return
      end subroutine
!
      subroutine read_nanop
       integer :: i
       open(unit=7,file="nanop.in",status="old",form="formatted")
        read (7,*)
        ! Read sphere info                
        call read_act(7)
        ! Read info                    
        read (7,*)
        read (7,*) molint
        read (7,*) pref    
        read (7,*) nmodes
        allocate(xmode(nmodes))
        allocate(occmd(nmodes))
        read (7,*) (xmode(i),i=1,nmodes) 
        read (7,*) (occmd(i),i=1,nmodes) 
       close(unit=7)
! SC create tessera for the NP
       call pedra_int('met')
       if(Fint.eq.'pcm') then
         pot_file=.false.
         inquire(file='ci_pot.inp',exist=pot_file)       
         if (Fint.eq.'pcm'.and.(.not.pot_file)) then
           call output_surf
           write(6,*) "Created output file with surface points"
           stop
         endif
       endif
! SC 07/02/16: initiate q0 and set to 0 (for solvent calculations is
!   set in the sister routine read_....gau 
       allocate (q0(nts_act))
       q0(:)=zero
       return
      end subroutine
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
