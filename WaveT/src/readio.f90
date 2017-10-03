      module readio
      use constants       
      implicit none
      save
!
      integer(i4b) :: n_f,n_ci,n_ci_read,n_step,n_out 
      !integer(i4b) :: imar !imar=0 Markvian, imar=1, nonMarkovian
      integer(i4b) :: i_sp=0,i_nr=0,i_de=0 !counters for quantum jump occurrences
      integer(i4b) :: nrnd !the time step for Euler-Maruyama is dt/nrnd
! SP 17/07/17: Changed to char flags
      integer(i4b) :: nr_typ !input integer for type of decay for the internal conversion
      integer(i4b) :: idep   !input integer for the dephasing operator
      integer(i4b) :: tdis   !input integer for Euler tdis=0, Matthews tdis=1 
      character(flg) :: Fdis_rel  !< Flag for decay for internal conversion, relaxation via dipole "dip" or matrix "mat"
      character(flg) :: Fdis_deph !< Flag for dephasing operator: exp(i delta_i)|i><i| "exp" or |i><i|-|0><0| "i-0"
! SP270917: added for merging to newer master
      integer(i4b) :: npulse !number of pulses
      real(dbl), allocatable :: mut_np2(:,:) !squared dipole from NP
      real(dbl) :: tdelay, pshift  ! time delay and phase shift with two pulses
      real(dbl) :: omega1, sigma1  ! for the second pulse
!
      real(dbl), allocatable :: c_i(:),e_ci(:)  ! coeff and energy from cis
      real(dbl), allocatable :: mut(:,:,:) !transition dipoles from cis
      real(dbl), allocatable :: nr_gam(:), de_gam(:) !decay rates for nonradiative and dephasing events
      real(dbl), allocatable :: sp_gam(:) !decay rate for spontaneous emission  
      real(dbl), allocatable :: sp_fact(:) !multiplicative factor for the decay rate for spontaneous emission
      real(dbl), allocatable :: tmom2(:) !square transition moments i->0
      real(dbl), allocatable :: delta(:) !phases randomly added during the propagation 
      real(dbl), allocatable :: de_gam1(:) !combined decay rates for idep=1
      real(dbl), allocatable :: tomega(:)  !generalized transition frequencies 
      real(dbl) :: dt,tau(2),start, krnd
! SP 17/07/17: Changed to char flags
      !logical :: dis !turns on the dissipation
      !logical :: qjump ! =.true. quantum jump, =.false. stochastic propagation
      !logical :: ernd=.false. !add normal number to E: E -> E + krnd*rnd()
! Global flags         
      character(flg) :: Fdis !< Flag for dissipation type: 
                             !! Markovian     : quantum jumps "mar-qjump", Euler-Maruyama "mar-EuMar", Leimkuhler-Matthews "mar-LeiMa"
                             !! Non-Markovian : quantum jumps "nma-qjump", Continuous stochastic propagator "nma-cstoc"
                             !! Random        : random energy term "ernd" 
      ! qjump works for the Markovian case
      ! Stochastic propagation from:
      ! Appl. Math. Res. Express vol. 2013 34-56 (2013)
      ! IMA J. Numer. Anal. vol. 36 13-79 (2016)
      real(dbl) :: t_mid,sigma,dir_ft(3),fmax(3),omega,mol_cc(3)
      character(flg) :: Ffld !< Field type 
      character(flg) :: Fmdm !< Flag for medium type, this will be defined in readio_medium after separation
      character(flg) :: Frad !< Flag for radiative damping 
      character(flg) :: Fres !< Flag for restart
      character(flg) :: Fful !< Flag for relaxation matrix
!      character(flg) :: Fpulse !< Flag for pulse             
! Flags read from input file
      character(flg) :: medium,radiative,dissipative,pulse
      character(flg) :: dis_prop
      character(flg) :: restart 
      character(flg) :: full ! full or only |e> -> |0> relaxation
      integer(i4b) :: iseed  ! seed for random number generator
      integer(i4b) :: nexc   ! number of excited states
      integer(i4b) :: nrel   ! number of relaxation channels
      integer(i4b) :: nf     ! number of relaxation channels (including |e> -> |0> terms) 
      integer(i4b) :: i,nspectra
! kind of surrounding medium and shape of the impulse
!     Fmdm=sol: solvent
!     Fmdm=nan: nanoparticle
!     Fmdm=vac: no medium
!
!     Ffld=gau: gaussian impulse
!SC
!     Frad: wheter or not to apply radiative damping
!     TO BE COMPLETED, SEE PROPAGATE.F90      
      real(dbl) :: eps_A,eps_gm,eps_w0,f_vel

      
      private
      public read_input,n_ci,n_ci_read,n_step,dt,           &
             Ffld,t_mid,sigma,omega,fmax,restart,           & 
             Fmdm,mol_cc,tau,start,c_i,e_ci,mut,            &
             Frad,n_out,iseed,n_f,dir_ft,full,              &
! SP 17/07/17: Changed to char flags
             tdis,nr_gam,de_gam,sp_gam,tmom2,nexc,delta, &
             deallocate_dis,i_sp,i_nr,i_de,nrnd,sp_fact,    &
!             nr_typ,idep,imar,de_gam1,krnd,ernd       
             de_gam1,krnd,Fdis,Fdis_deph,Fdis_rel,          &
             npulse,omega1,sigma1,tdelay,pshift,nrel            
             
!
      contains
!
      subroutine read_input

       !integer(i4b):: i,nspectra
       !character(3) :: medium,radiative,dissipative
       !character(5) :: dis_prop 
     
       !Molecular parameters 
       namelist /general/ n_ci_read,n_ci,mol_cc,n_f,medium,restart,full
       !External field paramaters
       namelist /field/ dt,n_step,n_out,Ffld,t_mid,sigma,omega,radiative,iseed,fmax, &
                        pulse,omega1,sigma1,tdelay,pshift
       !Namelist spectra
       namelist /spectra/ start,tau,dir_ft
       !Stochastic Schroedinger equation
       namelist /sse/ dissipative,idep,dis_prop,nrnd,tdis,nr_typ,krnd


       write(*,*) ''
       write(*,*) '*****************************************************'
       write(*,*) '*****************************************************'
       write(*,*) '**                                                 **' 
       write(*,*) '**                    WaveT                        **'
       write(*,*) '**     Evolution of molecular wavefunctions        **'
       write(*,*) '**   under external electromagnetic perturbations  **'
       write(*,*) '**                                                 **'
       write(*,*) '**                     by                          **'
       write(*,*) '**                Stefano Corni                    **'
       write(*,*) '**                Silvio Pipolo                    **'
       write(*,*) '**               Emanuele Coccia                   **'
       write(*,*) '**                 Gabriel Gil                     **'
       write(*,*) '**                Jacopo Fregoni                   **'
       write(*,*) '**                                                 **'
       write(*,*) '*****************************************************'
       write(*,*) '*****************************************************'    
       write(*,*) ''

       !Namelist mol 
       call init_nml_general()
       read(*,nml=general)
       call write_nml_general()

       !Namelist field
       call init_nml_field()
       read(*,nml=field)
       call write_nml_field() 

       !read gaussian output for CIS propagation
       call read_gau_out

       !Namelist spectra
       call init_nml_spectra()
       read(*,nml=spectra) 
       call write_nml_spectra() 
       !start, (tau(i),i=1,nspectra) !Start point for FT calculation &  Artificial damping 
       !(dir_ft(i),i=1,3) direction along which the field is oriented

       !Namelist sse
       call init_nml_sse() 
       read(*,nml=sse) 
       call write_nml_sse()
       if (Fdis(1:5).ne."nodis") call read_dis_params   

       return
      end subroutine
!
!
      subroutine read_gau_out
       integer(i4b) :: i,j
       character(4) :: junk
!    read ci energy
       open(7,file="ci_energy.inp",status="old")
       ! add ground state
       n_ci=n_ci+1
       n_ci_read=n_ci_read+1
       allocate (e_ci(n_ci))
       e_ci(1)=0.d0
       write (6,*) "Excitation energies read from input, in Hartree"
       do i=2,n_ci
         read(7,*) junk,junk,junk,e_ci(i)
         e_ci(i)=e_ci(i)*ev_to_au
         write(6,*) i-1,e_ci(i)
       enddo
       write (6,*) 
       close(7)
!    read transition dipoles (also for pcm: useful for analysis)
       open(7,file="ci_mut.inp",status="old")
       allocate (mut(3,n_ci,n_ci))
! SC 31/05/2016: now the GS data from gaussian is in debye and same format as others
!       read(7,*) junk,mut(1,1,1),mut(2,1,1),mut(3,1,1)
!       mut(:,1,1)=mut(:,1,1)*debye_to_au
       do i=1,n_ci_read
!         read(7,*) junk,mut(1,1,i),mut(2,1,i),mut(3,1,i)
        if (i.le.n_ci) then
         read(7,*)junk,junk,junk,junk,mut(1,1,i),mut(2,1,i),mut(3,1,i)
         mut(:,i,1)=mut(:,1,i)
        else
         read(7,*)
        endif
       enddo
       do i=2,n_ci_read
         do j=2,i   
          if (i.le.n_ci.and.j.le.n_ci) then
           read(7,*)junk,junk,junk,junk,mut(1,i,j),mut(2,i,j),mut(3,i,j)
           mut(:,j,i)=mut(:,i,j)
          else
           read(7,*)
          endif
         enddo
       enddo
!       write(6,*) "mut"
!       do i=1,n_ci
!        do j=1,n_ci
!         write(6,'(2i4,3f8.4)') i,j,mut(:,i,j)
!        enddo
!       enddo
       close(7)
!    read initial coefficients for the dynamics using the Slater determinants instead of the CIS_0 states.
       open(7,file="ci_ini.inp",status="old")
       allocate (c_i(n_ci))
       do i=1,n_ci
        read(7,*) c_i(i)
       enddo
       c_i=c_i/sqrt(dot_product(c_i,c_i))
       close(7)
       return
      end subroutine
!
      subroutine read_dis_params()
!------------------------------------------------------------------------
! @brief Read nonradiative and dephasing rates 
! Define spontaneous emission coefs from Einstein coefficients 
! Read phases randomly added during the propagation
!
! @date Created   : E. Coccia 21 Dec 2016
! Modified  :
!------------------------------------------------------------------------
     
       implicit none
        integer   :: i, j,k,idum, ierr0, ierr1, ierr2, ierr3, ierr4, err   
        real(dbl)  :: term   
        real(dbl)  :: rdum
        real(dbl)   :: rx,ry,rz
        real(dbl)   :: ix,iy,iz

       open(8,file='nr_rate.inp',status="old",iostat=ierr0,err=100)
       open(9,file='de_rate.inp',status="old",iostat=ierr1,err=101)
       open(10,file='de_phase.inp',status="old",iostat=ierr2,err=102)
       open(11,file='sp_rate.inp',status="old",iostat=ierr3,err=103)

       nexc = n_ci - 1
       if (Fful.eq.'Yesf') then
          nrel = nexc*(nexc-1)/2
       elseif (Fful.eq.'Nonf') then
          nrel = 0 
       endif
       nf=nexc+nrel 

       if (idep.ne.0.and.idep.ne.1) then
          write(*,*) 'Invalid value for idep, must be 0 or 1'
          stop
       endif

       if (nr_typ.ne.0.and.nr_typ.ne.1) then
          write(*,*) 'Invalid value for nr_typ, must be 0 or 1'
          stop
       endif 

       allocate(nr_gam(nf))
       if (idep.eq.0) then
          allocate(de_gam(nexc+1))
          allocate(delta(nexc+1))
       elseif (idep.eq.1) then
          allocate(de_gam(nexc))
          allocate(de_gam1(nexc+1))
       endif
       allocate(sp_gam(nf))
       allocate(tmom2(nf))
       allocate(sp_fact(nf))
       allocate(tomega(nf))

! Transition frequencies
       do i=1,nexc
          tomega(i) = e_ci(i+1)
       enddo
       if (Fful.eq.'Yesf') then 
          k=nexc
          do i=1,nexc
             do j=i+1,nexc
                k=k+1
                tomega(k) = abs(e_ci(i+1) - e_ci(j+1))
             enddo
          enddo
       endif
! Read dephasing terms
       do i=1,nexc
          read(9,*) idum, de_gam(i)
       enddo
       if (idep.eq.0) read(9,*) idum, de_gam(nexc+1)
! Read relaxation terms
       do i=1,nf
          read(8,*) idum, nr_gam(i)
          read(11,*) idum, sp_fact(i)
       enddo

! The sigma_z operator for dephasing has an extra factor 2
! If idep.eq.1 S_alpha = \sum_beta M(alpha, beta) |beta><beta|
! if alpha.eq. beta then M(alpha,beta) = -1
! otherwise M(alpha,beta) = 1
       if (idep.eq.1) then 
          de_gam=0.5d0*de_gam
          call define_gamma(de_gam,de_gam1,n_ci)
       endif
! Define spontaneous emission terms 
       term=4.d0/(3.d0*c**3)
       do i=1,nf
          sp_gam(i) = sp_fact(i)*term*tomega(i)**3
       enddo

! Define tmom2 =  <i|x|j>**2 + <i|y|j>**2 + <i|z|j>**2  
       do i=1,nexc
          tmom2(i) = mut(1,1,i+1)**2 + mut(2,1,i+1)**2 + mut(3,1,i+1)**2
       enddo
       if (Fful.eq.'Yesf') then 
          k=nexc
          do i=1,nexc
             do j=i+1,nexc 
                k=k+1
                tmom2(k) = mut(1,i+1,j+1)**2 + mut(2,i+1,j+1)**2 + mut(3,i+1,j+1)**2
             enddo
          enddo
       endif

       if (Fmdm.eq.'nan') then
          open(7,file="ci_mut_np.inp",status="old",iostat=ierr4,err=104)
          allocate(mut_np2(nf,3))
          do i=1,nf
             read(7,*) rdum, rx, ry, rz, ix, iy, iz, rdum
             mut_np2(i,1) = rx**2 + ix**2
             mut_np2(i,2) = ry**2 + iy**2
             mut_np2(i,3) = rz**2 + iz**2
          enddo
         close(7)
         do i=1,nf
            tmom2(i) = tmom2(i) + mut_np2(i,1) + mut_np2(i,2) + mut_np2(i,3)
         enddo

         write(*,*)"Contribution to the transition dipole <0|mu|i> ", &
                   "from NP"

104      if (ierr4.ne.0) then
           write(*,*)
           write(*,*) 'Dissipation and NP:'
           write(*,*) 'An error occurred during reading ci_mut_np.inp'
           write(*,*)
           stop
         endif
       endif
! Read phase randomly added during the propagation
       if (idep.eq.0) then
          do i=1,nexc+1
             read(10,*) idum, delta(i)
          enddo 
       endif

100    if (ierr0.ne.0) then
          write(*,*)
          write(*,*) 'Dissipation:'
          write(*,*) 'An error occurred during reading nr_rate.inp'
          write(*,*)
          stop
       endif  

101    if (ierr1.ne.0) then
          write(*,*)
          write(*,*) 'Dissipation:'
          write(*,*) 'An error occurred during reading de_rate.inp'
          write(*,*)
          stop
       endif  
     
102    if (ierr2.ne.0) then
          write(*,*)
          write(*,*) 'Dissipation:'
          write(*,*) 'An error occurred during reading de_phase.inp'
          write(*,*)
          stop
       endif

103    if (ierr3.ne.0) then
          write(*,*)
          write(*,*) 'Dissipation:'
          write(*,*) 'An error occurred during reading sp_rate.inp'
          write(*,*)
          stop
       endif
 
       close(8)
       close(9)
       close(10) 
       close(11)

       return

      end subroutine read_dis_params

      subroutine deallocate_dis()
!------------------------------------------------------------------------
! @brief Deallocate arrays for the dissipation 
!
! @date Created   : E. Coccia 22 Dec 2016
! Modified  :
!------------------------------------------------------------------------

       deallocate(nr_gam)
       deallocate(de_gam)
       if (idep.eq.1) deallocate(de_gam1)
       deallocate(sp_gam)
       deallocate(sp_fact)
       deallocate(tmom2)
       deallocate(tomega)
       if (idep.eq.0) deallocate(delta)
       if (Fdis.ne."nodis".and.Fmdm.eq.'nan') deallocate(mut_np2)

       return

      end subroutine deallocate_dis

      subroutine define_gamma(de_gam,de_gam1,nci)
!------------------------------------------------------------------------
! @brief Move from "physical" gammas to gammas 
! for keeping the population constant
! for each trajectory  
!
! @date Created   : E. Coccia 24 Mar 2017
! Modified  :
! @param de_gam(:), de_gam1(:)  
!------------------------------------------------------------------------

       implicit none
       integer                :: i
       integer, intent(in)    :: nci
       real(dbl), intent(in)    :: de_gam(nci-1)
       real(dbl), intent(inout) :: de_gam1(nci)

       de_gam1(1) = 0.d0
       do i=1,nexc
          de_gam1(i+1) = de_gam(i) - de_gam1(1) 
       enddo

       return

      end subroutine define_gamma

      subroutine init_nml_general()
!------------------------------------------------------------------------
! @brief Initialize variables in the namelist general 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param n_ci_read,n_ci,mol_cc,n_f,medium 
!------------------------------------------------------------------------

       ! Molecular center of charge
       mol_cc=0.d0
       ! Integer to be appended in all output files
       n_f=1
       ! Vacuum calculation
       medium='vac'
       ! Restart
       restart='n'
       ! Full or only |e> -> |0> relaxation
       full='n'

       return

      end subroutine init_nml_general

      subroutine init_nml_field()
!------------------------------------------------------------------------
! @brief Initialize variables in the namelist field 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param dt,n_step,n_out,Ffld,t_mid,sigma,omega,radiative,iseed,fmax 
!------------------------------------------------------------------------

       ! Time step in the propagation (a.u.)
       dt=0.05
       ! Number of steps
       n_step=10000
       ! Frequency in writing files
       n_out=1
       ! Type of field
       Ffld='gau'
       ! Center of the pulse
       t_mid=200
       ! Sigma for the pulse
       sigma=10 
       ! Frequency (no oscillation)
       omega=0.d0
       ! Stochastic field
       radiative='non'
       ! ieed for radiative and/or dissipation
       iseed=12345678
       ! Amplitude of the external field
       fmax=0.d0
       ! Number of pulses
       pulse='one'
       ! Time delay
       tdelay=2000.
       ! Phase shift
       pshift=0.d0
       ! Frequency second pulse
       omega1=0.d0
       ! Sigma for the second pulse
       sigma1=0.d0

       return

      end subroutine init_nml_field

      subroutine init_nml_spectra()
!------------------------------------------------------------------------
! @brief Initialize variables in the namelist spectra 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param start,tau,dir_ft 
!------------------------------------------------------------------------

       ! Parameter for computing spectra
       nspectra=1
       ! Parameter for computing spectra
       tau(:)=zero
       ! Direction fo the field (no field)
       dir_ft=0.d0

       return

      end subroutine init_nml_spectra

      subroutine init_nml_sse()
!------------------------------------------------------------------------
! @brief Initialize variables in the namelist sse 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param dissipative,idep,dis_prop,nrnd,tdis,nr_typ,krnd 
!------------------------------------------------------------------------

       ! Do dissipation 
       dissipative='non'
       ! Type of the dephasing operator
       idep=1
       ! Type of propagator
       dis_prop='qjump'
       ! If dis_prop='euler', steps for accumulating the Wiener process
       nrnd=1
       ! If dis_prop='euler', use the Euler-Maruyama algorithm
       tdis=0
       ! Type of nonradiative relaxation
       nr_typ=1
       ! Factor is dissipation='ernd'
       krnd=1.d0

       return

      end subroutine init_nml_sse

      subroutine write_nml_general()
!------------------------------------------------------------------------
! @brief Write variables in the namelist general and put conditions 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param n_ci_read,n_ci,mol_cc,n_f,medium 
!------------------------------------------------------------------------

       if (n_ci.le.0) then
          write(*,*) 'ERROR: number of excited states is wrong', n_ci
       endif 
       if (n_ci_read.le.0) then
          write(*,*) 'ERROR: number of excited states is wrong', n_ci
       endif
       write(6,*) "Number to append to dat file", n_f
       write (*,*) "Number of CIS states to be read",n_ci_read
       write (*,*) "Number of CIS states to be used",n_ci
       write(*,*) 'Molecular center of charge'
       write(*,*) mol_cc 
       write(*,*) ''
       select case (medium)
        case ('sol','Sol','SOL')
          write(*,*) "Solvent as external medium"
          Fmdm='Csol'
        case ('qso','Qso','QSO')
          write(*,*) "Quantum Solvent as external medium"
          Fmdm='Csol'
        case ('nan','Nan','NAN')
          write(*,*) "Nanoparticle as external medium"
          Fmdm='Cnan'
        case ('Qna','qna','QNA')
          write(*,*) "Quantum Nanoparticle as external medium"
          Fmdm='Qnan'
        case default
          write(*,*) "No external medium, vacuum calculation"
          Fmdm='vac'
       end select
       select case (restart)
        case ('n', 'N')
          write(*,*) 'No restart'
          Fres='Nonr'
        case ('y', 'Y')
          write(*,*) 'Restart from a previous calculation'
          Fres='Yesr'
       end select
       select case (full)
        case ('n', 'N') 
          write(*,*) 'Only |k> -> |0> relaxation'
          Fful='Nonf'
        case ('y', 'Y')
          write(*,*) 'Full relaxation'
          Fful='Yesf'
       end select
       write(*,*) ''

       return

      end subroutine write_nml_general

      subroutine write_nml_field()
!------------------------------------------------------------------------
! @brief Write variables in the namelist field and put conditions 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param dt,n_step,n_out,Ffld,t_mid,sigma,omega,radiative,iseed,fmax 
!------------------------------------------------------------------------

       write (*,*) "Time step (in au), number of steps, stride",dt, &
                 n_step,n_out
       write (*,*) "Time shape of the perturbing field",Ffld
       write (*,*) "time at the center of the pulse (au):",t_mid
       write (*,*) "Width of the pulse (time au):",sigma
       write (*,*) "Frequency (au):",omega
       write (*,*) "Maximum E field",fmax
       !SC
       select case (radiative)
        case ('rad','Rad','RAD')
         Frad='arl'
        case default
         Frad='non'
       end select
       select case (pulse)
        case ('one', 'ONE', 'One')
          write(*,*) 'One pulse'
          npulse=1
        case ('two', 'TWO', 'Two')
          write(*,*) 'Two pulses'
          write(*,*) 'Time delay =', tdelay
          write(*,*) 'Phase shift =', pshift
          write(*,*) 'Frequency second pulse (au) =', omega1
          write(*,*) 'Sigma second pulse =', sigma1
          npulse=2
        case default
          npulse=1
       end select
       write(*,*) ''

       return

      end subroutine write_nml_field

      subroutine write_nml_spectra()
!------------------------------------------------------------------------
! @brief Write variables in the namelist spectra and put conditions 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param start,tau,dir_ft 
!------------------------------------------------------------------------
 
       if(medium.ne.'vac') nspectra=2

       write(*,*) 'Starting point for FT calculation', start
       write(*,*) 'Artificial damping', (tau(i),i=1,nspectra)
       write(*,*) 'Direction along which the field is oriented', dir_ft
       write(*,*) ''

       return

      end subroutine write_nml_spectra

      subroutine write_nml_sse()
!------------------------------------------------------------------------
! @brief Write variables in the namelist sse and put conditions 
!
! @date Created   : E. Coccia 11 May 2017
! Modified  :
! @param dissipative,idep,dis_prop,nrnd,tdis,nr_typ,krnd 
!------------------------------------------------------------------------

       select case (dissipative)
        case ('mar', 'Mar', 'MAR')
          !dis=.true.
          !imar=0
          select case (idep)
           case (0)
            Fdis_deph="exp"
            write(*,*) 'Use sqrt(gamma_i) exp(i delta_i)|i><i|(i=0,nci)'
           case (1)
            Fdis_deph="i-0"
            write(*,*) 'Use sqrt(gamma_i) (|i><i|-|0><0|)'
          end select
          write(*,*) 'Markovian dissipation'
          select case (dis_prop)
           case ('qjump', 'Qjump', 'QJump')
             !qjump=.true. 
             Fdis="mar-qjump"
             write(*,*) 'Quantum jump algorithm'
           case ('Euler', 'EUler', 'EULER', 'euler')
             !qjump=.false.
             write(*,*) 'Continuous stochastic propagator'
             write(*,*) 'Time step for the Brownian motion is:', dt/nrnd
             select case (tdis)
              case (0)
                write(*,*) 'Euler-Maruyama algorithm'
                Fdis="mar-EuMar"
              case (1)
                write(*,*) 'Leimkuhler-Matthews algorithm'
                Fdis="mar-LeiMa"
             end select
          end select
          select case (nr_typ)
           case (0)
            Fdis_rel="dip"
            write(*,*) 'Internal conversion relaxation via dipole'
           case (1)
            Fdis_rel="mat"
            write(*,*) 'Internal conversion relaxation via given matrix'
          end select
        case ('nma', 'NMa', 'NMA', 'Nma')
          !dis=.true.
          !imar=1
          write(*,*) 'NonMarkovian dissipation'
          read(*,*) dis_prop
          select case (dis_prop)
           case ('qjump', 'Qjump', 'QJump')
             !qjump=.true.
             Fdis="nma-qjump"
             write(*,*) 'Quantum jump algorithm'
           case default
             !qjump=.false.
             Fdis="nma-cstoc"
             write(*,*) 'Continuous stochastic propagator'
          end select
        case ('rnd', 'RND', 'Rnd', 'RNd')
          !dis=.false.
          !ernd=.true.
          Fdis="ernd"
          write(*,*) "Normal random term added to CI energies"
        case default
          !dis=.false.
          Fdis="nodis"
          write(*,*) 'No dissipation'
       end select
       write(*,*) ''

       return

      end subroutine write_nml_sse


    end module readio
