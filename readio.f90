      module readio
      use cav_types       
      implicit none
      save
!
      real(8), parameter :: pi=3.141592653589793D+00
      real(8), parameter :: ev_to_au=0.0367493
      real(8), parameter :: debye_to_au=0.393456
      real(8), parameter :: slight=137.d0
      real(dbl), parameter :: zero=0.d0
      real(dbl), parameter :: one=1.d0
      real(dbl), parameter :: two=2.d0
      real(dbl), parameter :: twp=two*pi
      real(dbl), parameter :: pt5=0.5d0
      integer(i4b), parameter :: one_i=1
      complex(cmp), parameter :: onec=(one,zero)                
      complex(cmp), parameter :: twoc=(two,zero)                
      complex(16), parameter :: ui=(zero,one)
      integer(4) :: n_f,n_ci,n_ci_read,n_step,n_out,tdis !tdis=0 Euler, tdis=1 Matthews 
      integer(4) :: imar !imar=0 Markvian, imar=1, nonMarkovian
      integer(4) :: i_sp=0,i_nr=0,i_de=0 !counters for quantum jump occurrences
      integer(4) :: nrnd !the time step for Euler-Maruyama is dt/nrnd
      integer(4) :: nr_typ !type of decay for the internal conversion
      integer(4) :: idep !choose the dephasing operator
      real(8), allocatable :: c_i(:),e_ci(:)  ! coeff and energy from cis
      real(8), allocatable :: mut(:,:,:) !transition dipoles from cis
      real(8), allocatable :: nr_gam(:), de_gam(:) !decay rates for nonradiative and dephasing events
      real(8), allocatable :: sp_gam(:) !decay rate for spontaneous emission  
      real(8), allocatable :: sp_fact(:) !multiplicative factor for the decay rate for spontaneous emission
      real(8), allocatable :: tmom2_0i(:) !square transition moments i->0
      real(8), allocatable :: delta(:) !phases randomly added during the propagation 
      real(8), allocatable :: de_gam1(:) !combined decay rates for idep=1
      real(8) :: dt,tau(2),start, krnd
      logical :: dis !turns on the dissipation
      logical :: qjump ! =.true. quantum jump, =.false. stochastic propagation
      logical :: ernd=.false. !add normal number to E: E -> E + krnd*rnd()
      ! qjump works for the Markovian case
      ! Stochastic propagation from:
      ! Appl. Math. Res. Express vol. 2013 34-56 (2013)
      ! IMA J. Numer. Anal. vol. 36 13-79 (2016)
! SP 28/10/16: improved with a vector
!      integer(4) :: dir_ft
      real(8) :: t_mid,sigma,dir_ft(3),fmax(3),omega,mol_cc(3)
      character(3) :: mdm,Ffld,rad 
      integer(4) :: iseed  !seed for random number generator
      integer(4) :: nexc !number of excited states
! kind of surrounding medium and shape of the impulse
!     mdm=sol: solvent
!     mdm=nan: nanoparticle
!     mdm=vac: no medium
!
!     Ffld=gau: gaussian impulse
!SC
!     rad: wheter or not to apply radiative damping
!     TO BE COMPLETED, SEE PROPAGATE.F90      
      real(8) :: eps_A,eps_gm,eps_w0,f_vel

      
      private
      public read_input,n_ci,n_ci_read,n_step,dt, &
             Ffld,t_mid,sigma,omega,fmax, & 
             mdm,mol_cc,tau,start,c_i,e_ci,mut,ui,pi,zero,one,two,twp,&
             one_i,onec,twoc,pt5,rad,n_out,iseed,n_f,dir_ft, &
             dis,tdis,nr_gam,de_gam,sp_gam,tmom2_0i,nexc,delta, &
             deallocate_dis,qjump,i_sp,i_nr,i_de,nrnd,sp_fact,nr_typ, &
             idep, imar, de_gam1, krnd, ernd
!
      contains
!
      subroutine read_input
       integer(4):: i,nspectra
       character(3) :: medium,radiative,dissipative
       character(5) :: dis_prop 
       read(5,*)    
       read(5,*) n_f
       write(6,*) "Number to append to dat file", n_f
       read(5,*) n_ci_read,n_ci
       write (6,*) "Number of CIS states to be read",n_ci_read
       write (6,*) "Number of CIS states to be used",n_ci
       read(5,*) dt,n_step,n_out
       write (6,*) "Time step (in au), number of steps, stride",dt, &
                 n_step,n_out
       read(5,*) Ffld             
       write (6,*) "Time shape of the perturbing field",Ffld
       read(5,*) t_mid,sigma,omega
       write (6,*) "time at the center of the pulse (au):",t_mid
       write (6,*) "Width of the pulse (time au):",sigma
       write (6,*) "Frequency (au):",omega
       read(5,*) fmax
       write (6,*) "Maximum E field",fmax
!SC
       read(5,*) radiative
       select case (radiative)
        case ('rad','Rad','RAD')
         rad='arl'
         read(5,*) iseed
        case default
         rad='non'
       end select 
       ! Moleculalr center of charge: used?
       read(5,*) (mol_cc(i),i=1,3)
       write(*,*) 'Molecular center of charge'
       write(*,*) mol_cc
!    read gaussian output for CIS propagation
       call read_gau_out
!    read external medium type: sol, nan, vac
       read(5,*) 
       read(5,*) medium
       select case (medium)
        case ('sol','Sol','SOL')
          write(*,*) "Solvent as external medium" 
          mdm='sol'
        case ('nan','Nan','NAN')
          write(*,*) "Nanoparticle as external medium" 
          mdm='nan'
        case default
          write(*,*) "No external medium, vacuum calculation" 
          mdm='vac'
       end select
!    read spectra calculation parameters:
       nspectra=1
       tau(:)=zero
       if(medium.ne.'vac') nspectra=2 
       read(5,*) 
       read(5,*) start, (tau(i),i=1,nspectra) !Start point for FT calculation &  Artificial damping 
       read(5,*) (dir_ft(i),i=1,3)      ! direction along which the field is oriented
!    dissipation using SSE
       read(5,*)  
       read(5,*) dissipative
       select case (dissipative)
        case ('mar', 'Mar', 'MAR')
          dis=.true.
          imar=0
          read(*,*) idep
          select case (idep)
           case (0)
            write(*,*) 'Use sqrt(gamma_i) exp(i delta_i)|i><i|(i=0,nci)'
           case (1)
            write(*,*) 'Use sqrt(gamma_i) (|i><i|-|0><0|)'
          end select
          write(*,*) 'Markovian dissipation'
          if (rad.eq.'non') read(5,*) iseed
          read(*,*) dis_prop
          select case (dis_prop)
           case ('qjump', 'Qjump', 'QJump')
             qjump=.true.
             write(*,*) 'Quantum jump algorithm'
           case ('Euler', 'EUler', 'EULER', 'euler')
             qjump=.false.
             read(*,*) nrnd
             write(*,*) 'Continuous stochastic propagator'
             write(*,*) 'Time step for the Brownian motion is:', dt/nrnd
             read(*,*) tdis
             select case (tdis)
              case (0)
                write(*,*) 'Euler-Maruyama algorithm'
              case (1)
                write(*,*) 'Leimkuhler-Matthews algorithm'
             end select
          end select
          read(*,*) nr_typ
          select case (nr_typ)
           case (0)
            write(*,*) 'Internal conversion relaxation via dipole'
           case (1)
            write(*,*) 'Internal conversion relaxation via given matrix'
          end select
        case ('nma', 'NMa', 'NMA', 'Nma')
          dis=.true.
          imar=1
          write(*,*) 'NonMarkovian dissipation'
          if (rad.eq.'non') read(5,*) iseed
          read(*,*) dis_prop
          select case (dis_prop)
           case ('qjump', 'Qjump', 'QJump')
             qjump=.true.
             write(*,*) 'Quantum jump algorithm'
           case default
             qjump=.false.
             write(*,*) 'Continuous stochastic propagator'
          end select
        case ('rnd', 'RND', 'Rnd', 'RNd')
          dis=.false.
          ernd=.true.
          read(*,*) krnd
          write(*,*) "Normal random term added to CI energies"
          if (rad.eq.'non') read(5,*) iseed
        case default
          dis=.false.
          write(*,*) 'No dissipation'    
       end select
       if (dis) call read_dis_params   

       return
      end subroutine
!
!
      subroutine read_gau_out
       integer(4) :: i,j
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
! SC 16/02/2016: for efficiency reasons, it would be better
!   to define mut(3,n_ci,n_ci)
       allocate (mut(n_ci,n_ci,3))
! SC 31/05/2016: now the GS data from gaussian is in debye and same format as others
!       read(7,*) junk,mut(1,1,1),mut(1,1,2),mut(1,1,3)
!       mut(1,1,:)=mut(1,1,:)*debye_to_au
       do i=1,n_ci_read
!         read(7,*) junk,mut(1,i,1),mut(1,i,2),mut(1,i,3)
        if (i.le.n_ci) then
         read(7,*)junk,junk,junk,junk,mut(1,i,1),mut(1,i,2),mut(1,i,3)
         mut(i,1,:)=mut(1,i,:)
        else
         read(7,*)
        endif
       enddo
       do i=2,n_ci_read
         do j=2,i   
          if (i.le.n_ci.and.j.le.n_ci) then
           read(7,*)junk,junk,junk,junk,mut(i,j,1),mut(i,j,2),mut(i,j,3)
           mut(j,i,:)=mut(i,j,:)
          else
           read(7,*)
          endif
         enddo
       enddo
!       write(6,*) "mut"
!       do i=1,n_ci
!        do j=1,n_ci
!         write(6,'(2i4,3f8.4)') i,j,mut(i,j,:)
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
! Read nonradiative and dephasing rates 
! Define spontaneous emission coefs from Einstein coefficients 
! Read phases randomly added during the propagation
!
! Created   : E. Coccia 21 Dec 2016
! Modified  :
!------------------------------------------------------------------------
     
       implicit none
       integer   :: i, idum, ierr0, ierr1, ierr2, ierr3, err   
       real(8)  :: term  

       open(8,file='nr_rate.inp',status="old",iostat=ierr0,err=100)
       open(9,file='de_rate.inp',status="old",iostat=ierr1,err=101)
       open(10,file='de_phase.inp',status="old",iostat=ierr2,err=102)
       open(11,file='sp_rate.inp',status="old",iostat=ierr3,err=103)

       nexc = n_ci - 1

       if (idep.ne.0.and.idep.ne.1) then
          write(*,*) 'Invalid value for idep, must be 0 or 1'
          stop
       endif

       if (nr_typ.ne.0.and.nr_typ.ne.1) then
          write(*,*) 'Invalid value for nr_typ, must be 0 or 1'
          stop
       endif 

       allocate(nr_gam(nexc))
       if (idep.eq.0) then
          allocate(de_gam(nexc+1))
          allocate(delta(nexc+1))
       elseif (idep.eq.1) then
          allocate(de_gam(nexc))
          allocate(de_gam1(nexc+1))
       endif
       allocate(sp_gam(nexc))
       allocate(tmom2_0i(nexc))
       allocate(sp_fact(nexc))

! Read nonradiative and dephasing terms
       do i=1,nexc
          read(8,*) idum, nr_gam(i)
          read(9,*) idum, de_gam(i)
          read(11,*) idum, sp_fact(i)
       enddo
       if (idep.eq.0) read(9,*) idum, de_gam(nexc+1)
! The sigma_z operator for dephasing has an extra factor 2
! If idep.eq.1 S_alpha = \sum_beta M(alpha, beta) |beta><beta|
! if alpha.eq. beta then M(alpha,beta) = -1
! otherwise M(alpha,beta) = 1
       if (idep.eq.1) then 
          de_gam=0.5d0*de_gam
          call define_gamma(de_gam,de_gam1,n_ci)
       endif
! Define spontaneous emission terms 
       term=4.d0/(3.d0*slight)
       do i=1,nexc
          sp_gam(i) = sp_fact(i)*term*e_ci(i+1)**3
       enddo

! Define tmom2_0i =  <i|x|0>**2 + <i|y|0>**2 + <i|z|0>**2  
       do i=1,nexc 
          tmom2_0i(i) = mut(1,i+1,1)**2 + mut(1,i+1,2)**2 + mut(1,i+1,3)**2
       enddo

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
! Deallocate arrays for the dissipation 
!
! Created   : E. Coccia 22 Dec 2016
! Modified  :
!------------------------------------------------------------------------

       deallocate(nr_gam)
       deallocate(de_gam)
       if (idep.eq.1) deallocate(de_gam1)
       deallocate(sp_gam)
       deallocate(sp_fact)
       deallocate(tmom2_0i)
       if (idep.eq.0) deallocate(delta)

       return

      end subroutine deallocate_dis

      subroutine define_gamma(de_gam,de_gam1,nci)
!------------------------------------------------------------------------
! Move from "physical" gammas to gammas 
! for keeping the population constant
! for each trajectory  
!
! Created   : E. Coccia 24 Mar 2017
! Modified  :
!------------------------------------------------------------------------

       implicit none
       integer                :: i
       integer, intent(in)    :: nci
       real(8), intent(in)    :: de_gam(nci-1)
       real(8), intent(inout) :: de_gam1(nci)

       de_gam1(1) = 0.d0
       do i=1,nexc
          de_gam1(i+1) = de_gam(i) - de_gam1(1) 
       enddo

       return

      end subroutine define_gamma

    end module readio
