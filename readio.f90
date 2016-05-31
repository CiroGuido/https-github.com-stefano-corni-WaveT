      module readio
      use cav_types       
      implicit none
      save
!
      real(8), parameter :: pi=3.141592653589793D+00
      real(8), parameter :: ev_to_au=0.0367493
      real(8), parameter :: debye_to_au=0.393456
      real(dbl), parameter :: zero=0.d0
      real(dbl), parameter :: one=1.d0
      real(dbl), parameter :: two=2.d0
      real(dbl), parameter :: twp=two*pi
      real(dbl), parameter :: pt5=0.5d0
      integer(i4b), parameter :: one_i=1
      complex(cmp), parameter :: onec=(one,zero)                
      complex(cmp), parameter :: twoc=(two,zero)                
      complex(16), parameter :: ui=(zero,one)
      integer(4) :: n_ci,n_step
      real(8), allocatable :: c_i(:),e_ci(:)  ! coeff and energy from cis
      real(8), allocatable :: mut(:,:,:) !transition dipoles from cis
      real(8) :: dt,tau,start
      real(8) :: t_mid,sigma,fmax(3),omega,mol_cc(3)
      real(8), allocatable :: q0(:)
      character(3) :: mdm,tfield 
! kind of surrounding medium and shape of the impulse
!     mdm=sol: solvent
!     mdm=nan: nanoparticle
!     mdm=vac: no medium
!
!     tfield=gau: gaussian impulse
!     TO BE COMPLETED, SEE PROPAGATE.F90      
      real(8) :: eps_A,eps_gm,eps_w0,f_vel

      
      private
      public read_input,n_ci,n_step,dt,tfield,t_mid,sigma,omega,fmax, & 
             mdm,mol_cc,tau,start,c_i,e_ci,mut,ui,pi,zero,one,two,twp,&
             one_i,onec,twoc,pt5
!
      contains
!
      subroutine read_input
       integer(4):: i
       character(3) :: medium
       read(5,*) n_ci
       write (6,*) "Number of CIS states",n_ci
       read(5,*) dt,n_step
       write (6,*) "Time step (in au) and number of steps",dt,n_step
       read(5,*) tfield             
       write (6,*) "Time shape of the perturbing field",tfield
       read(5,*) t_mid,sigma,omega
       write (6,*) "time at the center of the pulse (au):",t_mid
       write (6,*) "Width of the pulse (time au):",sigma
       write (6,*) "Frequency (au):",omega
       read(5,*) fmax
!    read external medium type: sol, nan, vac
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
       ! Moleculalr center of charge: used?
       read(5,*) (mol_cc(i),i=1,3)
       write(*,*) 'Molecular center of charge'
       write(*,*) mol_cc
!    read spectra calculation parameters:
       read(5,*) tau       ! Artificial damping 
       read(5,*) start     ! Start point for FT calculation
!    read gaussian output for CIS propagation
       call read_gau_out
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
       allocate (e_ci(n_ci))
       e_ci(1)=0.d0
       do i=2,n_ci
         read(7,*) junk,junk,junk,e_ci(i)
         e_ci(i)=e_ci(i)*ev_to_au
       enddo
       close(7)
!    read transition dipoles (also for pcm: useful for analysis)
       open(7,file="ci_mut.inp",status="old")
! SC 16/02/2016: for efficiency reasons, it would be better
!   to define mut(3,n_ci,n_ci)
       allocate (mut(n_ci,n_ci,3))
! SC 31/05/2016: now the GS data from gaussian is in debye and same format as others
!       read(7,*) junk,mut(1,1,1),mut(1,1,2),mut(1,1,3)
!       mut(1,1,:)=mut(1,1,:)*debye_to_au
       do i=1,n_ci
!         read(7,*) junk,mut(1,i,1),mut(1,i,2),mut(1,i,3)
        read(7,*)junk,junk,junk,junk,mut(1,i,1),mut(1,i,2),mut(1,i,3)
        mut(i,1,:)=mut(1,i,:)
       enddo
       do i=2,n_ci
         do j=2,i   
           read(7,*)junk,junk,junk,junk,mut(i,j,1),mut(i,j,2),mut(i,j,3)
           mut(j,i,:)=mut(i,j,:)
         enddo
!! Gaussian subtract the dipole moment for the ground state
!         mut(i,i,:)=mut(i,i,:)+mut(1,1,:)
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
      end module
