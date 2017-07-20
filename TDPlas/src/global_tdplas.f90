      module global_tdplas
      implicit none
      INTEGER, PARAMETER :: dbl = selected_real_kind(14,200)
      INTEGER, PARAMETER :: sgl = selected_real_kind(6,30)
      INTEGER, PARAMETER :: i4b = selected_int_kind(9)
      INTEGER, PARAMETER :: cmp = dbl
!
      !LIGHT SPEED IN AU TO BE REVISED
      ! constants, conversion factors, and number literals
      real(8), parameter :: pi=3.141592653589793D+00
      real(8), parameter :: ev_to_au=0.0367493
      real(8), parameter :: debye_to_au=0.393456
      real(dbl), parameter :: TOANGS=0.52917724924D+00
      real(dbl), parameter :: ANTOAU=1.0D+00/TOANGS
      real(dbl), parameter :: c=1.37036d2
      real(8), parameter :: slight=137.d0
      real(dbl), parameter :: zero=0.d0
      real(dbl), parameter :: one=1.d0
      real(dbl), parameter :: two=2.d0
      real(dbl), parameter :: twp=two*pi
      real(dbl), parameter :: pt5=0.5d0
      integer(i4b), parameter :: one_i=1
      double complex, parameter :: zeroc=(zero,zero)      
      double complex, parameter :: onec=(one,zero)                
      double complex, parameter :: twoc=(two,zero)                
      double complex, parameter :: ui=(zero,one)

      ! block of global variables to be supplied by WaveT
      real(8) 		:: dt				! time step
      character(3) 	:: mdm				! kind of medium
      integer(4) 	:: n_ci,n_ci_read		! number of CIS states
      real(8), allocatable :: c_i(:),e_ci(:)		! CIS coefficients and energies
      real(8), allocatable :: mut(:,:,:)			! CIS transition dipoles
      real(8) 		:: mol_cc(3)			! molecule center
      real(8) 		:: fmax(3),omega			! field amplitude and frequency
      character(3) 	:: Ffld				! shape of impulse
      integer(4) 	:: n_out,n_f			! auxiliaries for output

      contains
  
      subroutine set_global_tdplas(this_dt,this_mdm,this_mol_cc,this_n_ci,this_n_ci_read,this_c_i,this_e_ci,this_mut,&
				   this_fmax,this_omega,this_Ffld,this_n_out,this_n_f)
        implicit none
        real(8)     , intent(in) :: this_dt				! time step
        character(3), intent(in) :: this_mdm				! kind of medium
        integer(4)  , intent(in) :: this_n_ci,this_n_ci_read		! number of CIS states
        real(8)     , intent(in) :: this_c_i(:),this_e_ci(:)		! CIS coefficients and energies
        real(8)     , intent(in) :: this_mut(:,:,:)			! CIS transition dipoles
        real(8)     , intent(in) :: this_mol_cc(3)			! molecule center
        real(8)     , intent(in) :: this_fmax(3),this_omega		! field amplitude and frequency
        character(3), intent(in) :: this_Ffld				! shape of impulse
        integer(4)  , intent(in) :: this_n_out,this_n_f			! auxiliaries for output
        
        dt=this_dt
        mdm=this_mdm
        mol_cc=this_mol_cc
        n_ci=this_n_ci
        n_ci_read=this_n_ci_read
        allocate(c_i(n_ci),e_ci(n_ci))
        allocate(mut(n_ci,n_ci,3))
        c_i=this_c_i
        e_ci=this_e_ci
        mut=this_mut
        fmax=this_fmax
        omega=this_omega
        Ffld=this_Ffld
        n_out=this_n_out
        n_f=this_n_f
        return
      end subroutine set_global_tdplas
      end module global_tdplas
