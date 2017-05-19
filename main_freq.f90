      program tdcis
      use readio  
      use readio_medium
      use td_contmed
      implicit none
      integer :: st,current,rate
      integer(4) :: n_omega,i
      real(8),allocatable :: omega_list(:)
      complex(16),allocatable :: eps_list(:)
      real(8) :: omega_ini,omega_end
      character(3) :: mode
!
!     read in the input parameter for the present evolution
      call system_clock(st,rate)
      read(5,*) mode
      select case(mode)
       case("Ext","EXT","ext")
        Ffreq="ext"
        call read_medium_freq
        read(5,*) n_omega,omega_ini,omega_end
        allocate(omega_list(n_omega))
        allocate(eps_list(n_omega))
        do i=1,n_omega
         omega_list(i)=(omega_end-omega_ini)/(n_omega-1)*(i-1)+omega_ini
         eps_list(i)=1.d0-eps_A/(omega_list(i)*(omega_list(i)+ui*eps_gm))
!         write(6,*) 'omega',omega_list(i)
        enddo
       case("Mol","mol","MOL")
         Ffreq="dip"
         read(5,*) eps_sol
         write(6,*) "Dynamical solvent constant",eps_sol
         read(5,*) n_omega
         n_ci=n_omega+1
         allocate(omega_list(n_omega))
         allocate(eps_list(n_omega))
         write(6,*) "frequency and eps read from input"
         do i=1,n_omega
          read(5,*) omega_list(i),eps_list(i)
          write(6,'(3d12.5)') omega_list(i),eps_list(i)
         enddo
         read(5,*) 
         read(5,*) mol_cc 
         allocate (mut(n_ci,1,3))
         write(6,*) "transition dipole read from input"
         do i=1,n_omega
          read(5,*) mut(i,1,:) 
          write(6,'(i6,3d12.5)') i,mut(i,1,:)
         enddo
! end of debug
       case("Gam","gam","GAM")
         Ffreq="gam"
         read(5,*) eps_sol
         write(6,*) "Dynamical solvent constant",eps_sol
         read(5,*) n_omega
         n_ci=n_omega+1
         allocate(omega_list(n_omega))
         allocate(eps_list(n_omega))
         write(6,*) "frequency and eps read from input"
         do i=1,n_omega
          read(5,*) omega_list(i),eps_list(i)
          write(6,'(3d12.5)') omega_list(i),eps_list(i)
         enddo
         call read_trans_out_medium
        case default
         write(6,*) "Pls select external field ext", &
               " or molecular dipole moli or potential from gamess gam"
         stop
        end select
      call system_clock(current)
      write(6,'("Done reading input, took", &
            F10.3,"s")') real(current-st)/real(rate)
!
!     diagonalise matrix
      call do_freq_mat(omega_list,eps_list,n_omega)    
      call system_clock(current)
      write(6,'("Done , total elapsed time", &
            F10.3,"s")') real(current-st)/real(rate)
!         
      deallocate(omega_list,eps_list)
      stop
      end
