      module global_wavet 
      use constants
      use td_contmed, only: set_charges,get_mdm_dip,get_gneq,init_mdm, &
                            prop_mdm,finalize_mdm
      use readio_medium, only: q0,Fmdm_relax,read_medium

      public set_q0charges,get_medium_dip,read_medium_input

      contains
  
      ! SP 270917: bridge subroutine to set charges qr_t to q0 during propagation, probably bad idea...
      subroutine set_q0charges
        implicit none
        call set_charges(q0)
        return
      end subroutine set_q0charges
      !
      ! SP 270917: small routine to set the dipole(t) in Sdip for spectra
      subroutine get_medium_dip(mdm_dip)
        implicit none
        real(dbl), intent(inout) :: mdm_dip(3)
        call get_mdm_dip(mdm_dip)
        return
      end subroutine get_medium_dip
      !
      ! Read medium input 
      subroutine read_medium_input
        implicit none
        call read_medium
        return
      end subroutine read_medium_input
      !
      ! get energies
      subroutine get_energies(e_vac,g_eq_t,g_neq_t,g_neq2_t)     
        implicit none
        real(dbl), intent(inout) :: e_vac,g_neq_t,g_neq2_t,g_eq_t
        call get_gneq(e_vac,g_eq_t,g_neq_t,g_neq2_t)
        return
      end subroutine get_energies
      !
      ! init medium
      subroutine init_medium(c,f,h,cm)     
        implicit none
        complex(cmp), intent(inout) :: c(:),cm(:),h(:,:)
        real(dbl), intent(inout) :: f(3)
        call init_mdm(c,f,h)
        return
      end subroutine init_medium
      !
      ! propagate medium
      subroutine prop_medium(i,c,f,h,cm)     
        implicit none
        complex(cmp), intent(inout) :: c(:),cm(:),h(:,:)
        real(dbl), intent(inout) :: f(3)
        integer(i4b), intent(inout) :: i
        call prop_mdm(i,c,f,h)
        return
      end subroutine prop_medium
      !
      ! deallocate medium
      subroutine finalize_medium
        implicit none
        call finalize_mdm
        return
      end subroutine finalize_medium
      

      end module global_wavet 
