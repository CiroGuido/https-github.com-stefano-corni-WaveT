program main_vib
  use readio
  use vib

  implicit none

  ! Input for vibrational levels
  call read_input_vib()
  ! Read electronic energies and dipoles
  call read_e_dip()
  ! Compute FC factors
  call compute_e_dip()
  write(*,*) 'ciao'
  ! Vibronic coupling for nonradiative relaxation
  if (coupling) call compute_coupling()
  ! Deallocate arrays 
  call deallocate_vib()
 
  stop

end program main_vib
