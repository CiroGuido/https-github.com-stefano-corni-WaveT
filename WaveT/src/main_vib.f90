program main_vib
 use readio
 use vib

 implicit none

  call read_input_vib()
  call read_e_dip()

  call compute_e_dip()

  call write_e_dip()

  ! Vibronic coupling for nonradiative relaxation
  if (coupling) call compute_coupling()
 
  call deallocate_vib()
 
  stop

end program main_vib
