      module cav_types
      implicit none
      INTEGER, PARAMETER :: dbl = selected_real_kind(14,200)
      INTEGER, PARAMETER :: sgl = selected_real_kind(6,30)
      INTEGER, PARAMETER :: i4b = selected_int_kind(9)
      INTEGER, PARAMETER :: cmp = dbl
!
      type tessera
       real(dbl) :: z
       real(dbl) :: phi
       real(dbl) :: fz
       real(dbl) :: dfdz
       real(dbl) :: dz
       real(dbl) :: area
      end type
!
      type tess_pcm
       real(dbl) :: x
       real(dbl) :: y
       real(dbl) :: z
       real(dbl) :: area
       real(dbl) :: n(3)
       real(dbl) :: rsfe
      end type
!
      type sfera
       real(dbl) :: x
       real(dbl) :: y
       real(dbl) :: z
       real(dbl) :: r
      end type
!
      end module
