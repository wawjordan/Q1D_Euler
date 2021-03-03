!============================= speed_of_sound ================================80
!>
!!
!<
!=============================================================================80
module speed_of_sound_mod
contains

function speed_of_sound(gamma,R_gas,Temp)

  use set_precision, only : prec
  implicit none

  real(prec), intent(in) :: gamma, R_gas, temperature
  real(prec)             :: speed_of_sound

  speed_of_sound = sqrt(gamma*R_gas*temperature)

end function speed_of_sound

end module speed_of_sound_mod
