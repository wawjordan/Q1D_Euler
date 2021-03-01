!=========================== calc_sound_speed ================================80
!>
!!
!<
!=============================================================================80

function calc_sound_speed(gamma,R_gas,Temp)

  use set_precision, only : prec
  implicit none

  real(prec), intent(in) :: gamma, R_gas, temperature
  real(prec)             :: speed_of_sound

  speed_of_sound = sqrt(gamma*R_gas*temperature)

end function calc_sound_speed
