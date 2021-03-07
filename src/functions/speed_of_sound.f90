!============================= speed_of_sound ================================80
!>
!!
!<
!=============================================================================80
function speed_of_sound(pressure,rho)

  use set_precision, only   : prec
  use fluid_constants, only : gamma

  implicit none

  real(prec) :: pressure, rho
  real(prec) :: speed_of_sound
  
  speed_of_sound = sqrt(gamma*pressure/rho)

end function speed_of_sound
