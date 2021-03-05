!============================= speed_of_sound ================================80
!>
!!
!<
!=============================================================================80
module speed_of_sound_mod

  contains

  function speed_of_sound(pressure,rho)

    use set_precision, only   : prec
    use fluid_constants, only : gamma

    implicit none

    real(prec), intent(in) :: pressure, rho
    real(prec)             :: speed_of_sound

    speed_of_sound = sqrt(gamma*pressure/rho)

  end function speed_of_sound

end module speed_of_sound_mod
