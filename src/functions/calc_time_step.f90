subroutine calc_time_step( lambda, dt,  V )
  
  use set_precision, only : prec
  use fluid_constants, only : gamma
  use set_inputs, only : CFL
  
  implicit none
  
  real(prec), dimension(:), intent(inout)   :: lambda
  real(prec), dimension(:), intent(inout)   :: dt
  real(prec), dimension(:,:), intent(inout) :: V
  real(prec), external :: speed_of_sound
  
  lambda(:) = abs(V(:,2)) + speed_of_sound(V(:,3),V(:,1))
  dt(:)     = CFL*dx/lambda(:)

end subroutine calc_time_step
