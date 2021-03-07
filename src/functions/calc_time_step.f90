subroutine calc_time_step( lambda, dt,  V )
  
  use set_precision,       only : prec
  use fluid_constants,     only : gamma
  use set_inputs,          only : CFL, dx
  use variable_conversion, only : speed_of_sound
  
  implicit none
  
  real(prec), dimension(:),   intent(inout) :: lambda
  real(prec), dimension(:),   intent(inout) :: dt
  real(prec), dimension(:,:), intent(inout) :: V
  real(prec), dimension(size(lambda))       :: sound_speed
  
  call speed_of_sound(V(:,3),V(:,1),sound_speed)
  
  lambda(:) = abs(V(:,2)) + sound_speed
  dt(:)     = CFL*dx/lambda(:)
  
end subroutine calc_time_step
