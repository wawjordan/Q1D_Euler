module time_integration
  
  use set_precision, only : prec
  use set_constants, only : one, half
  use grid_type
  use soln_type
  use variable_conversion
  use fluxes
  
  implicit none
  
contains

subroutine calc_time_step( lambda, dt,  V )

  use fluid_constants,     only : gamma
  use set_inputs,          only : CFL, dx

  real(prec), dimension(:),   intent(inout) :: lambda
  real(prec), dimension(:),   intent(inout) :: dt
  real(prec), dimension(:,:), intent(inout) :: V
  real(prec), dimension(size(lambda))       :: sound_speed

  call speed_of_sound(V(:,3),V(:,1),sound_speed)

  lambda(:) = abs(V(:,2)) + sound_speed
  dt(:)     = CFL*dx/lambda(:)

end subroutine calc_time_step

 
subroutine explicit_euler(grid,S,dt,F,U,R)
  
  type(grid_t), intent(in) :: grid
  real(prec), dimension(:,:), intent(inout) :: U
  real(prec), dimension(:,:), intent(in) :: F
  real(prec), dimension(:), intent(in) :: S,dt
  real(prec), dimension(:,:), intent(out) :: R
  integer :: i
  
  i = size(F,1)
  R = F(2:i+1,:)*grid%Ai(2:i+1) - F(1:i,:)*grid%Ai(1:i)
  R(:,2) = R(:,2) - S*grid%x
  
  U = U + dt/(grid%Ac*grid%dx)*R
  
end subroutine explicit_euler

subroutine residual_norms(R,Rnorm,pnorm)
  
  real(prec), dimension(:,:), intent(in) :: R
  real(prec), dimension(1,size(R,2)), intent(out) :: Rnorm
  integer :: pnorm
  
  if (pnorm == 0) then
    Rnorm = abs(maxval(R),1)
  elseif (pnorm == 1) then
    Rnorm = (one/real(size(R,1)))*sum(abs(R),1)
  elseif (pnorm == 2) then 
    Rnorm = sqrt((one/real(size(R,1)))*sum(R**2,1))
  else 
    Rnorm = abs(maxval(R),1)
  end if
  
end subroutine residual_norms

subroutine advance_solution(grid,soln)
  
  type(grid_t) :: grid
  type(soln_t) :: soln
  
  call calculate_time_step(soln%lambda,soln%dt,soln%V)

  call enforce_bndry(soln)

  call extrapolate_variables(soln%U,Uextrap,...)

  call central_flux(Uextrap,Fextrap)

  call jst_damping(soln%lambda,soln%D,soln%E,Fextrap,...)

  call calculate_sources(soln%S,...)

  call explicit_euler(grid,soln%S,soln%dt,Fextrap,soln%U,soln%R)

  call cons2prim(soln%U,soln%V)

  call isentropic_relations(soln%M,soln%V,soln%T) ! calc_mach_number?

  call output_soln()

  call enforce_bndry(soln)

end subroutine advance_solution

end module time_integration
