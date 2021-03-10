module time_integration
  
  use set_precision, only : prec
  use set_constants, only : one, half
  use set_inputs,    only : neq, i_low, i_high, ig_low, ig_high
  use grid_type
  use soln_type
  use variable_conversion
  use fluxes
  
  implicit none
  
contains

subroutine calc_time_step( dx, V, lambda, dt )

  use fluid_constants,     only : gamma
  use set_inputs,          only : CFL

  real(prec), dimension(:,:), intent(in)  :: V
  !type(grid_t),               intent(in)  :: grid
  real(prec), intent(in) :: dx
  real(prec), dimension(:),   intent(out) :: lambda
  real(prec), dimension(:),   intent(out) :: dt
  
  real(prec), dimension(size(lambda))     :: sound_speed

  call speed_of_sound(V(:,3),V(:,1),sound_speed)

  lambda(:) = abs(V(:,2)) + sound_speed
  !dt(:)     = CFL*grid%dx/lambda(:)
  dt(:) = CFL*dx/lambda(:)
end subroutine calc_time_step

 
subroutine explicit_euler(grid,S,dt,F,U,R)
  
  type(grid_t),               intent(in)    :: grid
  real(prec), dimension(ig_low:ig_high,neq), intent(inout) :: U
  real(prec), dimension(i_low-1:i_high,neq), intent(in)    :: F
  real(prec), dimension(ig_low:ig_high)  , intent(in)    :: S,dt
  real(prec), dimension(i_low:i_high,neq), intent(out)   :: R
  integer :: i
  
  R(:,1) = F(i_low:i_high,1)*grid%Ai(i_low:i_high) - F(i_low-1:i_high-1,1)*grid%Ai(i_low-1:i_high-1)
  R(:,2) = F(i_low:i_high,2)*grid%Ai(i_low:i_high) - F(i_low-1:i_high-1,2)*grid%Ai(i_low-1:i_high-1) &
                                                   - S(i_low:i_high)*grid%xc(i_low:i_high)
  R(:,3) = F(i_low:i_high,3)*grid%Ai(i_low:i_high) - F(i_low-1:i_high-1,3)*grid%Ai(i_low-1:i_high-1)
  
  U(i_low:i_high,1) = U(i_low:i_high,1) + minval(dt(i_low:i_high))/(grid%Ac(i_low:i_high)*grid%dx)*R(:,1)
  U(i_low:i_high,2) = U(i_low:i_high,2) + minval(dt(i_low:i_high))/(grid%Ac(i_low:i_high)*grid%dx)*R(:,2)
  U(i_low:i_high,3) = U(i_low:i_high,3) + minval(dt(i_low:i_high))/(grid%Ac(i_low:i_high)*grid%dx)*R(:,3)
  
end subroutine explicit_euler

subroutine residual_norms(R,Rnorm,pnorm)
  
  real(prec), dimension(:,:), intent(in) :: R
  real(prec), dimension(1,size(R,2)), intent(out) :: Rnorm
  integer :: pnorm
  
  if (pnorm == 0) then
    Rnorm(1,:) = abs(maxval(R,1))
  elseif (pnorm == 1) then
    Rnorm(1,:) = (one/real(size(R,1)))*sum(abs(R),1)
  elseif (pnorm == 2) then 
    Rnorm(1,:) = sqrt((one/real(size(R,1)))*sum(R**2,1))
  else 
    Rnorm(1,:) = abs(maxval(R,1))
  end if
  
end subroutine residual_norms

!subroutine advance_solution(grid,soln)
!  
!  type(grid_t) :: grid
!  type(soln_t) :: soln
!  
!  call calculate_time_step(soln%lambda,soln%dt,soln%V)
!
!  call enforce_bndry(soln)
!
!  call extrapolate_variables(soln%U,Uextrap,...)
!
!  call central_flux(Uextrap,Fextrap)
!
!  call jst_damping(soln%lambda,soln%D,soln%E,Fextrap,...)
!
!  call calculate_sources(soln%S,...)
!
!  call explicit_euler(grid,soln%S,soln%dt,Fextrap,soln%U,soln%R)
!
!  call cons2prim(soln%U,soln%V)
!
!  call isentropic_relations(soln%M,soln%V,soln%T) ! calc_mach_number?
!
!  call output_soln()
!
!  call enforce_bndry(soln)
!
!end subroutine advance_solution

end module time_integration
