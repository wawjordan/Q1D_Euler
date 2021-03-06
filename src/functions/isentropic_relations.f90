! Calculates primitive variable vector for given Mach number and stagnation
! conditions
subroutine isentropic_relations(M,V,T)

  use set_precision,   only : prec
  use set_constants,   only : one, half
  use fluid_constants, only : gamma, R_gas
  use set_inputs,      only : T0, p0
!  use speed_of_sound
  
  implicit none
  
  real(prec), external :: speed_of_sound
  real(prec), dimension(:,:), intent(inout) :: V
  real(prec), dimension(:),   intent(inout) :: M
  real(prec), dimension(:),   intent(inout) :: T

  
  real(prec), dimension(:), allocatable :: psi
  real(prec), dimension(:), allocatable :: p
  real(prec), dimension(:), allocatable :: rho
  real(prec), dimension(:), allocatable :: a
  real(prec), dimension(:), allocatable :: u
  
  allocate(psi(size(M,1)), &
           p(size(M,1)),   &
           rho(size(M,1)), &
           a(size(M,1)),   &
           u(size(M,1))    )
  
  psi = one + half*(gamma - one)*M**2
  T   = T0/psi
  p   = p0/psi**(gamma/(gamma - one))
  rho = p/(R_gas*T)
  a   = speed_of_sound(p,rho)
  u   = M*a

!  T = T0/(one + half*(gamma - one)*M**2)
!  V(:,3) = p0/(one + half*(gamma - one)*M**2)**(gamma/(gamma-1))
!  V(:,1) = V(:,3)/(R_gas*T)
!  V(:,2) = M*speed_of_sound(V(:,3),V(:,1))

  V(:,1) = rho
  V(:,2) = u
  V(:,3) = p
  
  deallocate(psi, p, rho, a, u )

end subroutine isentropic_relations
