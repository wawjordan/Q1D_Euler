! converts conserved variables to primitve variables
subroutine cons2prim(U,V)

  use set_precision,      only : prec
  use set_constants,      only : one, half
  use fluid_constants,    only : gamma

  implicit none

  real(prec), dimension(:,:), intent(in)    :: U
  real(prec), dimension(:,:), intent(inout) :: V

  V(:,1) = U(:,1)
  V(:,2) = U(:,2)/U(:,1)
  V(:,3) = (gamma - one)*U(:,3) - half*(gamma - one)*U(:,2)**2/U(:,1)

end subroutine cons2prim
