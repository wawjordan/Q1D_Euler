! converts primitive variable vector to conserved variable vector
function prim2cons(U,V)
  use set_precision, only : prec
  use set_constants, only : one, half
  use set_inputs,    only : gamma

  implicit none

  real(prec), dimension(:,:), intent(inout)  :: U
  real(prec), dimension(:,:), intent(in)     :: V

  U(:,1) = V(:,1)
  U(:,2) = V(:,1)*V(:,2)
  U(:,3) = ( V(:,3)/(gamma - one) ) + half*V(:,1)*V(:,2)**2

end function prim2cons
