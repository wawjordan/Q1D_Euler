module other_subroutines
  
  use set_precision, only : prec
  use set_constants, only : zero, one, two, three, half
  use variable_conversion
  use fluxes
  use soln_type
  use grid_type
  
  implicit none
  
contains

subroutine calculate_sources(V,dA,S)
  
  real(prec), dimension(:,:), intent(in) :: V
  real(prec), dimension(:), intent(in) :: dA
  real(prec), dimension(lbound(V,1):ubound(V,1)), intent(out) :: S
  
  S(:) = V(:,3)*dA(:)
  
end subroutine calculate_sources

subroutine jst_damping(lambda,U,V,d)
  
  use set_inputs, only : neq, k2, k4
  real(prec), dimension(:),   intent(in) :: lambda
  real(prec), dimension(:,:), intent(in) :: U,V
  real(prec), dimension(:,:), intent(out) :: d
  real(prec), dimension(lbound(U,1)+1:ubound(U,1)-2,neq) :: D1
  real(prec), dimension(lbound(U,1)+1:ubound(U,1)-2,neq) :: D3
  real(prec), dimension(lbound(U,1)+1:ubound(U,1)-2) :: e2
  real(prec), dimension(lbound(U,1)+1:ubound(U,1)-2) :: e4
  real(prec), dimension(lbound(U,1):ubound(U,1)) :: nu
  real(prec), dimension(lbound(U,1):ubound(U,1)) :: P
  
  integer :: i, i_low, i_high
  
  i_low  = lbound(U,1)
  i_high = ubound(U,1)
  
  P(:) = V(:,3)
  
  do i = i_low+1,i_high-1
    nu(i) = abs(P(i+1)-two*P(i)+P(i-1))/abs(P(i+1)+two*P(i)+P(i-1))
  end do
  nu(i_low)  = nu(i_low+1)
  nu(i_high) = nu(i_high-1)
  
  do i = i_low+1,i_high-2
    e2(i) = k2*max(nu(i-1),nu(i),nu(i+1),nu(i+2))
    e4(i) = max(zero,k4-e2(i))
  end do

  
  do i = i_low+1,i_high-2
    D1(i,:) = lambda(i)*e2(i)*(U(i+1,:)-U(i,:))
    D3(i,:) = lambda(i)*e4(i)*(U(i+2,:)-three*U(i+1,:)+three*U(i,:)-U(i-1,:))
  end do
  
  d(:,:) = D3(:,:) - D1(:,:)
  
end subroutine jst_damping

end module other_subroutines
