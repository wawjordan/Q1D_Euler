module other_subroutines
  
  use set_precision, only : prec
  use set_constants, only : zero, one, half
  use variable_conversion
  use fluxes
  use soln_type
  use grid_type
  
  implicit none
  
contains

subroutine calculate_sources(V,dA,S)
  
  real(prec), dimension(:,:), intent(in) :: V
  real(prec), dimension(:,:), intent(in) :: dA
  real(prec), dimension(size(V,1)), intent(out) :: S
  
  S = V(:,3)*dA
  
end subroutine calculate_sources

subroutine jst_damping(lambda,U,V,d)
  
  use set_inputs, only : k2, k4
  real(prec), dimension(:),   intent(in) :: lambda
  real(prec), dimension(:,:), intent(in) :: U,V
  real(prec), dimension(:,:), intent(out) :: d
  real(prec), dimension(size(d,1),neq) :: D1
  real(prec), dimension(size(d,1),neq) :: D3
  real(prec), dimension(size(d,1),1) :: e2
  real(prec), dimension(size(d,1),1) :: e4
  real(prec), dimension(size(d,1),1) :: nu
  real(prec), dimension(size(V,1),1) :: P
  
  integer :: i, N
  
  N = size(d,1)
  
  P = V(:,3)
  
  do i = 1,N
    nu(i) = abs(P(i+1)-two*P(i)+P(i-1))/abs(P(i+1)+two*P(i)+P(i-1))
  end do
  
  do i = 2,N-1
    e2(i) = k2*max(nu(i-1),nu(i),nu(i+1),nu(i+2))
    e4(i) = max(0,k4-e2(i))
  end do
  
  do i = 1,N
    D1(i,:) = lambda(i)*e2(i)*(U(i+1,:)-U(i,:))
  end do
  
  do i = 2,N-1
    D3(i,:) = lambda(i)*e4(i)*(U(i+2,:)-three*U(i+1,:)+three*U(i,:)-U(i-1,:))
  end do

  d = D3 - D1
  
end subroutine jst_damping

end module other_subroutines
