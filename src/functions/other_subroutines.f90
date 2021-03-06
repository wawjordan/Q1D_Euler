module other_subroutines
  
  use set_precision, only : prec
  use set_constants, only : zero, one, two, three, half
  use set_inputs, only : imax, neq, i_low, i_high, ig_low, ig_high
  use variable_conversion
  use fluxes
  use soln_type
  use grid_type
  
  implicit none
  
  contains
  
  !============================= calculate_sources ===========================80
  !>
  !! Description: 
  !!
  !! Inputs:      V  : 
  !!              dA : 
  !!
  !! Outputs:     S  : 
  !<
  !===========================================================================80
  subroutine calculate_sources(V,dA,S)
    
    real(prec), dimension(:,:), intent(in) :: V
    real(prec), dimension(:),   intent(in) :: dA
    real(prec), dimension(lbound(V,1):ubound(V,1)), intent(out) :: S
    
    S(:) = V(:,3)*dA(:)
    
  end subroutine calculate_sources
  
  
  !============================= jst_damping =================================80
  !>
  !! Description: 
  !!
  !! Inputs:      lambda : 
  !!              U      : 
  !!              V      :
  !!
  !! Outputs:     d      : 
  !<
  !===========================================================================80
  subroutine jst_damping(lambda,U,V,d)
    
    use set_inputs, only : neq, k2, k4
    
    real(prec), dimension(ig_low:ig_high),   intent(in) :: lambda
    real(prec), dimension(ig_low:ig_high,neq), intent(in) :: U,V
    real(prec), dimension(i_low-1:i_high,neq), intent(out) :: d
    real(prec), dimension(i_low-1:i_high) :: lambda_half
    real(prec), dimension(i_low-1:i_high,neq) :: D1
    real(prec), dimension(i_low-1:i_high,neq) :: D3
    real(prec), dimension(i_low-1:i_high) :: e2
    real(prec), dimension(i_low-1:i_high) :: e4
    real(prec), dimension(ig_low:ig_high) :: nu
    real(prec), dimension(ig_low:ig_high) :: P
    
    integer :: i
    
    P(:) = V(:,3)
    lambda_half = half*(lambda(i_low:i_high+1) + lambda(i_low-1:i_high))
    do i = i_low-1,i_high
      nu(i) = abs(P(i+1)-two*P(i)+P(i-1))/abs(P(i+1)+two*P(i)+P(i-1))
    end do
    nu(i_low-2)  = nu(i_low-1)
    nu(i_high+1) = nu(i_high)
    
    do i = i_low-1,i_high
      if (i == i_low-1) then
        e2(i) = k2*max(nu(i),nu(i+1),nu(i+2))
      elseif(i == i_high) then
        e2(i) = k2*max(nu(i-1),nu(i),nu(i+1))
      else
        e2(i) = k2*max(nu(i-1),nu(i),nu(i+1),nu(i+2))
      end if
      e4(i) = max(zero,k4-e2(i))
    end do
    
    do i = i_low-1,i_high
      D1(i,:) = lambda_half(i)*e2(i)*(U(i+1,:)-U(i,:))
      D3(i,:) = lambda_half(i)*e4(i)*(U(i+2,:)-three*U(i+1,:)+three*U(i,:)-U(i-1,:))
    end do
    
    d(:,:) = D3(:,:) - D1(:,:)
    
    d(i_low-1,:) = 2*d(i_low,:) - d(i_low+1,:)
    d(i_high,:) = 2*d(i_high-1,:) - d(i_high-2,:)

  end subroutine jst_damping
  
  !============================= output_soln =================================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid     : 
  !!              soln     : 
  !!              num_iter : 
  !<
  !===========================================================================80
  subroutine output_soln(grid,soln,num_iter)
    
    type( grid_t ), intent(in) :: grid
    type( soln_t ), intent(in) :: soln
    integer,        intent(in) :: num_iter
    
    integer :: i
    
    open(40,file='q1Dnozzle.dat',status='unknown')
    write(40,*) 'TITLE = "Quasi-1D Nozzle Solution"'
    write(40,*)' variables="x(m)""Area(m^2)""rho(kg/m^3)""u(m/s)"&
             & "Press(N/m^2)""Mach""U1""U2""U3"'
    ! Repeat the following each time you want to write out the solution
    write(40,*) 'zone T="',num_iter,'" '
    write(40,*) 'I=',imax
    write(40,*) 'DATAPACKING=POINT'
    write(40,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE &
             & DOUBLE DOUBLE )'
    do i = 1, imax
    write(40,*)grid%xc(i),grid%Ac(i),soln%V(i,1),soln%V(i,2),soln%V(i,3),&
             & soln%M(i),soln%U(i,1),soln%U(i,2),soln%U(i,3)
    enddo
    
  end subroutine output_soln
  
end module other_subroutines
