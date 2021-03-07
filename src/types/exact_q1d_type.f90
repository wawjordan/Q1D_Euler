module exact_q1d_type

  use set_precision, only : prec
  use set_constants, only : zero, one, two, half
  use fluid_constants, only : R_gas, gamma
  use set_inputs, only : imax, Astar, Aq, iSS, eps, max_newton_iter
  
  implicit none
  
  private
  
  public :: exact_q1d_t
  public :: allocate_exact_q1d, deallocate_exact_q1d
  public :: solve_exact_q1d
  
!============================== exact_q1d_t ==================================80
!>
!! Description: Derived type for analytic solution for quasi-1D nozzle flow.
!<
!=============================================================================80
  type exact_q1d_t

    real(prec), allocatable, dimension(:)   :: M
    real(prec), allocatable, dimension(:)   :: T
    real(prec), allocatable, dimension(:,:) :: V

  end type exact_q1d_t


contains

  !======================== allocate_exact_q1d ===============================80
  !>
  !! Description: Allocates exact_q1d_t type.
  !!
  !! Inputs:      soln:    exact_q1d_t type (unallocated).
  !!
  !! Outputs:     soln:    exact_q1d_t type (allocated).
  !<
  !===========================================================================80
  subroutine allocate_exact_q1d( soln )
    use set_inputs, only : imax, neq
    implicit none
    type(exact_q1d_t), intent(inout) :: soln
    integer :: i1, i2

    i1 = 1
    i2 = imax
    
    allocate( soln%M(i1:i2), soln%T(i1:i2), soln%V(i1:i2,neq) )
    soln%M = zero
    soln%T = zero
    soln%V = zero
    
  end subroutine allocate_exact_q1d

  !====================== deallocate_exact_q1d ===============================80
  !>
  !! Description: Deallocates exact_q1d_t type.
  !!
  !! Inputs:      soln:    exact_q1d_t type (allocated).
  !!
  !! Outputs:     soln:    exact_q1d_t type (deallocated).
  !<
  !===========================================================================80
  subroutine deallocate_exact_q1d( soln )
    
    implicit none
    
    type(exact_q1d_t), intent(inout) :: soln
    
    deallocate( soln%M, soln%T, soln%V )
    
  end subroutine deallocate_exact_q1d

  !============================ calc_variables ===============================80
  !>
  !! Description: Calculates primitive variables from isentropic relations.
  !!
  !! Inputs:      soln:    exact_q1d_t type.
  !!
  !! Outputs:     soln:    exact_q1d_t.
  !<
  !===========================================================================80
!  subroutine calc_variables( soln )
!    use set_inputs, only : gamma, 
!    use set_constants, only : half
!
!    soln%T = T0/( 1 + half*gm1*soln%M**2 )
!    soln%p = p0/( 1 + half*gm1*soln%M**2 )**(g/gm1)
!    soln%rho = soln%p/(R*soln%T)
!    soln%u = soln%M*sqrt(g*R*soln%T)
!
!  end subroutine calc_variables

  !=========================== solve_exact_q1d ===============================80
  !>
  !! Description: Calculates Mach number in nozzle from area relation.
  !!
  !! Inputs:      soln:    exact_q1d_t type.
  !!
  !! Outputs:     soln:    exact_q1d_t type.
  !<
  !===========================================================================80
  subroutine solve_exact_q1d(soln)

    use set_inputs, only : imax, Astar, Aq, iSS, eps, max_newton_iter

    !use isentropic_relations 
    implicit none

    type(exact_q1d_t), intent(inout) :: soln
    real(prec) :: x0
    real(prec) :: x1
    real(prec), dimension(2) :: xk
    real(prec) :: e
    integer :: i

    xk = -9999.99_prec
    e = -9999.99_prec
    
    write(*,*) lbound(soln%M), ubound(soln%M)
    write(*,*) lbound(soln%T), ubound(soln%T)
    
    x0 = eps
    x1 = one
    call newton_safe( Aq(1), f1, df1, x0, x1, soln%M(1), xk, e)
    
    write(*,*) lbound(soln%T), ubound(soln%T)
    call isentropic_relations( soln%M, soln%V, soln%T )
    do i = 2,imax
      if ( (iSS==1).and.(Aq(i) > Aq(i-1)) ) then
        x0 = one - eps
        x1 = 10.0_prec
      else
        x0 = eps
        x1 = one+eps
      endif
!      if (Aq(i)==Astar) then
!        soln%M(i) = one
!      else
        call newton_safe( Aq(i), f1, df1, x0, x1, soln%M(i), xk, e)
!      endif

      !call isentropic_relations( soln%M(i), soln%V(i,:), soln%T(i) )
      ! write(*,'(F20.14)') xk
      ! write(*,*) '_____________________________________________________________________'
    end do
    
!    call isentropic_relations( soln%M(i), soln%V(i,:), soln%T(i) )
    
  end subroutine solve_exact_q1d

  !=================================== f1 ====================================80
  !>
  !! Description: Fixed-point form of Area-Mach number relation as a function
  !!              of Mach number and area.
  !!
  !! Inputs:      M:    Mach number.
  !!              A1:   Local area.
  !!
  !! Outputs:     f:    f(M).
  !<
  !===========================================================================80
  function f1 (M,A1)
    real(prec) :: f1
    real(prec), intent (in) :: M
    real(prec), intent (in) :: A1
    f1 = ((two/(gamma+1))*(one+half*(gamma-1)*M**2))**((gamma+1)/(gamma-1)) - ((A1/Astar)**2)*M**2
    return
  end function f1

  !================================== df1 ====================================80
  !>
  !! Description: Derivative of fixed-point form of Area-Mach number relation
  !!              as a function of Mach number and area.
  !!
  !! Inputs:      M:    Mach number.
  !!              A1:   Local area.
  !!
  !! Outputs:     df1:    df(M)/dM.
  !<
  !===========================================================================80
  function df1 (M,A1)
    real(prec) :: df1
    real(prec), intent (in) :: M
    real(prec), intent (in) :: A1
    df1 = two*M*( ((two/(gamma+1))*(one+half*(gamma-1)*M**2))**(two/(gamma-1)) - (A1/Astar)**2 )
    return
  end function df1


  !============================== newton_safe3 ===============================80
  !>
  !! Description: Safe-guarded 1D Newton's method, but with hard-coded input
  !!              parameter 'A1' and reduced memory overhead because I'm dumb
  !!              and hard-coded array sizes in a previous version of the code.
  !!
  !! Inputs:      A1:   Area.
  !!              f:    two variable function f(x,A1).
  !!              df:   Derivative of f(x) w.r.t. x.
  !!              a:    Beginning of bracketing interval.
  !!              b:    End of bracketing interval.
  !!
  !! Outputs:     x:    Zero of f(x).
  !!              xk:   Array of last two iterates.
  !!              e:    approximate relative iterative error.
  !<
  !===========================================================================80
  subroutine newton_safe( A1, f1, df1, a, b, x, xk, e)

    use set_precision, only : prec
    use set_constants, only : half, one, zero
    use set_inputs

    implicit none

    integer :: k = 1
    real(prec), external :: f1, df1
    real(prec), intent(in) :: a, b, A1
    real(prec)             :: a2, b2
    real(prec), intent(out) :: x
    real(prec) :: x_new
    real(prec), intent(out), dimension(2) :: xk
    real(prec), intent(out) :: e

    xk = zero

    a2 = a
    b2 = b
    xk(1) = a2
    xk(2) = xk(1) - f1(xk(1),A1)/df1(xk(1),A1)

    if ( (xk(2) < a2).or.(xk(2) > b2) ) then
      xk(2) = a2 + half*(b2-a2)
      if ( sign(one,f1(a2,A1)) == sign(one,f1(xk(2),A1)) ) then
        a2 = xk(2)
      else
        b2 = xk(2)
      endif
    endif

    e = abs(xk(2)-xk(1))/abs(xk(2))
    do k = 2, max_newton_iter
      if (e < newton_tol) exit
      x_new = xk(2) - f1(xk(2),A1)/df1(xk(2),A1)
      if ( (x_new < a2).or.(x_new > b2) ) then
        x_new = a2 + half*(b2-a2)
        if ( sign(one,f1(a2,A1)) == sign(one,f1(x_new,A1)) ) then
          a2 = x_new
        else
          b2 = x_new
        endif
      endif
      xk(1) = xk(2)
      xk(2) = x_new
      e = abs(xk(2)-xk(1))/abs(xk(2))
    end do

    x = xk(2)

  end subroutine newton_safe

end module exact_q1d_type
