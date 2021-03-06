module exact_q1d_type

  use set_precision, only : prec
  use set_constants, only : half, one, two
  use set_inputs,    only : g, gp1, gm1, Aq, Astar, iSS, eps

  implicit none

!============================== exact_q1d_t ==================================80
!>
!! Description: Derived type for analytic solution for quasi-1D nozzle flow.
!<
!=============================================================================80
  type exact_q1d_t

    real(prec), allocatable, dimension(:) :: M
    real(prec), allocatable, dimension(:) :: T
    real(prec), allocatable, dimension(:) :: rho
    real(prec), allocatable, dimension(:) :: u
    real(prec), allocatable, dimension(:) :: p

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
    use set_constants, only : zero, one
    use set_inputs, only : iq
    implicit none
    type(exact_q1d_t), intent(inout) :: soln
    integer :: i1, i2

    i1 = 1
    i2 = iq

    allocate( soln%M(i1:i2), soln%T(i1:i2), soln%rho(i1:i2), &
              soln%u(i1:i2), soln%p(i1:i2) )
    soln%M = zero
    soln%T = zero
    soln%rho = zero
    soln%u = zero
    soln%p = zero
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
    deallocate( soln%M, soln%T, soln%rho, soln%u, soln%p )
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
  subroutine calc_variables( soln )
    use set_precision, only : prec
    use set_constants, only : half
    use set_inputs, only : g, gm1, R, p0, T0
    implicit none
    type(exact_q1d_t), intent(inout) :: soln

    soln%T = T0/( 1 + half*gm1*soln%M**2 )
    soln%p = p0/( 1 + half*gm1*soln%M**2 )**(g/gm1)
    soln%rho = soln%p/(R*soln%T)
    soln%u = soln%M*sqrt(g*R*soln%T)

  end subroutine calc_variables

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
    use set_precision, only : prec
    use set_constants, only : zero
    use set_inputs, only : iq, Astar, Aq, iSS, eps, max_newton_iter
    use subroutines
    implicit none

    type(exact_q1d_t), intent(inout) :: soln
    real(prec) :: x0
    real(prec) :: x1
    real(prec), dimension(max_newton_iter+1) :: xk
    real(prec), dimension(max_newton_iter) :: e
    integer :: i

    xk = -9999.99_prec
    e = -9999.99_prec

    x0 = eps
    x1 = one
    ! call newton_safe( f, df, x0, x1, soln%M(1), xk, e)
    call newton_safe2( Aq(1), f1, df1, x0, x1, soln%M(1), xk, e)
    do i = 2,iq
      if ( (iSS==1).and.(Aq(i) > Aq(i-1)) ) then
        x0 = one - eps
        x1 = 10.0_prec
      else
        x0 = eps
        x1 = one+eps
      endif
      if (Aq(i)==Astar) then
        soln%M(i) = one
      else
        call newton_safe2( Aq(i), f1, df1, x0, x1, soln%M(i), xk, e)
      endif

      call calc_variables( soln )
      ! write(*,'(F20.14)') xk
      ! write(*,*) '_____________________________________________________________________'
    end do
  end subroutine solve_exact_q1d

  !=================================== f =====================================80
  !>
  !! Description: Fixed-point form of Area-Mach number relation.
  !!
  !! Inputs:      M:    Mach number.
  !!
  !! Outputs:     f:    f(M).
  !<
  !===========================================================================80
  function f (M)
    real(prec) :: f
    real(prec), intent (in) :: M
    f = ((two/gp1)*(one+half*gm1*M**2))**(gp1/gm1) - ((one/Astar)**2)*M**2
    return
  end function f

  !================================== df =====================================80
  !>
  !! Description: Derivative of fixed-point form of Area-Mach number relation.
  !!
  !! Inputs:      M:    Mach number.
  !!
  !! Outputs:     df:    df(M)/dM.
  !<
  !===========================================================================80
  function df (M)
    real(prec) :: df
    real(prec), intent (in) :: M
    df = two*M*( ((two/gp1)*(one+half*gm1*M**2))**(two/gm1) - (one/Astar)**2 )
    return
  end function df

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
    f1 = ((two/gp1)*(one+half*gm1*M**2))**(gp1/gm1) - ((A1/Astar)**2)*M**2
    return
  end function f1

  !================================== df =====================================80
  !>
  !! Description: Derivative of fixed-point form of Area-Mach number relation
  !!              as a function of Mach number and area.
  !!
  !! Inputs:      M:    Mach number.
  !!              A1:   Local area.
  !!
  !! Outputs:     df:    df(M)/dM.
  !<
  !===========================================================================80
  function df1 (M,A1)
    real(prec) :: df1
    real(prec), intent (in) :: M
    real(prec), intent (in) :: A1
    df1 = two*M*( ((two/gp1)*(one+half*gm1*M**2))**(two/gm1) - (A1/Astar)**2 )
    return
  end function df1

end module exact_q1d_type
