module set_inputs

  use set_precision,   only : prec
  use set_constants,   only : zero, one, two, half, pi
  use fluid_constants, only : R_gas, gamma

  implicit none

  private

  public :: set_derived_inputs, xmin, xmax
  public :: imax, iSS, max_newton_iter, newton_tol, eps
  public :: p0, T0, Astar, a0, rho0, Aq, xq, neq

  integer :: max_newton_iter = 20
  integer :: imax = 5
  integer :: iSS = 1

  real(prec) :: newton_tol = 1.0e-15_prec
  real(prec) :: eps        = 1.0e-3_prec
  real(prec) :: p0         = 300.0_prec
  real(prec) :: T0         = 600.0_prec
  real(prec) :: Astar      = 0.2_prec
  real(prec) :: a0         = zero
  real(prec) :: rho0       = zero
  real(prec) :: xmin       = -one
  real(prec) :: xmax       = one
  integer :: neq        = 3
  real(prec), dimension(:), allocatable :: xq
  real(prec), dimension(:), allocatable :: Aq

  contains


  !=================================== area ==================================80
  !>
  !! Description: Calculates area distribution for nozzle.
  !!
  !! Inputs:      xq:   Position coordinate along x-axis.
  !!
  !! Outputs:     area: Area at specified coordinate.
  !<
  !===========================================================================80
  function area(xq)

    real(prec) :: area
    real(prec), intent(in) :: xq

    area = 0.2_prec + 0.4_prec*( one + sin(pi*(xq-0.5_prec)))

  end function area


  !=========================== set_derived_inputs ============================80
  !>
  !! Description: Sets derived quantities and prints to STDOUT.
  !!
  !<
  !===========================================================================80
  subroutine set_derived_inputs

    implicit none

    integer :: i, i1, i2
    ! real(prec), external :: area
    i  = 1
    i1 = 1
    i2 = imax
    allocate(xq(i1:imax), Aq(i1:imax))

    do i = i1,imax
      xq(i) = xmin + float(i-1)*(xmax-xmin)/float(imax-1)
    end do
    do i = i1,imax
      Aq(i) = area(xq(i))
    end do
    a0   = sqrt(gamma*R_gas*T0)
    rho0 = 1000.0_prec*p0/(R_gas*T0)

    write(*,'(A8,F20.14,A13)') 'R     = ', R_gas, ' [J/(kmol*K)]'
    write(*,'(A8,F20.14)')     'gamma = ', gamma
    write(*,'(A8,F20.14,A6)')  'a_0   = ', a0, ' [m/s]'
    write(*,'(A8,F20.14,A9)')  'rho_0 = ', rho0, ' [kg/m^3]'
    write(*,'(A8,F20.14,A6)')  'P_0   = ', p0, ' [kPa]'
    write(*,'(A8,F20.14,A4)')  'T_0   = ', T0, ' [K]'
    write(*,'(A8,F20.14,A6)')  'A*    = ', Astar, ' [m^2]'

  end subroutine set_derived_inputs

end module set_inputs
