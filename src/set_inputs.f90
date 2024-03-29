module set_inputs

  use set_precision,   only : prec
  use set_constants,   only : zero, one, pi
  use fluid_constants, only : R_gas, gamma

  implicit none

  private

  public :: imax, neq, xmin, xmax, n_ghost_cells
  public :: i_high, i_low, ig_high, ig_low
  public :: Astar, area, darea
  public :: CFL, k2, k4, eps, tol
  public :: max_iter, max_newton_iter, newton_tol
  public :: iSS, shock, ramp
  public :: p0, T0, a0, rho0, pb, p_ratio
  public :: set_derived_inputs, read_in
  
  integer :: imax    = 128
  integer :: i_low   = 10
  integer :: i_high  = 10
  integer :: ig_low  = 10
  integer :: ig_high = 10
  integer :: neq  = 3
  integer :: iSS  = 1
  integer :: shock = 0
  integer :: ramp = 0
  integer :: max_iter = 150000
  integer :: n_ghost_cells   = 2
  integer :: max_newton_iter = 1000

  real(prec) :: newton_tol = 1.0e-15_prec
  real(prec) :: tol = 1.0e-11_prec
  real(prec) :: eps        = 1.0e-3_prec
  real(prec) :: p0         = 300.0_prec
  real(prec) :: T0         = 600.0_prec
  real(prec) :: Astar      = 0.2_prec
  real(prec) :: a0         = zero
  real(prec) :: rho0       = zero
  real(prec) :: xmin       = -one
  real(prec) :: xmax       = one
  real(prec) :: CFL        = 0.1_prec
  real(prec) :: p_ratio    = 0.5_prec
  real(prec) :: pb         = 150000_prec
  real(prec) :: k2         = 1.0_prec/2.0_prec
  real(prec) :: k4         = 1.0_prec/64.0_prec

  contains


  subroutine read_in
      implicit none
      character(len=10) :: discard              ! ignore character string
      open(25,file='input.dat',status='old')
      read(25,*) discard
      read(25,*) discard
      read(25,*) discard, imax
      read(25,*) discard, shock
      read(25,*) discard, ramp
      read(25,*) discard, p_ratio
      read(25,*) discard, CFL
      read(25,*) discard, k2
      read(25,*) discard, k4
      close(25)
      pb = p_ratio*p0*1000_prec
  end subroutine read_in


  !=================================== area ==================================80
  !>
  !! Description: Calculates area distribution for nozzle.
  !!
  !! Inputs:      x:   Position coordinate along x-axis.
  !!
  !! Outputs:     area: Area at specified coordinate.
  !<
  !===========================================================================80
  function area(x)

    real(prec) :: area
    real(prec), intent(in) :: x

    area = 0.2_prec + 0.4_prec*( one + sin(pi*(x-0.5_prec)))

  end function area
  

  !=================================== darea =================================80
  !>
  !! Description: Calculates area derivative distribution for nozzle.
  !!
  !! Inputs:      x:   Position coordinate along x-axis.
  !!
  !! Outputs:     darea: analytic dA/dx evaluated at specified coordinate.
  !<
  !===========================================================================80
  function darea(x)

    real(prec) :: darea
    real(prec), intent(in) :: x

    darea = 0.4_prec*pi*cos(pi*(x-0.5_prec))

  end function darea

  !=========================== set_derived_inputs ============================80
  !>
  !! Description: Sets derived quantities and prints to STDOUT.
  !!
  !<
  !===========================================================================80
  subroutine set_derived_inputs

    a0   = sqrt(gamma*R_gas*T0)
    rho0 = 1000.0_prec*p0/(R_gas*T0)
    i_low = 1
    i_high = imax
    ig_low  = 1 - n_ghost_cells
    ig_high = imax + n_ghost_cells

    write(*,'(A8,F20.14,A13)') 'R     = ', R_gas, ' [J/(kmol*K)]'
    write(*,'(A8,F20.14)')     'gamma = ', gamma
    write(*,'(A8,F20.14,A6)')  'a_0   = ', a0, ' [m/s]'
    write(*,'(A8,F20.14,A9)')  'rho_0 = ', rho0, ' [kg/m^3]'
    write(*,'(A8,F20.14,A6)')  'P_0   = ', p0, ' [kPa]'
    write(*,'(A8,F20.14,A4)')  'T_0   = ', T0, ' [K]'
    write(*,'(A8,F20.14,A6)')  'A*    = ', Astar, ' [m^2]'

  end subroutine set_derived_inputs

end module set_inputs
