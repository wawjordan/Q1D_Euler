module grid_type
  
  use set_precision, only : prec
  
  implicit none
  
  type grid_t

    integer :: imax

    real(prec) :: xmin
    real(prec) :: xmax
    real(prec) :: dx
    
    real(prec), allocatable, dimension(:) :: xi
    real(prec), allocatable, dimension(:) :: xc
    real(prec), allocatable, dimension(:) :: Ai
    real(prec), allocatable, dimension(:) :: Ac
    real(prec), allocatable, dimension(:) :: dAc
    
  end type grid_t

contains
  
  subroutine allocate_grid( grid )
    
    use set_constants, only : zero, one, half
    use set_inputs   , only : imax, xmin, xmax, n_ghost_cells
    
    type( grid_t ), intent( inout ) :: grid
    integer :: i, i_low, i_high
    
    i_low  = 1 - n_ghost_cells
    i_high = imax + n_ghost_cells
    
    grid%dx = (xmax - xmin)/real(imax,prec)

    allocate( grid%xi(i_low-1:i_high), &
              grid%xc(i_low:i_high)  , &
              grid%Ai(i_low-1:i_high), &
              grid%Ac(i_low:i_high)  , &
              grid%dAc(i_low:i_high)   )
    
    grid%xi = [ (xmin + real(i,prec)/real(imax,prec)*(xmax-xmin),i=i_low-1,i_high) ]
    grid%xc = [ (half*(grid%xi(i) + grid%xi(i+1)),i=i_low-1,i_high-1) ]
    grid%Ai = one
    grid%Ac = one
    grid%dAc =  zero
    
  end subroutine allocate_grid
  
  subroutine deallocate_grid( grid )
    
    implicit none
    
    type( grid_t ), intent( inout ) :: grid
    
    deallocate( grid%xi, grid%xc, grid%Ai, grid%Ac, grid%dAc )
    
  end subroutine deallocate_grid
  
end module grid_type

