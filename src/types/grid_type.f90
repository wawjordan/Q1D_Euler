module grid_type
  
  use set_precision, only : prec
  
  implicit none
  
  type grid_t

    integer :: imax

    real(prec) :: x_min
    real(prec) :: x_max
    
    real(prec), allocatable, dimension(:) :: xi
    real(prec), allocatable, dimension(:) :: xc
    real(prec), allocatable, dimension(:) :: Ai
    real(prec), allocatable, dimension(:) :: Ac
    real(prec), allocatable, dimension(:) :: dAc
    
  end type grid_t
  
  subroutine allocate_grid( grid )
    
    use set_constants, only : zero, one
    use set_inputs   , only : imax, x_min, x_max
    
    implicit none
    
    type( grid_t ), intent( inout ) :: grid
    integer :: i, i_low, i_high
    
    i_low  = 1
    i_high = imax
    
    grid%dx = (x_max - x_min)/float(imax)

    allocate( grid%xi(i_low:i_high+1), &
              grid%xc(i_low:i_high)  , &
              grid%Ai(i_low:i_high+1), &
              grid%Ac(i_low:i_high)  , &
              grid%dAc(i_low:i_high)   )
    
    grid%xi = [ (x_min + float(i-1)/float(imax-1)*(x_max-x_min),i=1,i=imax) ]
    grid%xc = [ (half*(grid%xi(i) + grid%xi(i+1)),i=1,i=imax )
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

