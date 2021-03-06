module geometry

  use set_precision,   only : prec
  use set_constants,   only : zero, one, two, half, pi
  use fluid_constants, only : R_gas, gamma
  use set_inputs,      only : imax, area, darea, eps
  use soln_type
  use grid_type

  implicit none
  
  private
  
  public :: setup_geometry, teardown_geometry
  
  contains
  
  !============================= setup geometry  =============================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid : 
  !!              soln : 
  !!
  !! Outputs:     grid :
  !!              soln :
  !<
  !===========================================================================80
  subroutine setup_geometry( grid, soln )
    
    type( soln_t ), intent(inout) :: soln
    type( grid_t ), intent(inout) :: grid
    
    integer :: i, i_low, i_high
    
    call allocate_grid( grid )
    call allocate_soln( soln )
    
    i_low  = lbound(grid%xi,1)
    i_high = ubound(grid%xi,1)
    
    grid%Ai(i_low:0) = area( grid%xi(0) )
    
    do i = 1,imax
      grid%Ai(i) = area( grid%xi(i) )
    end do
    
    grid%Ai(imax:i_high) = area( grid%xi(imax) )
    
    i_low  = lbound(grid%xc,1)
    i_high = ubound(grid%xc,1)
    
    grid%Ac(i_low:1) = area( grid%xc(1) )
    do i = 2,imax-1
      grid%Ac(i) = area( grid%xc(i) )
    end do
    grid%Ac(imax:i_high) = area( grid%xc(imax) )
    
    
    grid%dAc(i_low:1) = darea( grid%xc(1) )
    do i = 2,imax-1
      grid%dAc(i) = darea( grid%xc(i) )
    end do
    grid%dAc(imax:i_high) = darea( grid%xc(imax) )
    !grid%dAc(i_low:1) = darea( grid%xc(1) )
    !do i = 1,imax-1
    !  grid%dAc(i) = ( grid%Ac(i+1) - grid%Ac(i) )/grid%dx
    !end do
    !grid%dAc(imax:i_high) = darea( grid%xc(imax) )
    
  end subroutine setup_geometry
  
  !========================== teardown geometry  =============================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid : 
  !!              soln : 
  !!
  !! Outputs:     grid :
  !!              soln :
  !<
  !===========================================================================80
  subroutine teardown_geometry( grid, soln )
    
    type( soln_t ) :: soln
    type( grid_t ) :: grid
    
    call deallocate_grid( grid )
    call deallocate_soln( soln )
    
  end subroutine teardown_geometry

end module geometry
