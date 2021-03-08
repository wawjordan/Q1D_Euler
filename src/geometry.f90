module geometry

  use set_precision,   only : prec
  use set_constants,   only : zero, one, two, half, pi
  use fluid_constants, only : R_gas, gamma
  use set_inputs,      only : area, darea
  use soln_type
  use grid_type

  implicit none
  
  private
  
  public :: setup_geometry, teardown_geometry
  


  contains
  
  subroutine setup_geometry( grid, soln )
    
    integer :: i
    type( soln_t ), intent(inout) :: soln
    type( grid_t ), intent(inout) :: grid
    
    call allocate_grid( grid )
    call allocate_soln( soln )
    
    do i = 1,size(grid%xi)
      grid%Ai(i) = area( grid%xi(i) )
    end do
    
    do i = 1,size(grid%xc)
      grid%Ac(i) = area( grid%xc(i) )
    end do
    
    do i = 1,size(grid%xc)
      grid%dAc(i) = darea( grid%xc(i) )
    end do
    
  end subroutine setup_geometry
  
  subroutine teardown_geometry( grid, soln )
    
    type( soln_t ) :: soln
    type( grid_t ) :: grid
    
    call deallocate_grid( grid )
    call deallocate_soln( soln )
    
  end subroutine teardown_geometry

end module geometry
