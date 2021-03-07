module geometry

  use set_precision,   only : prec
  use set_constants,   only : zero, one, two, half, pi
  use fluid_constants, only : R_gas, gamma
  use set_inputs,      only : area
  use soln_type
  use grid_type

  implicit none
  
  private
  
  public : setup_geometry, teardown_geometry
  
  


  contains
  
  subroutine setup_geometry
    
    implicit none
    
    type( soln_t ) :: soln
    type( grid_t ) :: grid
    
    call allocate_grid( grid )
    call allocate_soln( soln )
    
    grid%Ai = area( grid%xi )
    grid%Ac = area( grid%xc )
    
  end subroutine setup_geometry
  
  subroutine teardown_geometry
    
    implicit none
    
    type( soln_t ) :: soln
    type( grid_t ) :: grid
    
    call deallocate_grid( grid )
    call deallocate_soln( soln )
    
  end subroutine teardown_geometry

end module geometry
