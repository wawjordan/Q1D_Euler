module init_problem
  
  use set_precision,   only : prec
  use set_constants,   only : zero, one, two, half
  use fluid_constants, only : R_gas, gamma
  use set_inputs,      only : p0, T0
  use soln_type
  use grid_type
  
  implicit none
  
  private
  
  !public :: 
  
  contains
  
  subroutine initialize( grid, soln )
    
    implicit none
    
    integer :: i
    
    type( soln_t ), intent(inout) :: soln
    type( grid_t ), intent(inout) :: grid
    
    soln%M = 0.9_prec*grid%xc + one
    
  end subroutine initialize

end module init_problem
