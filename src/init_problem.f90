module init_problem
  
  use set_precision,   only : prec
  use set_constants,   only : zero, one, two, half
  use fluid_constants, only : R_gas, gamma
  use set_inputs,      only : p0, T0
  use soln_type
  use grid_type
  use variable_conversion
  
  implicit none
  
  private
  
  public :: initialize
  
  contains
  
  subroutine initialize( grid, soln )
    
    implicit none
    
    integer :: i
    
    type( soln_t ), intent(inout) :: soln
    type( grid_t ), intent(inout) :: grid
    
    soln%M(:) = 0.9_prec*grid%xc(:) + one
    call isentropic_relations(soln%M,soln%V,soln%T)
    
  end subroutine initialize

end module init_problem
