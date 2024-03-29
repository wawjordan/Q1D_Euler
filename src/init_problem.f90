module init_problem
  
  use set_precision,   only : prec
  use set_constants,   only : zero, one, two, half
  use fluid_constants, only : R_gas, gamma
  use set_inputs,      only : p0, T0, eps, ig_low, ig_high
  use soln_type
  use grid_type
  use variable_conversion
  
  implicit none
  
  private
  
  public :: initialize
  
  contains
  
  !================================ initialize  ==============================80
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
  subroutine initialize( grid, soln )
    
    implicit none
    
    integer :: i
    type( soln_t ), intent(inout) :: soln
    type( grid_t ), intent(inout) :: grid
    
    soln%mach = 0.9_prec*grid%xc + one
    
    do i = ig_low,ig_high
      if ( soln%mach(i) < eps ) then
        soln%mach(i) = eps
      end if
    end do
    
    call isentropic_relations(soln%mach,soln%V)
    call prim2cons(soln%U,soln%V)
    
    !do i = ig_low,ig_high
    !  write(*,*)">(U)", soln%U(i,:)
    !end do
    !write(*,*) 'grid%xc:  ','low = ',lbound(grid%xc,1),'  high= ',ubound(grid%xc,1)
    !write(*,*) 'grid%xi:  ','low = ',lbound(grid%xi,1),'  high= ',ubound(grid%xi,1)
    !write(*,*) 'soln%U:   ','low = ',lbound(soln%U,1), '  high= ',ubound(soln%U,1)

  end subroutine initialize

end module init_problem
