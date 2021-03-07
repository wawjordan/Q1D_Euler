program main_program
  
 ! use set_precision, only : prec
  use set_constants, only : set_derived_constants
  use fluid_constants, only : set_fluid_constants
  use set_inputs
  use soln_type
  use grid_type
  use exact_Q1D_type
  
  implicit none

  type( exact_q1d_t ) :: ex_soln
  !type( soln_t ) :: soln
  !type( grid_t ) :: grid
  character(len=100) :: header_str
  integer :: i
  
  call set_derived_constants
  call set_fluid_constants
  call set_derived_inputs
  !call allocate_grid ( grid )
  !call allocate_soln( soln )
  call allocate_exact_q1d( ex_soln )
  
  call solve_exact_q1d(ex_soln)
  !call isentropic_relations(ex_soln%M(1), ex_soln%V(1,1), ex_soln%T(1))
  write(header_str,*) '|    x   |    A    |         M         |'// &
  &  '        rho        |         u        |         p        |'
  100 format(2(F9.4),4(F20.14))
  write(*,*) trim(adjustl(header_str))
  do i = 1,imax
    write(*,100) xq(i), Aq(i), ex_soln%M(i), ex_soln%V(i,1), ex_soln%V(i,2), ex_soln%V(i,3)
  end do
  write(*,*) 'T = ', ex_soln%T
  call deallocate_exact_q1d( ex_soln )
  !call deallocate_soln( soln )
  !call deallocate_grid( grid )

end program main_program
