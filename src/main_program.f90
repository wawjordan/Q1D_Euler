program main_program
  
 ! use set_precision, only : prec
  use set_constants, only : set_derived_constants
  use fluid_constants, only : set_fluid_constants
  use set_inputs, only : set_derived_inputs
  use variable_conversion
  use geometry, only : setup_geometry, teardown_geometry
  use init_problem, only : initialize
  use grid_type
  use soln_type
  use exact_q1d_type
  implicit none
  
  character(len=100) :: header_str
  integer :: i
  type( grid_t )      :: grid
  type( soln_t )      :: soln
  type( exact_q1d_t ) :: ex_soln
  
  call set_derived_constants
  call set_fluid_constants
  call set_derived_inputs
  call setup_geometry(grid,soln)
  
  call initialize(grid,soln)
  call allocate_exact_q1d( ex_soln, grid )
  
  call solve_exact_q1d(ex_soln, grid)
  
  write(*,*) 'Exact solution at cell interfaces:'
  write(header_str,*) '|    x   |    A    |         M         |'// &
  &  '        rho        |         u        |         p        |'
  100 format(2(F9.4),4(F20.14))
  write(*,*) trim(adjustl(header_str))
  do i = 1,size(grid%xi)
    write(*,100) grid%xi(i), grid%Ai(i), ex_soln%Mi(i), &
                 ex_soln%Vi(i,1), ex_soln%Vi(i,2), ex_soln%Vi(i,3)/1000.0_prec
  end do
  
  write(*,*)
  
  write(*,*) 'Exact solution at cell centers:'
  write(header_str,*) '|    x   |    A    |         M         |'// &
  &  '        rho        |         u        |         p        |'
  write(*,*) trim(adjustl(header_str))
  do i = 1,size(grid%xc)
    write(*,100) grid%xc(i), grid%Ac(i), ex_soln%Mc(i), &
                 ex_soln%Vc(i,1), ex_soln%Vc(i,2), ex_soln%Vc(i,3)/1000.0_prec
  end do
  
  write(*,*)
  
  write(*,*) 'Initial solution values at cell centers:'
  write(header_str,*) '|    x   |    A    |         M         |'// &
  &  '        rho        |         u        |         p        |'
  write(*,*) trim(adjustl(header_str))
  do i = 1,size(grid%xc)
    write(*,100) grid%xc(i), grid%Ac(i), soln%M(i), soln%V(i,1), soln%V(i,2), soln%V(i,3)/1000.0_prec
  end do
  
  call deallocate_exact_q1d( ex_soln )
  call teardown_geometry(grid,soln)

end program main_program
