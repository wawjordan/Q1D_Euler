program main_program
  
 ! use set_precision, only : prec
  use set_constants, only : set_derived_constants
  use fluid_constants, only : set_fluid_constants
  use set_inputs, only : set_derived_inputs
  use geometry, only : setup_geometry, teardown_geometry
  use grid_type
  use soln_type
  implicit none
  
  !character(len=100) :: header_str
  integer :: i
  type( grid_t ) :: grid
  type( soln_t ) :: soln
  
  call set_derived_constants
  call set_fluid_constants
  call set_derived_inputs
  call setup_geometry(grid,soln)
  
  do i = 1, size(grid%xi)
    write(*,*) grid%xi(i), grid%Ai(i)
  enddo
  write(*,*)
  do i = 1, size(grid%xc)
    write(*,*) grid%xc(i), grid%Ac(i)
  enddo
  !call allocate_exact_q1d( ex_soln )
  
  !call solve_exact_q1d(ex_soln)
  !call isentropic_relations(ex_soln%M(1), ex_soln%V(1,1), ex_soln%T(1))
  !write(header_str,*) '|    x   |    A    |         M         |'// &
  !&  '        rho        |         u        |         p        |'
  !100 format(2(F9.4),4(F20.14))
  !write(*,*) trim(adjustl(header_str))
  !do i = 1,imax
  !  write(*,100) xq(i), Aq(i), ex_soln%M(i), ex_soln%V(i,1), ex_soln%V(i,2), ex_soln%V(i,3)
  !end do
  !write(*,*) 'T = ', ex_soln%T
  !call deallocate_exact_q1d( ex_soln )
  
  call teardown_geometry(grid,soln)

end program main_program
