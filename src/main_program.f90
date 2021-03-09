program main_program
  
 ! use set_precision, only : prec
  use set_constants, only : set_derived_constants
  use fluid_constants, only : set_fluid_constants
  use set_inputs, only : set_derived_inputs, imax, i_low, i_high
  use variable_conversion
  use time_integration
  use basic_boundaries
  use fluxes
  use other_subroutines
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
  
  write(*,*) 'grid%xc:  ','low = ',lbound(grid%xc,1),'  high= ',ubound(grid%xc,1)
  write(*,*) 'grid%xi:  ','low = ',lbound(grid%xi,1),'  high= ',ubound(grid%xi,1)
  write(*,*) 'soln%U:  ','low = ',lbound(soln%U,1),'  high= ',ubound(soln%U,1)
  write(*,*) 'soln%F:  ','low = ',lbound(soln%F,1),'  high= ',ubound(soln%F,1)
  call solve_exact_q1d( ex_soln, grid)
  
  !call calc_time_step(grid%dx,soln%V,soln%lambda,soln%dt)
  
  !call sub_in_bndry( soln%M, soln%U, soln%V )
  !call sup_out_bndry( soln%U, soln%V )
    
  call central_flux(soln%U, soln%F)
  call cons2prim(soln%U,soln%V)
  call prim2cons(soln%U,soln%V)
  write(*,*) 'soln%U:  ','low = ',lbound(soln%U,1),'  high= ',ubound(soln%U,1)
  call jst_damping(soln%lambda,soln%U,soln%V,soln%d)
   
  write(*,*) 'Exact solution at cell interfaces:'
  write(header_str,*) '|    x   |    A    |         M         |'// &
  &  '        rho        |         u        |         p        |'
  100 format(2(F9.4),4(F20.14))
  write(*,*) trim(adjustl(header_str))
  do i = 0,imax
    write(*,100) grid%xi(i), grid%Ai(i), ex_soln%Mi(i), &
                 ex_soln%Vi(i,1), ex_soln%Vi(i,2), ex_soln%Vi(i,3)/1000.0_prec
  end do
  
  write(*,*)
  
  write(*,*) 'Exact solution at cell centers:'
  write(header_str,*) '|    x   |    A    |         M         |'// &
  &  '        rho        |         u        |         p        |'
  write(*,*) trim(adjustl(header_str))
  do i = 1,imax
    write(*,100) grid%xc(i), grid%Ac(i), ex_soln%Mc(i), &
                 ex_soln%Vc(i,1), ex_soln%Vc(i,2), ex_soln%Vc(i,3)/1000.0_prec
  end do
  
  write(*,*)
  
!  write(*,*) 'Initial solution values at cell centers:'
!  write(header_str,*) '|    x   |    A    |         M         |'// &
!  &  '        rho        |         u        |         p        |'
!  write(*,*) trim(adjustl(header_str))
!  do i = i_low,i_high
!    write(*,100) grid%xc(i), grid%Ac(i), soln%M(i), soln%V(i,1), soln%V(i,2), soln%V(i,3)/1000.0_prec
!  end do
  
!  write(*,*)
  
!  write(*,*) 'Conserved values at cell centers:'
!  write(header_str,*) '|    x   |    A    |         M         |'// &
!  &  '        U(1)        |        U(2)       |        U(3)       |'
!  write(*,*) trim(adjustl(header_str))
!  do i = i_low,i_high
!    write(*,100) grid%xc(i), grid%Ac(i), soln%M(i), soln%U(i,1), soln%U(i,2), soln%U(i,3)
!  end do

!  write(*,*)
!  do i = i_low,i_high
!    write(*,*) 'x  ', grid%xc(i), 'dt  ', soln%dt(i)
!  end do
  
  call deallocate_exact_q1d( ex_soln )
  call teardown_geometry(grid,soln)

end program main_program
