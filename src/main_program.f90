program main_program
  
 ! use set_precision, only : prec
  use set_constants, only : set_derived_constants
  use fluid_constants, only : set_fluid_constants
  use set_inputs, only : set_derived_inputs, imax
  use set_inputs, only : i_high, i_low, ig_high, ig_low
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
  
  character(len=100) :: header_str1
  character(len=100) :: header_str2
  character(len=100) :: header_str3
  character(len=100) :: header_str4
  character(len=100) :: header_str5
  integer :: i, j
  type( grid_t )      :: grid
  type( soln_t )      :: soln
  type( exact_q1d_t ) :: ex_soln  

  call set_derived_constants
  call set_fluid_constants
  call set_derived_inputs
  call setup_geometry(grid,soln)
  
  call initialize(grid,soln)
  
  !call allocate_exact_q1d( ex_soln, grid )
  
  !call solve_exact_q1d( ex_soln, grid)
  
!  100 format(2(F9.4),4(G20.7))
!  200 format(1(F9.4),3(G20.7))
!  300 format(2(F9.4),3(G20.7))
!  write(header_str1,*) '|    x   |    A    |         M         |'// &
!  &  '        rho        |         u        |         p        |'
!  write(header_str2,*) '|    x   |        U(1)       |'// &
!  & '       U(2)       |       U(3)       |'
!  write(header_str3,*) '|    x   |        F(1)       |'// &
!  & '       F(2)       |       F(3)       |'
!  write(header_str4,*) '|    x   |    dA   |         S         |'// &
!  & '         dt       |       lambda     |' 
!  write(header_str5,*) '|    x   |        D(1)       |'// &
!  & '       D(2)       |       D(3)       |'

!  write(*,*) 'Initial Conditions: Primitives'
!  write(*,*) trim(adjustl(header_str1))
!  do i = ig_low,ig_high
!    write(*,100) grid%xc(i), grid%Ac(i), soln%mach(i), soln%V(i,1), soln%V(i,2), soln%V(i,3)/1000.0_prec
!  end do
!  write(*,*)
  
  call sub_in_bndry( soln%mach, soln%U, soln%V )

  call sub_out_bndry( soln%U, soln%V )
  
  call cons2prim( soln%U, soln%V )
 
  call speed_of_sound(soln%V(:,3),soln%V(:,1),soln%asnd)
  call calculate_sources(soln%V,grid%dAc,soln%src)
  call calc_time_step(grid%dx,soln%V,soln%asnd,soln%lambda,soln%dt)
!  write(*,*) 'Calculate Sources and Time step:'
!  write(*,*) trim(adjustl(header_str4))
!  do i = ig_low,ig_high
!    if( (i>=i_low).and.(i<=i_high) ) then
!      write(*,300) grid%xc(i),grid%dAc(i), soln%src(i), soln%dt(i), soln%lambda(i)
!    else 
!      write(*,300) grid%xc(i), zero, zero, zero, soln%lambda(i)
!    end if
!  end do
!  write(*,*)
  !call prim2cons(soln%U,soln%V)
  call central_flux( soln%U, soln%F )

    call prim2cons(soln%U,soln%V)

  call jst_damping(soln%lambda,soln%U,soln%V,soln%D) 
    soln%F = soln%F + soln%D
    call explicit_euler(grid,soln%src,soln%dt,soln%F,soln%U,soln%R)

  do j = 1,800
    
    call calculate_sources(soln%V(:,3),grid%dAc,soln%src)
    
    call speed_of_sound(soln%V(:,3),soln%V(:,1),soln%asnd)
    
    call calc_time_step(grid%dx,soln%V,soln%asnd,soln%lambda,soln%dt)
    
    call sub_in_bndry( soln%mach, soln%U, soln%V )
    
    call sub_out_bndry( soln%U, soln%V )
    
    call cons2prim( soln%U, soln%V )
    
    call central_flux(soln%U, soln%F)
    
    call prim2cons(soln%U,soln%V)
    
    call jst_damping(soln%lambda,soln%U,soln%V,soln%D)
    
    soln%F = soln%F + soln%D
    
    call explicit_euler(grid,soln%src,soln%dt,soln%F,soln%U,soln%R)
    
    call update_mach(soln%V,soln%mach )
    
    call cons2prim(soln%U,soln%V)
    
    if (mod(j,4)==0) then
      call output_soln(grid,soln,j)
    end if
  
  end do
  
  !call deallocate_exact_q1d( ex_soln )
  call teardown_geometry(grid,soln)

end program main_program
