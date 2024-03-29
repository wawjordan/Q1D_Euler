program main_program
  
  use set_precision, only : prec
  use set_constants, only : set_derived_constants
  use fluid_constants, only : set_fluid_constants
  use set_inputs, only : set_derived_inputs, read_in
  use set_inputs, only : max_iter, neq, tol!, shock
  use variable_conversion
  use time_integration
  use basic_boundaries, only : enforce_bndry
  use fluxes
  use other_subroutines
  use geometry, only : setup_geometry, teardown_geometry
  use init_problem, only : initialize
  use grid_type
  use soln_type
  use exact_q1d_type
  implicit none
  
  character(len=100) :: header_str1
  integer :: j, pnorm
  !real(prec), dimension(3) :: rnorm
  !real(prec), dimension(3) :: rinit
  !real(prec), dimension(3) :: DEnorm
  type( grid_t )      :: grid
  type( soln_t )      :: soln
  type( exact_q1d_t ) :: ex_soln  
  
  pnorm = 2
  
  open(50,file='temp.txt',status='unknown')
  
  100 format(1(I0.8),3(G20.7))
  
  write(header_str1,*) " Iter  |   ||density||    |"// &
  & "   ||velocity||    |   ||pressure||    |"
  call set_derived_constants
  call set_fluid_constants
  call read_in
  call output_file_headers
  call set_derived_inputs
  call setup_geometry(grid,soln)
  
  call initialize(grid,soln)
  
  call allocate_exact_q1d( ex_soln )
  call solve_exact_q1d( ex_soln, grid)
  
  if (shock.eq.0) then
    call calc_de( soln, ex_soln, soln%DE, soln%DEnorm, pnorm )
  end if

  call output_soln(grid,soln,ex_soln,0)
  
  call enforce_bndry( soln )
  
  call cons2prim( soln%U, soln%V )
 
  call speed_of_sound(soln%V(:,3),soln%V(:,1),soln%asnd)
  
  call calculate_sources(soln%V,grid%dAc,soln%src)
  
  call calc_time_step(grid%dx,soln%V,soln%asnd,soln%lambda,soln%dt)
  
  call central_flux( soln%U, soln%F )
  
  call prim2cons(soln%U,soln%V)
  
  call jst_damping(soln%lambda,soln%U,soln%V,soln%D) 
  
  soln%F = soln%F + soln%D
  
  call explicit_euler(grid,soln%src,soln%dt,soln%F,soln%U,soln%R)
  
  call residual_norms(soln%R,soln%rinit,pnorm,(/one,one,one/))
   
  write(*,*) 'Residual Norms: Iteration 0'
  write(*,*) header_str1
  write(*,100) j, soln%rinit(1), soln%rinit(2), soln%rinit(3)
  write(*,*)
  write(*,*) 'Relative Residual Norms:'
  
  do j = 1,max_iter
    
    call calculate_sources(soln%V(:,3),grid%dAc,soln%src)
    
    call speed_of_sound(soln%V(:,3),soln%V(:,1),soln%asnd)
    
    call calc_time_step(grid%dx,soln%V,soln%asnd,soln%lambda,soln%dt)
    
    call enforce_bndry( soln )
    
    call cons2prim( soln%U, soln%V )
    
    call central_flux(soln%U, soln%F)
    
    call prim2cons(soln%U,soln%V)
    
    call jst_damping(soln%lambda,soln%U,soln%V,soln%D)
    
    soln%F = soln%F + soln%D
    
    call explicit_euler(grid,soln%src,soln%dt,soln%F,soln%U,soln%R)
    
    call update_mach(soln%V,soln%mach )
    
    call cons2prim(soln%U,soln%V)
    
    if (mod(j,10000)==0) then
      if (shock.eq.0) then
        call calc_de( soln, ex_soln, soln%DE, soln%DEnorm, pnorm )
      end if
      call output_soln(grid,soln,ex_soln,j)
    end if
    
    call residual_norms(soln%R,soln%rnorm,pnorm,soln%rinit)
    
    if (all(soln%rnorm<tol) ) then
      exit
    end if
    if (mod(j,10000)==0) then
      write(*,*) header_str1
    end if
    if (mod(j,1000)==0) then
      if (shock.eq.0) then
        call calc_de( soln, ex_soln, soln%DE, soln%DEnorm, pnorm )
      end if
      write(*,100) j, soln%rnorm(1), soln%rnorm(2), soln%rnorm(3)
    end if
    if (mod(j,10)==0) then
      if (shock.eq.0) then
        call calc_de( soln, ex_soln, soln%DE, soln%DEnorm, pnorm )
      end if
      call output_res(soln,j)
    end if
  end do
  
  write(50,*) imax, j, soln%rnorm(1), soln%rnorm(2), soln%rnorm(3)
  if (shock.eq.0) then 
    call calc_de( soln, ex_soln, soln%DE, soln%DEnorm, 1 )
    write(50,*) imax, 1,soln%DEnorm
    call calc_de( soln, ex_soln, soln%DE, soln%DEnorm, 2 )
    write(50,*) imax, 2,soln%DEnorm
    call calc_de( soln, ex_soln, soln%DE, soln%DEnorm, 0 )
    write(50,*) imax, 0,soln%DEnorm
    call calc_de( soln, ex_soln, soln%DE, soln%DEnorm, pnorm )
  end if
  

  call output_soln(grid,soln,ex_soln,j+1)
  call output_res(soln,j)
  
  call deallocate_exact_q1d( ex_soln )
  call teardown_geometry(grid,soln)
  close(50)
end program main_program
