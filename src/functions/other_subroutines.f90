module other_subroutines
  
  use set_precision, only : prec
  use set_constants, only : zero, one, two, three, half
  use set_inputs, only : imax, neq, i_low, i_high, ig_low, ig_high
  use fluid_constants, only : gamma
  use variable_conversion
  use fluxes
  use soln_type
  use grid_type
  
  implicit none
  
  contains
  
  !============================= calculate_sources ===========================80
  !>
  !! Description: 
  !!
  !! Inputs:      V  : 
  !!              dA : 
  !!
  !! Outputs:     S  : 
  !<
  !===========================================================================80
  subroutine calculate_sources(P,dA,S)
    
    !real(prec), dimension(ig_low:ig_high,neq), intent(in) :: V
    real(prec), dimension(i_low:i_high),   intent(in) :: dA
    real(prec), dimension(ig_low:ig_high),   intent(in) :: P
    real(prec), dimension(i_low:i_high),   intent(out) :: S
    integer :: i
    !S(i_low:i_high) = V(i_low:i_high,3)*dA(i_low:i_high)/1000.0_prec
    S(i_low:i_high) = P(i_low:i_high)*dA(i_low:i_high)
    !do i = i_low,i_high
    !write(*,*) i, "dA : ",dA(i), "V: ",V(i,:)," S: ",S(i)
    !end do
  end subroutine calculate_sources
  
  
  !============================= jst_damping =================================80
  !>
  !! Description: 
  !!
  !! Inputs:      lambda : 
  !!              U      : 
  !!              V      :
  !!
  !! Outputs:     d      : 
  !<
  !===========================================================================80
  subroutine jst_damping(lambda,U,V,d)
    
    use set_inputs, only : k2, k4
    
    real(prec), dimension(ig_low:ig_high),     intent(in) :: lambda
    real(prec), dimension(ig_low:ig_high,neq), intent(in) :: U,V
    real(prec), dimension(i_low-1:i_high,neq), intent(out) :: d
    real(prec), dimension(i_low-1:i_high) :: lambda_half
    real(prec), dimension(i_low-1:i_high,neq) :: D1
    real(prec), dimension(i_low-1:i_high,neq) :: D3
    real(prec), dimension(i_low-1:i_high) :: e2
    real(prec), dimension(i_low-1:i_high) :: e4
    real(prec), dimension(ig_low:ig_high) :: nu
    real(prec), dimension(ig_low:ig_high) :: P
    
    integer :: i
    
    P(:) = V(:,3)
    lambda_half = half*(lambda(i_low:i_high+1) + lambda(i_low-1:i_high))
    do i = i_low-1,i_high
      nu(i) = abs(P(i+1)-two*P(i)+P(i-1))/abs(P(i+1)+two*P(i)+P(i-1))
    end do
    nu(i_low-2)  = abs(2*nu(i_low-1) - nu(i_low))
    nu(i_high+1) = abs(2*nu(i_high) - nu(i_high-1))
    
    do i = i_low-1,i_high
      !if (i == i_low-1) then
      !  e2(i) = k2*max(nu(i),nu(i+1),nu(i+2))
      !elseif(i == i_high) then
      !  e2(i) = k2*max(nu(i-1),nu(i),nu(i+1))
      !else
        e2(i) = k2*max(nu(i-1),nu(i),nu(i+1),nu(i+2))
      !end if
      e4(i) = max(zero,k4-e2(i))
    end do
    
    do i = i_low-1,i_high
      D1(i,:) = lambda_half(i)*e2(i)*(U(i+1,:)-U(i,:))
      D3(i,:) = lambda_half(i)*e4(i)* &
              & ( U(i+2,:) - three*U(i+1,:) + three*U(i,:) - U(i-1,:) )
    end do
    
    d(:,:) = D3(:,:) - D1(:,:)
    
    !d(i_low-1,:) = 2*d(i_low,:) - d(i_low+1,:)
    !d(i_high-1,:) = 2*d(i_high-2,:) - d(i_high-3,:)
    !d(i_high,:) = 2*d(i_high-1,:) - d(i_high-2,:)
    

  end subroutine jst_damping
  
  

  !========================== output_file_headers ============================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid     : 
  !!              soln     : 
  !!              num_iter : 
  !<
  !===========================================================================80
subroutine output_file_headers

    use set_inputs, only : imax, CFL, k2, k4, shock, ramp, p0, pb

    character(len=1024) :: dirname
    character(len=1024) :: filename
    character(len=64) ::   shock_str
    character(len=64) ::   ncells_str
    character(len=64) ::   CFL_str
    character(len=64) ::   kappa2_str
    character(len=64) ::   kappa4_str

    ! Set up output directories
    write (ncells_str  , "(A1,I0.3,A1)") "N"  , imax   , "/"
    if (shock.eq.1) then
      write (shock_str, "(A16,I0.3)") "normal-shock-pb-",int(1000*pb/p0)
    else
      write (shock_str, "(A10)") "isentropic"
    end if
    if (ramp.eq.1) then
      write (CFL_str  , "(A4,I0.3,A5)") "CFL-", int(1000*cfl),"-ramp"
    else
      write (CFL_str  , "(A4,I0.3,A8)") "CFL-", int(1000*cfl),"-no-ramp"
    end if
    write (kappa2_str, "(A4,I0.4)") "_K2-"  , int(1000*k2)
    write (kappa4_str, "(A4,I0.4)") "_K4-"  , int(1000*k4)
    write (dirname, *) adjustl(trim(shock_str)),"/"//  &
    &                   adjustl(trim(ncells_str))
    write (filename,*) trim(CFL_str)//     &
    &                   trim(kappa2_str)// &
    &                   trim(kappa4_str)
    write(*,*) trim(adjustl(filename))
    write(*,*) trim(adjustl(dirname))
    !call execute_command_line ('mkdir -p results/' // adjustl(trim(dirname)))

  ! Set up output files (history and solution)
  ! open(30,file='history.dat',status='unknown')
    !open(30,file= 'results/'//trim(adjustl(dirname))//  &
    !&               trim(adjustl(filename))//'_history.dat',status='unknown')
    !open(30,file='q1Dnozzle_hist.dat',status='unknown')
    !write(30,*) 'TITLE = "Quasi-1D Nozzle Iterative Residual History"'
    !write(30,*) 'variables="Iteration""Time(s)""Res1""Res2""Res3"'

    !open(40,file= 'results/'//trim(adjustl(dirname))//  &
    !&               trim(adjustl(filename))//'_field.dat',status='unknown')
    !open(40,file='q1Dnozzle.dat',status='unknown')
    !write(40,*) 'TITLE = "Quasi-1D Nozzle Solution"'
    !if(shock.eq.0) then
    !  write(40,*) 'variables="x(m)""A(m^2)""rho(kg/m^3)""u(m/s)""p(N/m^2)"  &
    !  & "M""U1""U2""U3""rho-exact""u-exact""p-exact""DE-rho""DE-u""DE-p"'
    !elseif(shock.eq.1) then
    !  write(40,*) 'variables="x(m)""A(m^2)""rho(kg/m^3)""u(m/s)""p(N/m^2)"&
    !  &           "M""U1""U2""U3"'
    !else
    !  write(*,*) 'ERROR! shock must equal 0 or 1!!!'
    !  stop
    !endif
  
  end subroutine output_file_headers



  !============================= output_soln =================================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid     : 
  !!              soln     : 
  !!              num_iter : 
  !<
  !===========================================================================80
  subroutine output_soln(grid,soln,num_iter)
    
    type( grid_t ), intent(in) :: grid
    type( soln_t ), intent(in) :: soln
    integer,        intent(in) :: num_iter
    
    integer :: i
    
    ! Repeat the following each time you want to write out the solution
    write(30,*) 'zone T="',num_iter,'" '
    write(30,*) 'I=',imax
    write(30,*) 'DATAPACKING=POINT'
    write(30,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE &
             & DOUBLE DOUBLE)'
    do i = 1, imax
    write(30,*) grid%xc(i),grid%Ac(i),soln%V(i,1),soln%V(i,2),soln%V(i,3),&
             & soln%mach(i),soln%U(i,1),soln%U(i,2),soln%U(i,3)
    enddo
    
  end subroutine output_soln




  !============================= output_res ==================================80
  !>
  !! Description: 
  !!
  !! Inputs:      grid     : 
  !!              soln     : 
  !!              num_iter : 
  !<
  !===========================================================================80
  subroutine output_res(rnorm,num_iter)
    
    real(prec), dimension(neq), intent(in) :: rnorm
    integer, intent(in) :: num_iter
    integer :: i
    ! Repeat the following each time you want to write out the solution
    write(40,*) num_iter,(rnorm(i),i=1,neq)
    
  end subroutine output_res
  
end module other_subroutines
