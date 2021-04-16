module other_subroutines
  
  use set_precision, only : prec
  use set_constants, only : zero, one, two, three, half, fourth
  use set_inputs, only : imax, neq, i_low, i_high, ig_low, ig_high, shock
  use set_inputs, only : epsM, kappaM
  use fluid_constants, only : gamma
  use variable_conversion
  use limiter_calc, only : limiter_fun
  use soln_type, only : soln_t
  use exact_q1d_type, only : exact_q1d_t
  use grid_type, only : grid_t
  
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
    !integer :: i
    !S(i_low:i_high) = V(i_low:i_high,3)*dA(i_low:i_high)/1000.0_prec
    S(i_low:i_high) = P(i_low:i_high)*dA(i_low:i_high)
    !do i = i_low,i_high
    !write(*,*) i, "dA : ",dA(i), "V: ",V(i,:)," S: ",S(i)
    !end do
  end subroutine calculate_sources
  
  subroutine MUSCL_extrap( V, left, right )
    
    real(prec), dimension(ig_low:ig_high,neq), intent(in)  :: V
    real(prec), dimension(i_low-1:i_high,neq), intent(out) :: left, right
    real(prec), dimension(i_low-1:i_high,neq) :: r_plus, r_minus
    real(prec), dimension(i_low-1:i_high,neq) :: psi_plus, psi_minus
    real(prec), dimension(neq) :: den
    integer :: i
    
    do i = i_low-1,i_high
      den = V(i+1,:) - V(i,:)
      den = sign(one,den)*max(abs(den),1e-6_prec)
      r_plus(i,:)   = ( V(i+2,:) - V(i+1,:) )/den
      r_minus(i,:)  = ( V(i,:) - V(i-1,:) )/den
      !write(*,*) i, r_plus(i,1), r_plus(i,2), r_plus(i,3), &
      !         &    r_minus(i,1),r_minus(i,2), r_minus(i,3)
    end do
    !write(*,*)
    !stop
    !r_plus  = ( V(i_low+1:i_high+2,:) - V(i_low:i_high+1,:) )/&
    !        & ( V(i_low:i_high+1,:) - V(i_low-1:i_high,:) )
    !r_minus = ( V(i_low-1:i_high,:) - V(i_low-2:i_high-1,:) )/&
    !        & ( V(i_low:i_high+1,:) - V(i_low-1:i_high,:) )
    
    call limiter_fun(r_plus,psi_plus)
    call limiter_fun(r_minus,psi_minus)
    
    do i = i_low,i_high-1
      left(i,:) = V(i,:) + fourth*epsM*( &
         & (one-kappaM)*psi_plus(i-1,:)*(V(i,:)-V(i-1,:)) + &
         & (one+kappaM)*psi_minus(i,:)*(V(i+1,:)-V(i,:)) )
      right(i,:) = V(i+1,:) - fourth*epsM*( &
         & (one+kappaM)*psi_minus(i+1,:)*(V(i+1,:)-V(i,:)) + &
         & (one-kappaM)*psi_plus(i,:)*(V(i+2,:)-V(i+1,:)) )
    end do
    !left(i_low-1,:)  = two*left(i_low,:) - left(i_low+1,:)
    !left(i_high,:)   = two*left(i_high-1,:) - left(i_high-2,:)
    !right(i_low-1,:) = two*right(i_low,:) - right(i_low+1,:)
    !right(i_high,:)  = two*right(i_high-1,:) - right(i_high-2,:)
    left(i_low-1,:)  = V(i_low-1,:)
    left(i_high,:)   = V(i_high,:)
    right(i_low-1,:) = V(i_low,:)
    right(i_high,:)  = V(i_high+1,:)
    !write(*,*)
    !do i = i_low-1,i_high
    !  write(*,*) i, left(i,1), left(i,2), left(i,3), right(i,1), right(i,2), right(i,3)
    !end do

  end subroutine MUSCL_extrap
  
  
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
  
  subroutine calc_de( soln, exact_soln, DE, DEnorm, pnorm )
    
    type(soln_t), intent(inout) :: soln
    type(exact_q1d_t), intent(in) :: exact_soln
    real(prec), dimension(i_low:i_high,1:neq), intent(out) :: DE
    real(prec), dimension(1,1:neq), intent(out) :: DEnorm
    integer, intent(in) :: pnorm
    real(prec) :: Linv
    Linv = one/real(i_high-i_low)
    DE = soln%V(i_low:i_high,1:neq) - exact_soln%Vc(i_low:i_high,1:neq)
    
    if (pnorm == 0) then
      DEnorm(1,1:neq) = maxval(abs(DE),1)
    elseif (pnorm == 1) then
      DEnorm(1,1) = Linv*sum(abs(DE(:,1)),1)
      DEnorm(1,2) = Linv*sum(abs(DE(:,2)),1)
      DEnorm(1,3) = Linv*sum(abs(DE(:,3)),1)
    elseif (pnorm == 2) then
      DEnorm(1,1) = sqrt(Linv*sum(DE(:,1)**2,1))
      DEnorm(1,2) = sqrt(Linv*sum(DE(:,2)**2,1))
      DEnorm(1,3) = sqrt(Linv*sum(DE(:,3)**2,1))
    else
      DEnorm(1,1) = maxval(abs(DE(:,1)))
      DEnorm(1,2) = maxval(abs(DE(:,2)))
      DEnorm(1,3) = maxval(abs(DE(:,3)))
    end if
    
  end subroutine calc_de

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

    use set_inputs, only : imax, CFL, k2, k4, shock, ramp, p_ratio

    character(len=1024) :: dirname
    character(len=1024) :: filename
    character(len=64) ::   shock_str
    character(len=64) ::   ncells_str
    character(len=64) ::   CFL_str
    character(len=64) ::   kappa2_str
    character(len=64) ::   kappa4_str

    ! Set up output directories
    write (ncells_str  , "(A1,I0.3,A1)") "N"  , imax   , "_"
    if (shock.eq.1) then
      write (shock_str, "(A16,I0.2)") "normal-shock-pb-",int(100*p_ratio)
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
    !write (dirname, *) adjustl(trim(shock_str)),"/"//  &
    !&                   adjustl(trim(ncells_str))
    write (dirname, *) adjustl(trim(shock_str)),"/"
    write (filename,*)  trim(ncells_str)//  &
    &                   trim(CFL_str)//     &
    &                   trim(kappa2_str)//  &
    &                   trim(kappa4_str)
    
    call execute_command_line ('mkdir -p ../results/' // adjustl(trim(dirname)))

  ! Set up output files (history and solution)
  ! open(30,file='history.dat',status='unknown')
    open(30,file= '../results/'//trim(adjustl(dirname))//  &
    &               trim(adjustl(filename))//'_history.dat',status='unknown')
    write(30,*) 'TITLE = "Quasi-1D Nozzle Iterative Residual History"'
    if(shock.eq.0) then
      write(30,*) 'variables="Iteration""Res1""Res2""Res3""DE1""DE2""DE3"'
    else
      write(30,*) 'variables="Iteration""Res1""Res2""Res3"'
    end if


    open(40,file= '../results/'//trim(adjustl(dirname))//  &
    &               trim(adjustl(filename))//'_field.dat',status='unknown')
    write(40,*) 'TITLE = "Quasi-1D Nozzle Solution"'
    if(shock.eq.0) then
      write(40,*) 'variables="x(m)""A(m^2)""rho(kg/m^3)""u(m/s)""p(N/m^2)"  &
      & "M""U1""U2""U3""M-exact""rho-exact""u-exact""p-exact""DE1""DE2""DE3"'
    elseif(shock.eq.1) then
      write(40,*) 'variables="x(m)""A(m^2)""rho(kg/m^3)""u(m/s)""p(N/m^2)"&
      &"M""U1""U2""U3"'
    else
      write(*,*) 'ERROR! shock must equal 0 or 1!!!'
      stop
    endif
  
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
  subroutine output_soln(grid,soln,ex_soln,num_iter)
    
    type( grid_t ), intent(in) :: grid
    type( soln_t ), intent(in) :: soln
    type( exact_q1d_t ), intent(in) :: ex_soln
    integer,        intent(in) :: num_iter
    
    integer :: i


    open(40,status='unknown')
    !write(40,*) 'TITLE = "Quasi-1D Nozzle Solution"'
    !write(40,*)' variables="x(m)""Area(m^2)""rho(kg/m^3)""u(m/s)"&
    !         & "Press(N/m^2)""Mach""U1""U2""U3""F1""F2""F3"'
    ! Repeat the following each time you want to write out the solution
    !write(40,*) 'zone T="',num_iter,'" '
    !write(40,*) 'I=',imax
    !write(40,*) 'DATAPACKING=POINT'
    !write(40,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE &
    !         & DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)'
    !do i = 1, imax
    !write(40,*)grid%xc(i),grid%Ac(i),soln%V(i,1),soln%V(i,2),soln%V(i,3),&
    !         & soln%mach(i),soln%U(i,1),soln%U(i,2),soln%U(i,3),&
    !         & soln%F(i,1),soln%F(i,2),soln%F(i,3)
    !enddo

    
    ! Repeat the following each time you want to write out the solution
    write(40,*) 'zone T="',num_iter,'" '
    write(40,*) 'I=',imax
    if(shock.eq.0) then
      write(40,*) 'DATAPACKING=POINT'
      write(40,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE &
             & DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE &
             & DOUBLE)'
      do i = i_low,i_high
        write(40,*) grid%xc(i),grid%Ac(i),soln%V(i,1),soln%V(i,2),soln%V(i,3),&
             & soln%mach(i),soln%U(i,1),soln%U(i,2),soln%U(i,3), &
             & ex_soln%Mc(i), ex_soln%Vc(i,1), ex_soln%Vc(i,2),  &
             & ex_soln%Vc(i,3), soln%DE(i,1), soln%DE(i,2), soln%DE(i,3)
      end do
    elseif(shock.eq.1) then
      write(40,*) 'DATAPACKING=POINT'
      write(40,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE &
             & DOUBLE DOUBLE)'
      do i = i_low,i_high
        write(40,*) grid%xc(i),grid%Ac(i),soln%V(i,1),soln%V(i,2),soln%V(i,3),&
             & soln%mach(i),soln%U(i,1),soln%U(i,2),soln%U(i,3)
      end do
    endif
    
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
  subroutine output_res(soln,num_iter)
    
    !real(prec), dimension(neq), intent(in) :: rnorm
    type( soln_t), intent(in) :: soln
    integer, intent(in) :: num_iter
    integer :: i
    ! Repeat the following each time you want to write out the solution
    if(shock.eq.0) then
      write(30,*) num_iter,(soln%rnorm(i),i=1,neq),(soln%DEnorm(i),i=1,neq)
    else
      write(30,*) num_iter,(soln%rnorm(i),i=1,neq)
    end if
    
  end subroutine output_res
  
end module other_subroutines
