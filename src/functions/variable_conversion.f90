module variable_conversion
   
  use set_precision,   only : prec
  use set_constants,   only : one, half
  use fluid_constants, only : gamma, R_gas
  use set_inputs,      only : p0,T0, neq, ig_low, ig_high

  implicit none

  contains
  
  !============================== speed_of_sound  ============================80
  !>
  !! Description:
  !!
  !! Inputs:      pressure :
  !!              rho      :
  !!
  !! Outputs:     sound_speed :
  !<
  !===========================================================================80
  subroutine speed_of_sound( pressure, rho, sound_speed, i_start, i_stop )
    integer, intent(in) :: i_start, i_stop
    real(prec), dimension(:), intent(in)  :: pressure, rho
    real(prec), dimension(:), intent(out) :: sound_speed
    
    sound_speed = sqrt(gamma*pressure/rho)
    
  end subroutine speed_of_sound
  
  !=============================== prim2cons =================================80
  !>
  !! Description: 
  !!
  !! Inputs:      U :
  !!
  !! Outputs:     U : 
  !!              V :
  !<
  !===========================================================================80
  subroutine prim2cons( U, V, i_start, i_stop )
    
    integer, intent(in) :: i_start, i_stop
    real(prec), dimension(ig_low:ig_high,neq), intent(inout)  :: U
    real(prec), dimension(ig_low:ig_high,neq), intent(in)     :: V
    
    U(:,1) = V(:,1)
    U(:,2) = V(:,1)*V(:,2)
    U(:,3) = ( V(:,3)/(gamma - one) ) + half*V(:,1)*V(:,2)**2
    
  end subroutine prim2cons
  
  !================================ update_mach ==============================80
  !>
  !! Description: 
  !!
  !! Inputs:      V : 
  !!              M : 
  !!
  !! Outputs:     M :
  !<
  !===========================================================================80
  subroutine update_mach( V, M, i_start, i_stop )
    
    integer, intent(in) :: i_start, i_stop
    real(prec), dimension(i_start:i_stop,neq),  intent(in) :: V 
    real(prec), dimension(i_start:i_stop), intent(inout) :: M
    
    real(prec), dimension(i_start:i_stop) :: a
    
    call speed_of_sound(V(:,3),V(:,1),a(:),i_start,i_stop)
    
    M = abs(V(:,2))/a(:)
    
  end subroutine update_mach
  
  !=============================== cons2prim =================================80
  !>
  !! Description: 
  !!
  !! Inputs:      U : 
  !!              V : 
  !!
  !! Outputs:     U : 
  !!              V : 
  !<
  !===========================================================================80
  subroutine cons2prim( U, V, i_start, i_stop )
    
    integer, intent(in) :: i_start, i_stop
    real(prec), dimension(i_start:i_stop,neq), intent(inout) :: U
    real(prec), dimension(i_start:i_stop,neq), intent(inout) :: V
    
    V(:,1) = U(:,1)
    V(:,2) = U(:,2)/U(:,1)
    V(:,3) = (gamma - one)*U(:,3) - half*(gamma - one)*U(:,2)**2/U(:,1)
    
    call limit_primitives(U,V,i_start,i_stop)
    
  end subroutine cons2prim
  
  !========================= limit_primitives ================================80
  !>
  !! Description: 
  !!
  !! Inputs:      U : 
  !!              V : 
  !!
  !! Outputs:     U : 
  !!              V : 
  !<
  !===========================================================================80
  subroutine limit_primitives(U,V, i_start, i_stop)
    
    integer, intent(in) :: i_start, i_stop
    real(prec), dimension(i_start:i_stop,neq), intent(inout) :: U
    real(prec), dimension(i_start:i_stop,neq), intent(inout) :: V
    integer :: i
    !logical, dimension(lbound(U,1):ubound(U,1)) :: mask
    
    !mask = .false.
    
    !where ( (V(:,1)<0.001_prec).or.(V(:,2)<10.0_prec).or.(V(:,3)<500.0_prec) )
    !  mask = .true.
    !end where
    
    do i = i_start,i_stop
      if (V(i,1)<0.001_prec) then
        V(i,1) = 0.001_prec
      end if
      if (V(i,2)<10.0_prec) then
        V(i,2) = 10.0_prec
      end if
      if (V(i,3)<500.0_prec) then
        V(i,3) = 500.0_prec
      end if
      !if (mask(i)) then
      !  call prim2cons(U,V)
      !end if
    end do
    call prim2cons(U,V,i_start,i_stop)
    
  end subroutine limit_primitives
  
  !========================= isentropic_relations ============================80
  !>
  !! Description: 
  !!
  !! Inputs:      M : 
  !!              V : 
  !!
  !! Outputs:     V : 
  !<
  !===========================================================================80
  subroutine isentropic_relations(M,V,i_start,i_stop)
    
    integer, intent(in) :: i_start, i_stop
    real(prec), dimension(i_start:i_stop,neq), intent(inout) :: V
    real(prec), dimension(i_start:i_stop),   intent(in) :: M
    
    real(prec), dimension(i_start:i_stop)  :: T
    
    T(:) = T0/(one + half*(gamma - one)*M(:)**2)
    V(:,3) = 1000.0_prec*p0/(one + half*(gamma - one)*M(:)**2)**(gamma/(gamma-1))
    V(:,1) = V(:,3)/(R_gas*T(:))
    call speed_of_sound(V(:,3),V(:,1),V(:,2),i_start,i_stop)
    V(:,2) = M(:)*V(:,2)
    
  end subroutine isentropic_relations
  
end module variable_conversion
