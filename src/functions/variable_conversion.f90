module variable_conversion
   
  use set_precision,   only : prec
  use set_constants,   only : one, half
  use fluid_constants, only : gamma, R_gas
  use set_inputs,      only : p0,T0

  implicit none

contains

subroutine speed_of_sound(pressure,rho,sound_speed)

  real(prec), dimension(:), intent(in)  :: pressure, rho
  real(prec), dimension(:), intent(out) :: sound_speed
  
  sound_speed = sqrt(gamma*pressure/rho)

end subroutine speed_of_sound


subroutine prim2cons(U,V)

  real(prec), dimension(:,:), intent(inout)  :: U
  real(prec), dimension(:,:), intent(in)     :: V

  U(:,1) = V(:,1)
  U(:,2) = V(:,1)*V(:,2)
  U(:,3) = ( V(:,3)/(gamma - one) ) + half*V(:,1)*V(:,2)**2

end subroutine prim2cons


subroutine cons2prim(U,V)
  
  real(prec), dimension(:,:), intent(inout) :: U
  real(prec), dimension(:,:), intent(inout) :: V

  V(:,1) = U(:,1)
  V(:,2) = U(:,2)/U(:,1)
  V(:,3) = (gamma - one)*U(:,3) - half*(gamma - one)*U(:,2)**2/U(:,1)
  
  call limit_primitives(U,V)
  
end subroutine cons2prim


subroutine limit_primitives(U,V)
  
  real(prec), dimension(:,:), intent(inout) :: U
  real(prec), dimension(:,:), intent(inout) :: V
  integer :: i
  logical, dimension(lbound(U,1):ubound(U,1)) :: mask
  
  mask = .false.
  
  where ( (V(:,1)<0.001_prec).or.(V(:,2)<10.0_prec).or.(V(:,3)<500.0_prec) )
    mask = .true.
  end where
  
  do i = lbound(U,1),ubound(U,1)
    if (V(i,1)<0.001_prec) then
      V(i,1) = 0.001_prec
    end if
    if (V(i,2)<10.0_prec) then
      V(i,2) = 10.0_prec
    end if
    if (V(i,3)<500.0_prec) then
      V(i,3) = 500.0_prec
    end if
    write(*,*) 'mask: ',mask(i)
    if (mask(i)) then
      call prim2cons(U,V)
    end if
  end do
  
end subroutine limit_primitives

subroutine isentropic_relations(M,V)

  real(prec), dimension(:,:), intent(inout) :: V
  real(prec), dimension(:),   intent(in) :: M
  !real(prec), dimension(:),   intent(out) :: T
  real(prec), dimension(size(V,1))  :: T
  
  T(:) = T0/(one + half*(gamma - one)*M(:)**2)
  V(:,3) = 1000.0_prec*p0/(one + half*(gamma - one)*M(:)**2)**(gamma/(gamma-1))
  V(:,1) = V(:,3)/(R_gas*T(:))
  call speed_of_sound(V(:,3),V(:,1),V(:,2))
  V(:,2) = M(:)*V(:,2)

end subroutine isentropic_relations


end module variable_conversion
