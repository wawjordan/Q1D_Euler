module fluxes
  
  use set_precision, only : prec
  use set_constants, only : one, three, half
  use fluid_constants, only : gamma
  use set_inputs, only : neq, i_low, i_high, ig_low, ig_high
  
  implicit none
  
contains
  
subroutine central_flux(U,F)
  
  real(prec), dimension(ig_low:ig_high,neq), intent(inout)  :: U 
  real(prec), dimension(i_low-1:i_high,neq), intent(inout)  :: F
  real(prec), dimension(i_low-1:i_high,neq)  :: Ui

  Ui(i_low-1:i_high,:) = half*(U(i_low:i_high+1,:) + U(i_low-1:i_high,:))
  
  F(:,1) = Ui(:,2)
  F(:,2) = half*(three-gamma)*( Ui(:,2)**2 )/Ui(:,1) + (gamma-one)*Ui(:,3)
  F(:,3) = Ui(:,3)*Ui(:,2)/Ui(:,1) + Ui(:,2)/Ui(:,1)*( (gamma-one)*Ui(:,3) - half*(gamma-one)*Ui(:,2)**2/Ui(:,1) )
  
end subroutine central_flux

end module fluxes
