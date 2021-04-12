module fluxes
  
  use set_precision, only : prec
  use set_constants, only : one, three, half, fourth
  use fluid_constants, only : gamma
  use set_inputs, only : neq, i_low, i_high, ig_low, ig_high
  
  implicit none
  
  contains
   
  !================================= central_flux ============================80
  !>
  !! Description: 
  !!
  !! Inputs:      U : 
  !!              F : 
  !!
  !! Outputs:     U :
  !!              F :
  !<
  !===========================================================================80
  subroutine central_flux(U,F)
     
    real(prec), dimension(i_low-1:i_high,neq), intent(out)  :: F
    real(prec), dimension(i_low-1:i_high,neq)   :: Ui 
    real(prec), dimension(ig_low:ig_high,neq), intent(in)  :: U
    
    Ui(i_low-1:i_high,:) = half*(U(i_low:i_high+1,:) + U(i_low-1:i_high,:))
    
    F(:,1) = Ui(:,2)
    F(:,2) = half*(three-gamma)*( Ui(:,2)**2 )/Ui(:,1) &
           + (gamma-one)*Ui(:,3)
    F(:,3) = Ui(:,3)*Ui(:,2)/Ui(:,1) &
           + Ui(:,2)/Ui(:,1)*( (gamma-one)*Ui(:,3) &
           - half*(gamma-one)*Ui(:,2)**2/Ui(:,1) )
    
  end subroutine central_flux
  
  subroutine van_leer1(V,F,asnd,mach)
  
    real(prec), dimension(i_low-1:i_high,neq), intent(out) :: F
    real(prec), dimension(ig_low:ig_high,neq), intent(in)  :: V
    real(prec), dimension(ig_low:ig_high), intent(in)  :: asnd
    real(prec), dimension(ig_low:ig_high), intent(in)  :: mach
    real(prec), dimension(i_low-1:i_high)  :: c_plus, c_minus, d_plus, d_minus
    real(prec) :: Fc, Fp
    real(prec) :: M_plus, M_minus, p_plus, p_minus
    real(prec) :: beta_left, beta_right, alpha_plus, alpha_minus
    
    integer :: i
    
    do i = i_low-1,i_high
      M_plus  =  fourth*(mach(i)+1)**2
      M_minus = -fourth*(mach(i+1)-1)**2
      beta_left  = -max(0,one-int(abs(mach(i))))
      beta_right = -max(0,one-int(abs(mach(i+1))))
      alpha_plus  = half*(1+sign(one,mach(i)))
      alpha_minus = half*(1-sign(one,mach(i+1)))
      c_plus(i)  = alpha_plus*(1+beta_left)*mach(i)-beta_left*M_plus
      c_minus(i) = alpha_minus*(1+beta_right)*mach(i+1)-beta_right*M_minus
      p_plus  = M_plus*(-mach(i)+2)
      p_minus = M_minus*(-mach(i+1)-2)
      d_plus(i)  = alpha_plus*(1+beta_left) - beta_left*p_plus
      d_minus(i) = alpha_minus*(1+beta_right) - beta_right*p_minus
    end do
    F(:,1) = V(i_low-1:i_high,1)*asnd(i_low-1:i_high)*c_plus + &
         & + V(i_low:i_high+1,1)*asnd(i_low:i_high+1)*c_minus
    F(:,2) = V(i_low-1:i_high,1)*asnd(i_low-1:i_high)*c_plus*V(i_low-1:i_high,2) + &
         & + V(i_low:i_high+1,1)*asnd(i_low:i_high+1)*c_minus*V(i_low:i_high+1,2) 
    F(:,3) = V(i_low-1:i_high,1)*asnd(i_low-1:i_high)*c_plus + &
         & + V(i_low:i_high+1,1)*asnd(i_low:i_high+1)*c_minus
  end subroutine van_leer1

end module fluxes
