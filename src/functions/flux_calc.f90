module flux_calc
  
  use set_precision, only : prec
  use set_constants, only : zero, one, two, three, half, fourth
  use fluid_constants, only : gamma
  use set_inputs, only : neq, imax, i_low, i_high, ig_low, ig_high
  use variable_conversion, only : cons2prim, speed_of_sound
  
  implicit none
  
  private
  
  public :: flux_fun, select_flux
  
  procedure( calc_flux ), pointer :: flux_fun
   
  abstract interface
    
    subroutine calc_flux(left_state, right_state, F)
      
      import :: prec, i_low, i_high, neq
      real(prec), dimension(:,:), intent(in) :: left_state, right_state
      real(prec), dimension(i_low-1:i_high,neq), intent(out) :: F
      
    end subroutine calc_flux
    
  end interface
  
contains
  
  subroutine select_flux()
    
    use set_inputs, only : flux_scheme
    
    flux_fun => null()
    
    select case(flux_scheme)
    
    case(1)
      flux_fun => central_flux
    case(2)
      flux_fun => van_leer_flux
    !case(3)
    !  flux_fun => roe_flux
    case default
    
    end select
  
  end subroutine select_flux
  
  subroutine central_flux(left, right, F)
    
    real(prec), dimension(:,:), intent(in) :: left, right
    real(prec), dimension(i_low-1:i_high,neq), intent(out)   :: F
    real(prec), dimension(i_low-1:i_high,neq)   :: Ui
    
    Ui(i_low-1:i_high,:) = half*(left + right)
    
    F(:,1) = Ui(:,2)
    F(:,2) = half*(three-gamma)*( Ui(:,2)**2 )/Ui(:,1) &
           + (gamma-one)*Ui(:,3)
    F(:,3) = Ui(:,3)*Ui(:,2)/Ui(:,1) &
           + Ui(:,2)/Ui(:,1)*( (gamma-one)*Ui(:,3) &
           - half*(gamma-one)*Ui(:,2)**2/Ui(:,1) )
  
  end subroutine central_flux
  
  subroutine van_leer_flux(left, right, F)
    
    real(prec), dimension(:,:), intent(in) :: left, right
    real(prec), dimension(i_low-1:i_high,neq), intent(out)   :: F
    real(prec), dimension(i_low-1:i_high,neq)   :: VL, VR
    real(prec), dimension(i_low-1:i_high)   :: aL, aR
    real(prec), dimension(i_low-1:i_high)   :: ML, MR
    real(prec), dimension(i_low-1:i_high)   :: M_plus, M_minus
    real(prec), dimension(i_low-1:i_high)   :: alpha_plus, alpha_minus
    real(prec), dimension(i_low-1:i_high)   :: beta_L, beta_R
    real(prec), dimension(i_low-1:i_high)   :: c_plus, c_minus
    real(prec), dimension(i_low-1:i_high)   :: d_plus, d_minus
    real(prec), dimension(i_low-1:i_high)   :: p_plus, p_minus
    integer :: i
    
    VL = left
    VR = right
    
    call cons2prim(VL,VL)
    call cons2prim(VR,VR)
    
    !call cons2prim(left,VL)
    !call cons2prim(right,VR)
    !do i = i_low-1,i_high
    !  write(*,*) i, VL(i,1), VL(i,2), VL(i,3)
    !end do
    call speed_of_sound(VL(:,3),VL(:,1),aL)
    call speed_of_sound(VR(:,3),VR(:,1),aR)
    
    ML = VL(:,2)/aL
    MR = VR(:,2)/aR
    M_plus  = fourth*(ML+one)**2
    M_minus = -fourth*(MR-one)**2
    
    !write(*,*)
    !write(*,*) 'i      ','ML      ','MR      ','Mplus      ','Mminus      '
    !do i = i_low-1,i_high
    !  write(*,*) i, ML(i), MR(i), M_plus(i), M_minus(i)
    !end do
    
    beta_L = -max(zero,one-int(abs(ML)))
    beta_R = -max(zero,one-int(abs(MR)))
    alpha_plus  = half*(one+sign(one,ML))
    alpha_minus = half*(one-sign(one,MR))
    c_plus  = alpha_plus*(one+beta_L)*ML - beta_L*M_plus
    c_minus = alpha_minus*(one+beta_R)*MR - beta_R*M_minus
    
    p_plus  = M_plus*(- ML + two)
    p_minus = M_minus*(- MR - two)
    d_plus  = alpha_plus*(one+beta_L) - beta_L*p_plus
    d_minus = alpha_minus*(one+beta_R) - beta_R*p_minus
    
    do i = i_low-1,i_high
      write(*,*) i, d_plus(i), d_minus(i)
    end do
    
    F(:,1) = VL(:,1)*aL*c_plus + VR(:,1)*aR*c_minus
    F(:,2) = VL(:,1)*aL*c_plus*VL(:,2) + VR(:,1)*aR*c_minus*VR(:,2) + &
          & d_plus*VL(:,3) + d_minus*VR(:,3)
    F(:,3) = VL(:,1)*aL*c_plus*( aL**2/(gamma-1) + half*VL(:,2)**2) + &
          & VR(:,1)*aR*c_minus*( aR**2/(gamma-1) + half*VR(:,2)**2)
    !write(*,*)
    !write(*,*)
    !do i = i_low-1,i_high
    !  write(*,*) i, F(i,1), F(i,2), F(i,3)
    !end do
    
  end subroutine van_leer_flux
  
end module flux_calc
