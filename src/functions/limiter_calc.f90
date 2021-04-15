module limiter_calc
  
  use set_precision, only : prec
  use set_constants, only : zero, one, two, three, four, half, fourth
  use set_inputs, only : neq, imax, i_low, i_high, ig_low, ig_high
  
  implicit none
  
  private
  
  public :: limiter_fun, select_limiter
  
  procedure( calc_limiter ), pointer :: limiter_fun
   
  abstract interface
    
    subroutine calc_limiter( r, psi )
      
      import :: prec, i_low, i_high, neq
      real(prec), dimension(:,:), intent(in) :: r
      real(prec), dimension(i_low-1:i_high,neq), intent(out) :: psi
      
    end subroutine calc_limiter
    
  end interface
  
contains
  
  subroutine select_limiter()
    
    use set_inputs, only : limiter_scheme
    
    limiter_fun => null()
    
    select case(limiter_scheme)
    
    case(1)
      limiter_fun => van_leer_limiter
    case(2)
      limiter_fun => van_albada_limiter
    case(3)
      limiter_fun => minmod_limiter
    case(4)
      limter_fun => beta_limiter
    case default
    
    end select
  
  end subroutine select_limiter
  
  subroutine van_leer_limiter( r, psi )
    
    real(prec), dimension(:,:), intent(in) :: left, right
    real(prec), dimension(i_low-1:i_high,neq), intent(out)   :: psi
    
    psi = (r + abs(r))/(one + r)
    
  end subroutine van_leer_limiter
  
  subroutine van_albada_limiter( r, psi )
    
    real(prec), dimension(:,:), intent(in) :: r
    real(prec), dimension(i_low-1:i_high,neq), intent(out)   :: psi
    
    psi = (r**2 + r)/(one + r**2)
    
  end subroutine van_albada_limiter

  subroutine minmod_limiter( r, psi )
    
    real(prec), dimension(:,:), intent(in) :: r
    real(prec), dimension(i_low-1:i_high,neq), intent(out)   :: psi
    
    psi = half*(one + sign(one,r))*min(r,one)
    
  end subroutine minmod_limiter
  
  subroutine beta_limiter( r, psi )
    
    real(prec), dimension(:,:), intent(in) :: r
    real(prec), dimension(i_low-1:i_high,neq), intent(out)   :: psi
    
    psi = maxval((/ zero, min(beta*r,one), min(r,beta) /))
    
  end subroutine beta_limiter
  
end module limiter_calc
