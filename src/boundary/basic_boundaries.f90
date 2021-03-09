module basic_boundaries

use set_precision, only : prec
use set_constants, only : one, two
use set_inputs,    only : imax, neq, eps
use variable_conversion, only : prim2cons, isentropic_relations

implicit none

private

public :: sub_in_bndry
public :: sub_out_bndry, sup_out_bndry

contains

  subroutine sub_in_bndry( M, U, V )
    
    real(prec), dimension(:),   intent(inout) :: M
    real(prec), dimension(:,:), intent(inout) :: U, V
    !ireal(prec), dimension(:)  :: T
    integer :: i, i_low
    
    i_low  = lbound(M,1)
    
    if ( i_low < 1 ) then
      do i = 0,-1,i_low
        M(i) = 2*M(i+1) - M(i+2)
      end do
      call isentropic_relations(M(i_low:0),V(i_low:0,:))
      call prim2cons(U(i_low:0,:),V(i_low:0,:))
    end if
    
  end subroutine sub_in_bndry


  subroutine sub_out_bndry( U, V )
    
    use set_inputs, only : pb
    real(prec), dimension(:,:), intent(inout) :: U, V
    integer :: i, i_high
    
    i_high = ubound(V,1)
    
    if ( i_high > imax ) then
      do i = imax+1,i_high
        if (i == imax+1) then
          V(i,1) = two*V(i-1,1) - V(i-2,1)
          V(i,2) = two*V(i-1,2) - V(i-2,2)
          V(i,3) = two*pb - V(i-1,3)
        else
          V(i,1) = V(i-1,1)
          V(i,2) = V(i-1,2)
          V(i,3) = V(i-1,3)
        end if
      end do
      call prim2cons(U(imax+1:i_high,:),V(imax+1:i_high,:))
    end if
    
  end subroutine sub_out_bndry
  
  
  subroutine sup_out_bndry( U, V )
    
    use set_inputs, only : p0
    
    real(prec), dimension(:,:), intent(inout) :: U, V
    integer :: i, j, i_high
    
    i_high = ubound(V,1)
    
    if ( i_high > imax ) then
      do i = imax+1,i_high
        V(i,:) = two*V(i-1,:) - V(i-2,:)
        if ( V(i,1) < eps) then
            V(i,1) = eps
        elseif ( V(i,2) < eps ) then
            V(i,2) = eps
        elseif ( V(i,3)/(1000.0_prec*p0) < 1.0e-5_prec ) then
            V(i,3) = 0.01_prec*1000.0_prec*p0
        end if
      end do
      call prim2cons(U(imax+1:i_high,:),V(imax+1:i_high,:))
    end if
    
  end subroutine sup_out_bndry
  
 ! subroutine enforce_bndry()   
 ! end subroutine sub_in_bndry
  
end module basic_boundaries
