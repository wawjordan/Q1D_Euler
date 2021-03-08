module basic_boundaries

use set_precision, only : prec
use set_constants, only : one, two
use set_inputs,    only : neq, eps
use variable_conversion, only : prim2cons, isentropic_relations

implicit none

contains

  subroutine sub_in_bndry( M, U1, V1 )
    
    real(prec), dimension(:), intent(in) :: M
    real(prec), dimension(1,neq), intent(out) :: U1, V1
    real(prec), dimension(1)  :: M0, T1
    T1 = one
    M0 = two*M(1) - M(2)
    
    if (M0<eps) then
      M0 = eps
    end if
    
    call isentropic_relations(M0,V1,T1)
    call prim2cons(U1,V1)
    
  end subroutine sub_in_bndry


  subroutine sub_out_bndry( V, Un, Vn )
    
    use set_inputs, only : pb
    real(prec), dimension(:,:), intent(in) :: V
    real(prec), dimension(1,neq), intent(out) :: Un, Vn
    integer :: i
        
    i = size(V,1)
    
    Vn(1,1) = two*V(i,1) - V(i-1,1)
    Vn(1,2) = two*V(i,2) - V(i-1,2)
    Vn(1,3) = two*pb - V(i,3)
    
    call prim2cons(Un,Vn)
    
  end subroutine sub_in_bndry
  
  
  subroutine sup_out_bndry( V, Un, Vn )
    
    real(prec), dimension(:,:), intent(in) :: V
    real(prec), dimension(1,:), intent(out) :: Un, Vn
    integer :: i
     
    i = size(V,1)
    
    Vn = two*V(i,:) - V(i-1,:)
    
    call prim2cons(Un,Vn)
    
  end subroutine sub_in_bndry
  
  
  subroutine enforce_bndry()
    
    
    
  end subroutine sub_in_bndry
  
end module basic_boundaries
