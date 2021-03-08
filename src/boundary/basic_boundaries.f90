module basic_boundaries

use set_precision, only : prec
use set_constants, only : one, two
use set_inputs,    only : neq, eps
use variable_conversion, only : prim2cons, isentropic_relations

implicit none

contains

  subroutine sub_in_bndry( M, U1, V1 )
    
    real(prec), dimension(:), intent(in) :: M
    real(prec), dimension(2,neq), intent(out) :: U1, V1
    real(prec), dimension(2)  :: M0, T1
    T1 = one
    M0(1) = two*M(1) - M(2)
    if (M0(1)<eps) then
      M0(1) = eps
    end if
    M0(2) = two*M0(1) - M(2)
    if (M0(2)<eps) then
      M0(2) = eps
    end if
    
    call isentropic_relations(M0,V1,T1)
    call prim2cons(U1,V1)
    
  end subroutine sub_in_bndry


  subroutine sub_out_bndry( V, Un, Vn )
    
    use set_inputs, only : pb
    real(prec), dimension(:,:), intent(in) :: V
    real(prec), dimension(2,neq), intent(out) :: Un, Vn
    integer :: i
        
    i = size(V,1)
    
    Vn(1,1) = two*V(i,1) - V(i-1,1)
    Vn(1,2) = two*V(i,2) - V(i-1,2)
    Vn(1,3) = two*pb - V(i,3)

    Vn(2,1) = two*Vn(1,1) - V(i,1)
    Vn(2,2) = two*Vn(1,2) - V(i,2)
    Vn(2,3) = two*Vn(1,3) - pb
    
    call prim2cons(Un,Vn)
    
  end subroutine sub_in_bndry
  
  
  subroutine sup_out_bndry( V, Un, Vn )
    
    real(prec), dimension(:,:), intent(in) :: V
    real(prec), dimension(2,:), intent(out) :: Un, Vn
    integer :: i
     
    i = size(V,1)
    
    Vn(1,:) = two*V(i,:) - V(i-1,:)
    Vn(2,:) = two*Vn(1,:) - V(i,:)
    
    call prim2cons(Un,Vn)
    
  end subroutine sub_in_bndry
  
  
  subroutine enforce_bndry()
    
    
    
  end subroutine sub_in_bndry
  
end module basic_boundaries
