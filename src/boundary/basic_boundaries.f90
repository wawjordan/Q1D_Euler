module basic_boundaries

  use set_precision, only : prec
  use set_constants, only : one, two
  use set_inputs,    only : imax, neq, eps, i_low, i_high, ig_low, ig_high
  use variable_conversion, only : prim2cons, isentropic_relations, update_mach
  
  implicit none
  
  private
  
  public :: sub_in_bndry
  public :: sub_out_bndry, sup_out_bndry

  contains
  
  !============================= sub_in_bndry  ===============================80
  !>
  !! Description: 
  !!
  !! Inputs:      M : 
  !!              U : 
  !!              V :
  !!
  !! Outputs:     M :
  !!              U : 
  !!              V :
  !<
  !===========================================================================80
  subroutine sub_in_bndry( M, U, V )
    
    real(prec), dimension(ig_low:ig_high), intent(inout) :: M
    real(prec), dimension(ig_low:ig_high,1:neq), intent(inout) :: U, V
    !real(prec), intent(inout) :: M(ig_low:ig_high)
    !real(prec), intent(inout) :: U(ig_low:ig_high,:), V(ig_low:ig_high,:)
    
    integer :: i
    
    call update_mach(V,M)
    
    do i = i_low-1,-1,ig_low
      M(i) = 2*M(i+1) - M(i+2)
    end do
    
    !call isentropic_relations(M,V)
    !call prim2cons(U,V)
    call isentropic_relations(M(ig_low:i_low-1),V(ig_low:i_low-1,:))
    call prim2cons(U(ig_low:i_low-1,:),V(ig_low:i_low-1,:))
    
  end subroutine sub_in_bndry
  
  !============================ sub_out_bndry ================================80
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
  subroutine sub_out_bndry( U, V )
    
    use set_inputs, only : pb
    
    real(prec), dimension(:,:), intent(inout) :: U, V
    
    integer :: i, i_high
    
    i_high = ubound(V,1)
    
    do i = i_high+1,ig_high
      
      if (i == i_high+1) then
        V(i,1) = two*V(i-1,1) - V(i-2,1)
        V(i,2) = two*V(i-1,2) - V(i-2,2)
        V(i,3) = two*pb - V(i-1,3)
      else
        V(i,1) = V(i-1,1)
        V(i,2) = V(i-1,2)
        V(i,3) = V(i-1,3)
      end if
      
    end do
    
    call prim2cons(U(i_high+1:ig_high,:),V(i_high+1:ig_high,:))
      
  end subroutine sub_out_bndry
  
  !============================= sup_out_bndry ===============================80
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
  subroutine sup_out_bndry( U, V )
    
    use set_inputs, only : p0
    
    real(prec), dimension(ig_low:ig_high,neq), intent(inout) :: U, V
    
    integer :: i
    
      
    do i = i_high+1,ig_high
      V(i,:) = two*V(i-1,:) - V(i-2,:)
      
      !if ( V(i,1) < eps) then
      !    V(i,1) = eps
      !elseif ( V(i,2) < eps ) then
      !    V(i,2) = eps
      !elseif ( V(i,3)/(1000.0_prec*p0) < 1.0e-5_prec ) then
      !    V(i,3) = 0.01_prec*1000.0_prec*p0
      !end if
      
    end do
    !call prim2cons(U,V)
    call prim2cons( U(i_high+1:ig_high,:), V(i_high+1:ig_high,:) )
    
  end subroutine sup_out_bndry
  
 ! subroutine enforce_bndry()   
 ! end subroutine sub_in_bndry
  
end module basic_boundaries
