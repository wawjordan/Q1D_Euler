module basic_boundaries

  use set_precision, only : prec
  use set_constants, only : one, two
  use set_inputs,    only : imax, neq, eps, i_low, i_high, ig_low, ig_high
  use variable_conversion, only : prim2cons, cons2prim, isentropic_relations, update_mach
  
  implicit none
  
  private
  
  public :: sub_in_bndry
  public :: sub_out_bndry, sup_out_bndry
  public :: enforce_bndry
  
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
    
    do i = i_low-1,-1,ig_low
      M(i) = 2*M(i+1) - M(i+2)
    end do
    
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
    
    real(prec), dimension(ig_low:ig_high,neq), intent(inout) :: U, V
    
    integer :: i
    
    do i = i_high+1,ig_high
      
      if (i == i_high+1) then
        V(i,1) = two*V(i-1,1) - V(i-2,1)
        V(i,2) = two*V(i-1,2) - V(i-2,2)
        V(i,3) = two*pb - V(i-1,3)
      else
        V(i,1) = two*V(i-1,1) - V(i-2,1)
        V(i,2) = two*V(i-1,2) - V(i-2,2)
        V(i,3) = two*V(i-1,3) - V(i-2,3)
      end if
      
    end do
    
    call prim2cons(U(i_high+1:ig_high,:),V(i_high+1:ig_high,:))
    call cons2prim(U(i_high+1:ig_high,:),V(i_high+1:ig_high,:))
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
    
    real(prec), dimension(ig_low:ig_high,neq), intent(inout) :: U, V
    
    integer :: i
    
      
    !do i = i_high+1,ig_high
    !  V(i,:) = two*V(i-1,:) - V(i-2,:)
    !end do
    do i = i_high+1,ig_high
      U(i,:) = two*U(i-1,:) - U(i-2,:)
      !U(i,:) = U(i-1,:)
    end do
    !U(i_high+1,:) = two*U(i_high,:)-U(i_high-1,:)
    !U(i_high+2,:) = U(i_high+1,:)
    !call prim2cons(U,V)
    !call prim2cons( U(i_high+1:ig_high,1:neq), V(i_high+1:ig_high,1:neq) )
    call cons2prim( U(i_high+1:ig_high,1:neq), V(i_high+1:ig_high,1:neq) )
    !call isentropic_relations(
  end subroutine sup_out_bndry
  
 subroutine enforce_bndry( soln )
  
  use set_inputs, only : shock
  use soln_type,  only : soln_t

  type( soln_t ), intent(inout) :: soln
  
  if (shock.eq.1) then
    call sub_in_bndry( soln%mach, soln%U, soln%V )
    call sub_out_bndry( soln%U, soln%V )
  elseif (shock.eq.0) then
    call sub_in_bndry( soln%mach, soln%U, soln%V )
    call sup_out_bndry( soln%U, soln%V )
  else
    write(*,*) 'ERROR! shock must equal 0 or 1!!!'
    stop
  end if

 end subroutine enforce_bndry
  
end module basic_boundaries
