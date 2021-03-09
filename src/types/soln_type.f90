module soln_type

  use set_precision, only : prec
  use set_constants, only : one
  
  implicit none

  type soln_t
    
    real(prec), allocatable, dimension(:,:) :: V
    real(prec), allocatable, dimension(:,:) :: U
    real(prec), allocatable, dimension(:,:) :: F
    real(prec), allocatable, dimension(:,:) :: d
    real(prec), allocatable, dimension(:)   :: M
    real(prec), allocatable, dimension(:)   :: T
    real(prec), allocatable, dimension(:,:) :: R
    real(prec), allocatable, dimension(:)   :: dt
    !real(prec), allocatable, dimension(:,:) :: eps
    real(prec), allocatable, dimension(:)   :: lambda

  end type soln_t

contains

  subroutine allocate_soln( soln )

    use set_constants, only : zero
    use set_inputs,    only : imax, neq, n_ghost_cells
    
    implicit none
    
    type(soln_t), intent(inout) :: soln
    integer :: i_low1, i_high1
    integer :: i_low2, i_high2
    
    i_low1  = 1 - n_ghost_cells
    i_high1 = imax + n_ghost_cells
    i_low2  = 1
    i_high2 = imax

    allocate(soln%U(     i_low1 : i_high1,neq),      &
             soln%V(     i_low1 : i_high1,neq),      &
             soln%M(     i_low1 : i_high1),      &
             soln%T(     i_low1 : i_high1),      &
             soln%dt(    i_low1 : i_high1),     &
             soln%lambda(i_low1 : i_high1),     &
             soln%R(     i_low2 : i_high2,neq),      &
             soln%F(     i_low2-1 : i_high2,neq),      &
             soln%D(     i_low2-1 : i_high2,neq),      )


!    allocate(soln%U(i_low:i_high,neq),      &
!             soln%V(i_low:i_high,neq),      &
!             soln%F(i_low:i_high,neq),      &
!             soln%D(i_low:i_high,neq),      &
!             soln%M(i_low:i_high),      &
!             soln%T(i_low:i_high),      &
!             soln%R(i_low:i_high,neq),      &
!             soln%dt(i_low:i_high),     &
!             soln%eps(i_low:i_high,neq),    &
!             soln%lambda(i_low:i_high)  )

    soln%U   = zero
    soln%V   = zero
    soln%F   = zero
    soln%D   = zero
    soln%M   = zero
    soln%T   = zero
    soln%R   = zero
    soln%dt  = zero
    soln%lambda = zero

  end subroutine allocate_soln
  
  subroutine deallocate_soln( soln )
  
    implicit none
    
    type(soln_t), intent(inout) :: soln
    
    deallocate(soln%U,      &
               soln%V,      &
               soln%F,      &
               soln%D,      &
               soln%M,      &
               soln%T,      &
               soln%R,      &
               soln%dt,     &
               soln%lambda  )
    
  end subroutine deallocate_soln

end module soln_type
