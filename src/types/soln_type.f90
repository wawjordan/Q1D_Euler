module soln_type

  use set_precision, only : prec
  use set_constants, only : zero, one
  use set_inputs,    only : neq, max_iter, i_low, i_high, ig_low, ig_high
  
  implicit none

  type soln_t
    
    real(prec), allocatable, dimension(:,:) :: V
    real(prec), allocatable, dimension(:,:) :: U
    real(prec), allocatable, dimension(:,:) :: R
    real(prec), allocatable, dimension(:,:) :: F
    real(prec), allocatable, dimension(:,:) :: D
    real(prec), allocatable, dimension(:)   :: asnd
    real(prec), allocatable, dimension(:)   :: mach
    real(prec), allocatable, dimension(:)   :: temp
    real(prec), allocatable, dimension(:)   :: dt
    real(prec), allocatable, dimension(:)   :: src
    real(prec), allocatable, dimension(:)   :: lambda
    real(prec), allocatable, dimension(:,:) :: rnorm 
    real(prec), allocatable, dimension(:,:) :: rinit

  end type soln_t

  contains
  
  
  !============================= allocate_soln ===============================80
  !>
  !! Description: 
  !!
  !! Inputs:      soln : 
  !!
  !! Outputs:     soln : 
  !<
  !===========================================================================80
  subroutine allocate_soln( soln )

    type(soln_t), intent(inout) :: soln
    
    allocate( soln%V( ig_low:ig_high, neq ), &
              soln%U( ig_low:ig_high, neq ), &
              soln%R( i_low:i_high,   neq ), &
              soln%F( i_low-1:i_high, neq ), &
              soln%D( i_low-1:i_high, neq ), &
              soln%asnd( ig_low:ig_high ),      &
              soln%mach( ig_low:ig_high ),      &
              soln%temp( ig_low:ig_high ),      &
              soln%Src( ig_low:ig_high ),      &
              soln%dt( ig_low:ig_high ),     &
              soln%lambda( ig_low:ig_high ), &
               soln%rinit( 1, neq ), &
              soln%rnorm( 1:max_iter, neq ) )

    soln%V   = zero
    soln%U   = zero
    soln%R   = zero
    soln%F   = zero
    soln%D   = zero
    soln%asnd   = zero
    soln%mach   = zero
    soln%temp   = zero
    soln%src   = zero
    soln%dt  = zero
    soln%lambda = zero
    soln%rinit = zero
    soln%rnorm = zero

  end subroutine allocate_soln
  
  
  !=========================== deallocate_soln ===============================80
  !>
  !! Description: 
  !!
  !! Inputs:      soln : 
  !!
  !! Outputs:     soln : 
  !<
  !===========================================================================80
  subroutine deallocate_soln( soln )
  
    implicit none
    
    type(soln_t), intent(inout) :: soln
    
    deallocate(soln%V,      &
               soln%U,      &
               soln%R,      &
               soln%F,      &
               soln%D,      &
               soln%asnd,      &
               soln%mach,      &
               soln%temp,      &
               soln%src,      &
               soln%dt,     &
               soln%lambda, &
               soln%rinit,  &
               soln%rnorm   )
    
  end subroutine deallocate_soln

end module soln_type
