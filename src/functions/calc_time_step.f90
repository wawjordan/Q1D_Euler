subroutine calc_time_step( lambda, V )
  
  use set_precision, only : prec
  use fluid_constants, only : gamma
  use set_inputs, only : CFL
