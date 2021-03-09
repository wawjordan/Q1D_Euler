# Propagate lib_srcs list to parent
#add_lib_sources(
#    speed_of_sound.f90
#    cons2prim.f90
#    prim2cons.f90
#    isentropic_relations.f90
#    calc_time_step.f90
#)
add_lib_sources(
     variable_conversion.f90
     time_integration.f90
     fluxes.f90
     other_subroutines.f90
)
