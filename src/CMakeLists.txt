add_subdirectory( constants )
add_subdirectory( functions )
add_subdirectory( types )
add_subdirectory( boundary )
# Propagate lib_srcs list to parent
add_lib_sources( set_precision.f90 )
add_lib_sources( set_inputs.f90 )
add_lib_sources( geometry.f90 )
add_lib_sources( init_problem.f90 )

add_exe_sources( test_program main_program.f90 )
