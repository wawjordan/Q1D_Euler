add_subdirectory( constants )
add_subdirectory( functions )
add_subdirectory( types )
# Propagate lib_srcs list to parent
add_lib_sources( set_precision.f90 )
add_lib_sources( set_inputs.f90 )

add_exe_sources( test_program main_program.f90 )
