add_test_sources( functions test_speed_of_sound.pf )
add_test_sources( functions test_cons2prim.pf )
add_test_sources( functions test_prim2cons.pf )
add_test_sources( functions test_isentropic_relations.pf )
add_test_sources( functions test_calc_time_step.pf )
add_test_sources( functions test_central_flux.pf )

add_test_target( functions unit ${functions_sources} )
