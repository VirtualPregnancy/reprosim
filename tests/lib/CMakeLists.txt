
# Add all the files that make a single test, we could have multiple files testing
# the same module.  Don't add test files into the same test that test different modules.
# These are all .pf files.

set(DIAGNOSTICS_TEST_SRCS
    test_diagnostics.pf)
set(GEOMETRY_TEST_SRCS 
    test_geometry.pf
)
set(PRF_INLET_PRESSURE_SRCS
    test_prf_inlet_pressure.pf
)
   
#test_define_1d_elements.pf 
#test_define_node_geometry.pf
# Make use of the pFUnit helper function to create a test.
# Arguments    : - test_package_name: Name of the test package
#                - test_sources     : List of pf-files to be compiled
#                - extra_sources    : List of extra Fortran source code used for testing (if none, input empty string "")
#                - extra_sources_c  : List of extra C/C++ source code used for testing (if none, input empty string "")
add_pfunit_test(diagnostics_test "${DIAGNOSTICS_TEST_SRCS}" "" "")
add_pfunit_test(geometry_test "${GEOMETRY_TEST_SRCS}" "" "")
add_pfunit_test(prf_inlet_pressure_test "${PRF_INLET_PRESSURE_SRCS}" "" "")
# Link the test to the reprosim library target.
target_link_libraries (diagnostics_test reprosim)
target_link_libraries (geometry_test test_utils)
target_link_libraries (prf_inlet_pressure_test test_utils)
# We shouldn't require this line, but we do WHY???????
target_include_directories(geometry_test PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../utils)
target_include_directories(prf_inlet_pressure_test PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../utils)