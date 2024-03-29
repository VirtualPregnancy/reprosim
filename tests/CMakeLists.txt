

# Unit testing
set(
  tests
  "geometry"
  "diagnostics"
)
set(
  test-srcs
  "testdriver.f90"
)

set(
  TEST_SUITE_USE_STATEMENT
)
set(
  TEST_SUITE_LIST_STATEMENT
)
 
foreach(t IN LISTS tests)
  string(MAKE_C_IDENTIFIER ${t} t) 
  list(APPEND test-srcs "test_${t}.f90")
  set(TEST_SUITE_USE_STATEMENT "${TEST_SUITE_USE_STATEMENT}  use test_${t}, only : collect_${t}\n")
  set(TEST_SUITE_LIST_STATEMENT "${TEST_SUITE_LIST_STATEMENT}    new_testsuite(\"${t}\", collect_${t}), &\n")
endforeach()

string(REGEX REPLACE ", &\n$" " &" TEST_SUITE_LIST_STATEMENT "${TEST_SUITE_LIST_STATEMENT}")
message(STATUS "${TEST_SUITE_LIST_STATEMENT}")
set(test-main "${CMAKE_CURRENT_BINARY_DIR}/testmain.f90")
configure_file(testmain.in.f90 ${test-main})
list(APPEND test-srcs ${test-main})

add_executable(
  "${PROJECT_NAME}-tester"
  "${test-srcs}"
)

target_compile_definitions(
  "${PROJECT_NAME}-tester"
  PRIVATE
  "WITH_QP=$<BOOL:${WITH_QP}>"
  "WITH_XDP=$<BOOL:${WITH_XDP}>"
)

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    set_source_files_properties(testdriver.f90 PROPERTIES COMPILE_FLAGS -cpp)
endif ()

if (MSVC)
  target_compile_options("${PROJECT_NAME}-tester" PUBLIC -fpp)
  set(_REPROSIM_LINK_LIBRARY reprosim_obj)
else()
  set(_REPROSIM_LINK_LIBRARY reprosim)
endif()

set(UTILS_SOURCES data_generator.f90)

add_library(utils ${UTILS_SOURCES})
target_link_libraries(utils PRIVATE reprosim)

target_link_libraries(
  "${PROJECT_NAME}-tester"
  PRIVATE
  ${_REPROSIM_LINK_LIBRARY} 
  utils
)

foreach(t IN LISTS tests)
  add_test("${PROJECT_NAME}/${t}" "${PROJECT_NAME}-tester" "${t}")
endforeach()

