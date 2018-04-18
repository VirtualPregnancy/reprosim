set(TEST_UTILS_SOURCES test_data.f90)

add_library(test_utils ${TEST_UTILS_SOURCES})

target_link_libraries(test_utils reprosim)
