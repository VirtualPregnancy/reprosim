
add_subdirectory(pFUnit)

# Include the config file so we can use the test creation function directly.
include("${CMAKE_CURRENT_BINARY_DIR}/pFUnit/pFUnitConfig.cmake")


add_subdirectory(utils)

add_subdirectory(lib)
