
add_subdirectory(pFUnit-4.8.0)

# Include the config file so we can use the test creation function directly.
include("${CMAKE_CURRENT_BINARY_DIR}/pFUnit-4.8.0/pFUnitConfig.cmake")


add_subdirectory(utils)

add_subdirectory(lib)
