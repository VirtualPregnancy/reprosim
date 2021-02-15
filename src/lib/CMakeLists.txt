
# The list of all source files that are required to build the library
set(LIB_SRCS
  arrays.f90
  diagnostics.f90
  geometry.f90
  indices.f90
  math_utilities.f90
  other_consts.f90
  pressure_resistance_flow.f90
  repro_exports.f90
  solve.f90

)

set(REPROSIM_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Modules")

# Create an object library so we can use the objects in two different libraries
add_library(reprosim_obj OBJECT ${LIB_SRCS})
set_target_properties(reprosim_obj PROPERTIES
  Fortran_MODULE_DIRECTORY "${REPROSIM_MODULE_DIRECTORY}"
  POSITION_INDEPENDENT_CODE TRUE)

if (MSVC)
  # Visual Studio requires a dummy file as it will not build a library from only objects.
  set(dummy_file "${CMAKE_CURRENT_BINARY_DIR}/dummy.f90")
  file(WRITE ${dummy_file} "module dummy
implicit none
contains
end module dummy\n")
endif ()

# Create the reprosim fortran library
add_library(reprosim $<TARGET_OBJECTS:reprosim_obj> ${dummy_file})

set_target_properties(reprosim PROPERTIES
  Fortran_VISIBILITY_PRESET hidden
  C_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN 1)

target_include_directories(reprosim PUBLIC "${REPROSIM_MODULE_DIRECTORY}")

