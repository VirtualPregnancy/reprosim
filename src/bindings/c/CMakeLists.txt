
project(ReprosimCBindings LANGUAGES C Fortran)

set(C_FORTRAN_LIB_SRCS
    arrays.f90
    diagnostics.f90
    exports.f90
    geometry.f90
    indices.f90
    nutrient_transport.f90
    pressure_resistance_flow.f90
    utils.f90
)
set(C_C_LIB_SRCS
    arrays.c
    diagnostics.c
    exports.c
    geometry.c
    indices.c
    nutrient_transport.c
    pressure_resistance_flow.c
    utils.c
)

set(C_LIB_HDRS
  symbol_export.h
  arrays.h
  diagnostics.h
  exports.h
  geometry.h
  indices.h
  nutrient_transport.h
  pressure_resistance_flow.h
  utils.h

)

add_library(reprosim_c_fortran STATIC ${C_FORTRAN_LIB_SRCS} $<TARGET_OBJECTS:reprosim_obj>)
set_target_properties(reprosim_c_fortran PROPERTIES
    POSITION_INDEPENDENT_CODE True
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Modules"
    FOLDER bindings/c)
target_include_directories(reprosim_c_fortran PRIVATE $<TARGET_PROPERTY:reprosim_obj,Fortran_MODULE_DIRECTORY>)
add_library(reprosim_c SHARED ${C_C_LIB_SRCS} ${C_LIB_HDRS}) #  reprosim to be replaced with whatever the final name of the code is
set_target_properties(reprosim_c PROPERTIES
    FOLDER bindings/c
    C_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1)
# Have to set this manually, as it is only set automatically if the file extension has an
# uppercase 'F' and not lowercase as we have here for GNU compiler.
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    set_source_files_properties(${C_FORTRAN_LIB_SRCS} PROPERTIES COMPILE_FLAGS -cpp)
endif ()
target_link_libraries(reprosim_c PUBLIC reprosim_c_fortran)

add_custom_target(cbindings ALL)
add_dependencies(cbindings reprosim_c)

