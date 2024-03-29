
project(ReprosimPythonBindings LANGUAGES C)


set(ADDITIONAL_GENERATED_FILES)
# SWIG has already been found at this point.
include(${SWIG_USE_FILE})

set(PYTHONLIBS_SWIG_FLAGS -relativeimport)

set(INTERFACE_SRCS
  ../interface/arrays.i
  ../interface/diagnostics.i
  ../interface/fetal.i
  ../interface/geometry.i
  ../interface/indices.i
  ../interface/pressure_resistance_flow.i
  ../interface/repro_exports.i

)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../c)

if(MSVC)
  set(PYTHON_PACKAGE_DIR ${CMAKE_CURRENT_BINARY_DIR}/$(Configuration)/reprosim)
  set(PYTHON_PACKAGE_INIT_FILE ${CMAKE_CURRENT_BINARY_DIR}/$(Configuration)/reprosim/__init__.py)
  set(PYTHON_SETUP_PY_FILE_DIR ${CMAKE_CURRENT_BINARY_DIR}/$(Configuration))
else()
  set(PYTHON_PACKAGE_DIR ${CMAKE_CURRENT_BINARY_DIR}/reprosim)
  set(PYTHON_PACKAGE_INIT_FILE ${CMAKE_CURRENT_BINARY_DIR}/reprosim/__init__.py)
  set(PYTHON_SETUP_PY_FILE_DIR ${CMAKE_CURRENT_BINARY_DIR})
endif()
list(APPEND ADDITIONAL_GENERATED_FILES ${PYTHON_PACKAGE_INIT_FILE})
set(CMAKE_SWIG_OUTDIR "${PYTHON_PACKAGE_DIR}")

set(CMAKE_SWIG_FLAGS ${PYTHONLIBS_SWIG_FLAGS})

set(MODULE_TARGETS)
foreach(SWIG_INTERFACE ${INTERFACE_SRCS})
  get_filename_component(MODULE_NAME ${SWIG_INTERFACE} NAME_WE)
  swig_add_library(${MODULE_NAME} LANGUAGE python SOURCES ${SWIG_INTERFACE})
  
  # The Python debug libraries on windows also require the Python release libraries
  # (this appears to be the case for the libraries built with Visual Studio 14 2015),
  # so we add the release libraries here for all platforms.  This probably means that
  # when we try and link to the debug version of the Python libraires on OSX and GNU/Linux 
  # we will not succeed.  But as this is rarely done it hopefully won't become an issue.
  swig_link_libraries(${MODULE_NAME} ${PYTHON_LIBRARY_RELEASE} reprosim_c Python::Module Python::NumPy)
  set(MODULE_TARGET ${SWIG_MODULE_${MODULE_NAME}_REAL_NAME})
  set_target_properties(${MODULE_TARGET} PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${PYTHON_PACKAGE_DIR}
    RUNTIME_OUTPUT_DIRECTORY ${PYTHON_PACKAGE_DIR}
    FOLDER bindings/python
    C_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1)
  list(APPEND MODULE_TARGETS ${MODULE_TARGET})

  if (MSVC)
    if (${PYTHON_LIBRARY_DEBUG})
      # Again from the above comment we need the debug Python library and release Python library when
      # linking with link.exe on Windows.
      target_link_libraries(${MODULE_TARGET} $<$<CONFIG:Debug>:PYTHON_LIBRARY_DEBUG>)
    else()
      # Disable use of Debug Python libraries when not available and using MSVC.
      # See http://stackoverflow.com/questions/11311877/creating-a-dll-from-a-wrapped-cpp-file-with-swig
      target_compile_definitions(${MODULE_TARGET} PRIVATE $<$<CONFIG:Debug>:SWIG_PYTHON_INTERPRETER_NO_DEBUG>)
    endif()
    set_target_properties(${MODULE_TARGET} PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY_DEBUG ${PYTHON_PACKAGE_DIR}
      LIBRARY_OUTPUT_DIRECTORY_RELEASE ${PYTHON_PACKAGE_DIR})
  endif()
endforeach()

# Alternative option is to configure the setup.py file to the appropriate location.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py ${PYTHON_SETUP_PY_FILE_DIR}/setup.py @ONLY)

add_custom_target(pybindings ALL)
add_dependencies(pybindings ${MODULE_TARGETS})
# Create the Python __init__.py file.
add_custom_command( TARGET pybindings POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E touch ${PYTHON_PACKAGE_INIT_FILE})

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${ADDITIONAL_GENERATED_FILES}")
