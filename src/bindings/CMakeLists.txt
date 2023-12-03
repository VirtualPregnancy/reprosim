
add_subdirectory(c)

find_package(SWIG QUIET)
mark_as_advanced(SWIG_EXECUTABLE)

if(SWIG_FOUND)
  option(REPROSIM_BUILD_PYTHON_BINDINGS "Build Python bindings for ${PROJECT_NAME}" YES)


  if(REPROSIM_BUILD_PYTHON_BINDINGS)
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
      set(_FIND_PYTHON_DEVELOPMENT_TYPE Development.Module)
    else()
      set(_FIND_PYTHON_DEVELOPMENT_TYPE Development)
    endif()
    find_package(Python COMPONENTS Interpreter ${_FIND_PYTHON_DEVELOPMENT_TYPE})
    if(Python_FOUND AND Python_Development_FOUND)
      add_subdirectory(python)
    else()
      message(WARNING "Python bindings requested but a Python interpreter and library were not found.")
    endif()
  endif()
endif()
