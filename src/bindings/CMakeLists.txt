
add_subdirectory(c)

find_package(SWIG QUIET)
mark_as_advanced(SWIG_EXECUTABLE)

if(SWIG_FOUND)
  option(REPROSIM_BUILD_PYTHON_BINDINGS "Build Python bindings for ${PROJECT_NAME}" YES)


  if(REPROSIM_BUILD_PYTHON_BINDINGS)
    find_package(Python COMPONENTS Interpreter Development.Module NumPy)
    if(Python_FOUND AND Python_Development.Module_FOUND AND Python_NumPy_FOUND)
      add_subdirectory(python)
    else()
      message(WARNING "Python bindings requested but at least one Python requirement was not found: ${Python_FOUND}, ${Python_Development.Module_FOUND}, or ${Python_NumPy_FOUND}.")
    endif()
  endif()
endif()
