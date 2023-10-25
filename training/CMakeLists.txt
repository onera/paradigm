if(PDM_ENABLE_TRAINING)
  include(check_python_module)

  if (PDM_ENABLE_FORTRAN)
    check_python_module(sphinxfortran REQUIRED)
  endif()

  find_package(Jupyter REQUIRED)
  find_package(Jupytext REQUIRED)

  set(TRAINING_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/training)
  set(TRAINING_BUILD  ${CMAKE_CURRENT_BINARY_DIR}/training)

  # TODO : Run jupytext --to ipynb on all *.md notebooks

