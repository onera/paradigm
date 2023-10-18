if(PDM_ENABLE_NOTEBOOKS)

  message("Coucou!!!")

  # Look for juypter notebooks (pip install jupyter)
  # Look for jupytext (pip install jupytext)
  # TODO : find_package(...)

  set(TRAINING_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/training)
  set(TRAINING_BUILD  ${CMAKE_CURRENT_BINARY_DIR}/training)

  # TODO : Run jupytext --to ipynb on all *.md notebooks

endif()
