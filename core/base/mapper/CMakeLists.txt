ttk_add_base_library(mapper
  SOURCES
    Mapper.cpp
  HEADERS
    Mapper.h
  DEPENDS
    triangulation
    unionFind
    dimensionReduction
    lDistanceMatrix
  )

target_compile_options(mapper PUBLIC -ggdb -O0)
#target_compile_options(mapper PUBLIC -fsanitize=address)
#target_link_options(mapper PUBLIC -fsanitize=address)
