ttk_add_base_library(topologicalMapper
  SOURCES
    TopologicalMapper.cpp
  HEADERS
    TopologicalMapper.h
  DEPENDS
    geometry
    unionFind
    Boost::boost
)
target_compile_options(topologicalMapper PUBLIC -ggdb -O0)

