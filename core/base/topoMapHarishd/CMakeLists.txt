ttk_add_base_library(topoMapHarishd
  SOURCES
    TopoMapHarishd.cpp
    geomutils.cpp
  HEADERS
    TopoMapHarishd.h
    DisjointSets.h
    geomutils.h
  DEPENDS
    geometry
    dimensionReduction
    unionFind
)
target_compile_options(topoMapHarishd PUBLIC -ggdb -O0)

