ttk_add_base_library(topoMap
  SOURCES
    TopoMap.cpp
  HEADERS
    TopoMap.h
  DEPENDS
    geometry
    unionFind
    Boost::boost
)

