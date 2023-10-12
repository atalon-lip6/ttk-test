ttk_add_base_library(topoMap
  SOURCES
    TopoMap.cpp
  HEADERS
    TopoMap.h
  DEPENDS
    geometry
    unionFind
)

if(TTK_ENABLE_QHULL)
  target_compile_definitions(topoMap PUBLIC TTK_ENABLE_QHULL)
  target_link_libraries(topoMap PRIVATE Qhull::qhullcpp)
  target_link_libraries(topoMap PRIVATE Qhull::qhull_r)
endif()

