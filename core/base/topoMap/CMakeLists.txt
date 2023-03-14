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
SET(qhullLibs qhullcpp qhull_r)
  target_compile_definitions(topoMap PUBLIC TTK_ENABLE_QHULL)
  target_include_directories(topoMap PUBLIC ${Qhull_DIR})
  target_link_libraries(topoMap PUBLIC ${qhullLibs})
endif()

