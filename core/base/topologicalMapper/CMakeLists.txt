ttk_add_base_library(topologicalMapper
  SOURCES
    TopologicalMapper.cpp
  HEADERS
    TopologicalMapper.h
  DEPENDS
    triangulation
    dimensionReduction
    unionFind
)

if(TTK_ENABLE_QHULL)
  target_compile_definitions(topologicalCompression PUBLIC TTK_ENABLE_QHULL)
  target_include_directories(topologicalCompression SYSTEM PRIVATE ${QHULL_INCLUDE_DIR})
  target_link_libraries(topologicalCompression PRIVATE ${QHULL_LIBRARIES})
endif()

