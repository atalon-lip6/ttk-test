ttk_add_base_library(topologicalMapper
  SOURCES
    TopologicalMapper.cpp
  HEADERS
    TopologicalMapper.h
  DEPENDS
    triangulation
    dimensionReduction
    unionFind
    distanceMatrixDistorsion
)
target_compile_options(topologicalMapper PUBLIC -ggdb -O0)
if(TTK_ENABLE_QHULL)
SET(qhullLibs qhullcpp qhull_r)
  target_compile_definitions(topologicalMapper PRIVATE TTK_ENABLE_QHULL)
  target_include_directories(topologicalMapper PRIVATE ${Qhull_DIR})
  target_link_libraries(topologicalMapper PRIVATE ${qhullLibs})
endif()

