# register a new filter to build in the TTK plugin
# deduce the location of the corresonding vtk.module file
# also register the xml file if given
macro(ttk_register_pv_filter vtkModuleDir xmlFile)
  if(NOT EXISTS "${VTKWRAPPER_DIR}/${vtkModuleDir}/vtk.module")
    message(FATAL_ERROR
      "Register a paraview module without the corresponding vtk.module: "
      ${VTKWRAPPER_DIR}/${vtkModuleDir}
      )
  endif()

  ttk_parse_module_file("${VTKWRAPPER_DIR}/${vtkModuleDir}/ttk.module")
  cmake_parse_arguments("TTK" "" "NAME" "SOURCES;HEADERS;DEPENDS" ${moduleFileContent})

  ttk_get_target(${TTK_NAME} TTK_TARGET)
  message(STATUS "Target: ${TTK_TARGET}")
  if(NOT "${TTK_TARGET}" STREQUAL "")
    message(STATUS "Append: ${vtkModuleDir}")
    list(APPEND TTK_MODULES ${vtkModuleDir})
    list(APPEND TTK_VTK_MODULE_FILES ${VTKWRAPPER_DIR}/${vtkModuleDir}/vtk.module)
    if(${xmlFile})
      list(APPEND TTK_XMLS ${CMAKE_CURRENT_LIST_DIR}/${xmlFile})
    endif()
  endif()
endmacro()
