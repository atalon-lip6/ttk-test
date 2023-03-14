#include <TopoMap.h>

ttk::TopoMap::TopoMap() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("TopoMap");
}
ttk::TopoMap::~TopoMap() = default;
