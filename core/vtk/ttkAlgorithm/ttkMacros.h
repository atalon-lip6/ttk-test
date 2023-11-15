#pragma once

class vtkIdTypeArray;
class vtkIntArray;

#define TTK_COMMA ,

#ifdef TTK_ENABLE_64BIT_IDS
using ttkSimplexIdTypeArray = vtkIdTypeArray;
#else
using ttkSimplexIdTypeArray = vtkIntArray;
#endif

#ifndef vtkSetEnumMacro
#define vtkSetEnumMacro(name, enumType)                                        \
  virtual void Set##name(enumType _arg) {                                      \
    vtkDebugMacro(<< this->GetClassName() << " (" << this                      \
                  << "): setting " #name " to "                                \
                  << static_cast<std::underlying_type<enumType>::type>(_arg)); \
    if(this->name != _arg) {                                                   \
      this->name = _arg;                                                       \
      this->Modified();                                                        \
    }                                                                          \
  }
#endif

#ifndef vtkGetEnumMacro
#define vtkGetEnumMacro(name, enumType)                                      \
  virtual enumType Get##name() const {                                       \
    vtkDebugMacro(<< this->GetClassName() << " (" << this << "): returning " \
                  << #name " of "                                            \
                  << static_cast<std::underlying_type<enumType>::type>(      \
                       this->name));                                         \
    return this->name;                                                       \
  }
#endif

#define ttkSetEnumMacro(name, enumType)                   \
  virtual void Set##name(int _arg) {                      \
    vtkDebugMacro(<< this->GetClassName() << " (" << this \
                  << "): setting " #name " to " << _arg); \
    if(this->name != static_cast<enumType>(_arg)) {       \
      this->name = static_cast<enumType>(_arg);           \
      this->Modified();                                   \
    }                                                     \
  }                                                       \
  vtkSetEnumMacro(name, enumType);

#ifdef TTK_REDUCE_TEMPLATE_INSTANTIATIONS
// reduced list of template instantiations by redefining vtkTemplateMacro
#include <vtkSetGet.h>
#ifdef vtkTemplateMacro
#undef vtkTemplateMacro
#define vtkTemplateMacro(call)                    \
  vtkTemplateMacroCase(VTK_DOUBLE, double, call); \
  vtkTemplateMacroCase(VTK_FLOAT, float, call);   \
  vtkTemplateMacroCase(VTK_INT, int, call);       \
  vtkTemplateMacroCase(VTK_LONG_LONG, long long, call);
#endif // vtkTemplateMacro
#endif // TTK_REDUCE_TEMPLATE_INSTANTIATIONS

#define ttkVtkTemplateMacroCase(                         \
  dataType, triangulationType, triangulationClass, card, call) \
case triangulationType: {                                               \
    switch(card) { \
      case 0: \
        {          \
        using TTK_TT = triangulationClass<0>;                                 \
        switch(dataType) { vtkTemplateMacro((call)); };      \
        } \
        break; \
      case 1: \
        {          \
        using TTK_TT = triangulationClass<1>;                                 \
        switch(dataType) { vtkTemplateMacro((call)); };      \
        } \
        break; \
      case 2: \
        {          \
        using TTK_TT = triangulationClass<2>;                                 \
        switch(dataType) { vtkTemplateMacro((call)); };      \
        } \
        break; \
      case 3: \
        {          \
        using TTK_TT = triangulationClass<3>;                                 \
        switch(dataType) { vtkTemplateMacro((call)); };      \
        } \
        break; \
    } \
  }; break



#define ttkVtkTemplateMacro(dataType, triangulationType, card, call)            \
  switch(triangulationType) {                                             \
    ttkVtkTemplateMacroCase(dataType, ttk::Triangulation::Type::EXPLICIT, \
                            ttk::ExplicitTriangulation, card, call);            \
    ttkVtkTemplateMacroCase(dataType, ttk::Triangulation::Type::IMPLICIT, \
                            ttk::ImplicitNoPreconditions, card, call);          \
    ttkVtkTemplateMacroCase(dataType,                                     \
                            ttk::Triangulation::Type::HYBRID_IMPLICIT,    \
                            ttk::ImplicitWithPreconditions, card, call);        \
    ttkVtkTemplateMacroCase(dataType, ttk::Triangulation::Type::PERIODIC, \
                            ttk::PeriodicNoPreconditions, card, call);          \
    ttkVtkTemplateMacroCase(dataType,                                     \
                            ttk::Triangulation::Type::HYBRID_PERIODIC,    \
                            ttk::PeriodicWithPreconditions, card, call);        \
    ttkVtkTemplateMacroCase(dataType, ttk::Triangulation::Type::COMPACT,  \
                            ttk::CompactTriangulation, card, call);             \
  }

#define ttkTemplate2IdMacro(call)                                           \
  vtkTemplate2MacroCase1(VTK_LONG_LONG, long long, call);                   \
  vtkTemplate2MacroCase1(VTK_UNSIGNED_LONG_LONG, unsigned long long, call); \
  vtkTemplate2MacroCase1(VTK_ID_TYPE, vtkIdType, call);                     \
  vtkTemplate2MacroCase1(VTK_LONG, long, call);                             \
  vtkTemplate2MacroCase1(VTK_UNSIGNED_LONG, unsigned long, call);           \
  vtkTemplate2MacroCase1(VTK_INT, int, call);                               \
  vtkTemplate2MacroCase1(VTK_UNSIGNED_INT, unsigned int, call);

#ifndef vtkTemplate2MacroCase1
#define vtkTemplate2MacroCase1(type1N, type1, call)                            \
  vtkTemplate2MacroCase2(type1N, type1, VTK_DOUBLE, double, call);             \
  vtkTemplate2MacroCase2(type1N, type1, VTK_FLOAT, float, call);               \
  vtkTemplate2MacroCase2(type1N, type1, VTK_LONG_LONG, long long, call);       \
  vtkTemplate2MacroCase2(                                                      \
    type1N, type1, VTK_UNSIGNED_LONG_LONG, unsigned long long, call);          \
  vtkTemplate2MacroCase2(type1N, type1, VTK_ID_TYPE, vtkIdType, call);         \
  vtkTemplate2MacroCase2(type1N, type1, VTK_LONG, long, call);                 \
  vtkTemplate2MacroCase2(                                                      \
    type1N, type1, VTK_UNSIGNED_LONG, unsigned long, call);                    \
  vtkTemplate2MacroCase2(type1N, type1, VTK_INT, int, call);                   \
  vtkTemplate2MacroCase2(type1N, type1, VTK_UNSIGNED_INT, unsigned int, call); \
  vtkTemplate2MacroCase2(type1N, type1, VTK_SHORT, short, call);               \
  vtkTemplate2MacroCase2(                                                      \
    type1N, type1, VTK_UNSIGNED_SHORT, unsigned short, call);                  \
  vtkTemplate2MacroCase2(type1N, type1, VTK_CHAR, char, call);                 \
  vtkTemplate2MacroCase2(type1N, type1, VTK_SIGNED_CHAR, signed char, call);   \
  vtkTemplate2MacroCase2(type1N, type1, VTK_UNSIGNED_CHAR, unsigned char, call)
#endif

#ifndef vtkTemplate2MacroCase2
#define vtkTemplate2MacroCase2(type1N, type1, type2N, type2, call) \
  case vtkTemplate2PackMacro(type1N, type2N): {                    \
    typedef type1 VTK_T1;                                          \
    typedef type2 VTK_T2;                                          \
    call;                                                          \
  }; break
#endif

// -----------------------------------------------------------------------------

#define ttkTypeMacroErrorCase(idx, type)                          \
  default: {                                                      \
    this->printErr("Unsupported " #idx "-th Template Data Type: " \
                   + std::to_string(static_cast<int>(type)));     \
  } break;

#define ttkTypeMacroCase(enum, type, card, number, call) \
  case enum: {                                     \
      switch(card) { \
      case 0: \
        {          \
        using T##number = type<0>;                                 \
        call;                                          \
        } \
        break; \
      case 1: \
        {          \
        using T##number = type<1>;                                 \
        call;                                          \
        } \
        break; \
      case 2: \
        {          \
        using T##number = type<2>;                                 \
        call;                                          \
        } \
        break; \
      case 3: \
        {          \
        using T##number = type<3>;                                 \
        call;                                          \
        } \
        break; \
    } \
  }; break



#define ttkTypeMacroT(group, card, call)                                            \
  switch(group) {                                                             \
    ttkTypeMacroCase(ttk::Triangulation::Type::EXPLICIT,                      \
                     ttk::ExplicitTriangulation, card, 0, call);                    \
    ttkTypeMacroCase(                                                         \
      ttk::Triangulation::Type::COMPACT, ttk::CompactTriangulation, card, 0, call); \
    ttkTypeMacroCase(ttk::Triangulation::Type::IMPLICIT,                      \
                     ttk::ImplicitNoPreconditions, card, 0, call);                  \
    ttkTypeMacroCase(ttk::Triangulation::Type::HYBRID_IMPLICIT,               \
                     ttk::ImplicitWithPreconditions, card, 0, call);                \
    ttkTypeMacroCase(ttk::Triangulation::Type::PERIODIC,                      \
                     ttk::PeriodicNoPreconditions, card, 0, call);                  \
    ttkTypeMacroCase(ttk::Triangulation::Type::HYBRID_PERIODIC,               \
                     ttk::PeriodicWithPreconditions, card, 0, call);                \
    ttkTypeMacroErrorCase(0, group);                                          \
  }

#define ttkTypeMacroCaseSimple(enum, type, number, call) \
  case enum: {                                     \
        using T##number = type;                                 \
        call;                                          \
  }; break



#define ttkTypeMacroR(group, call)                 \
  switch(group) {                                  \
    ttkTypeMacroCaseSimple(VTK_FLOAT, float, 0, call);   \
    ttkTypeMacroCaseSimple(VTK_DOUBLE, double, 0, call); \
    ttkTypeMacroErrorCase(0, group);               \
  }

#define ttkTypeMacroI(group, call)                                         \
  switch(group) {                                                          \
    ttkTypeMacroCaseSimple(VTK_INT, int, 0, call);                               \
    ttkTypeMacroCaseSimple(VTK_UNSIGNED_INT, unsigned int, 0, call);             \
    ttkTypeMacroCaseSimple(VTK_CHAR, char, 0, call);                             \
    ttkTypeMacroCaseSimple(VTK_SIGNED_CHAR, signed char, 0, call);               \
    ttkTypeMacroCaseSimple(VTK_UNSIGNED_CHAR, unsigned char, 0, call);           \
    ttkTypeMacroCaseSimple(VTK_LONG, long, 0, call);                             \
    ttkTypeMacroCaseSimple(VTK_LONG_LONG, long long, 0, call);                   \
    ttkTypeMacroCaseSimple(VTK_UNSIGNED_LONG, unsigned long, 0, call);           \
    ttkTypeMacroCaseSimple(VTK_UNSIGNED_LONG_LONG, unsigned long long, 0, call); \
    ttkTypeMacroCaseSimple(VTK_ID_TYPE, vtkIdType, 0, call);                     \
    ttkTypeMacroErrorCase(0, group);                                       \
  }

#define ttkTypeMacroA(group, call)                                         \
  switch(group) {                                                          \
    ttkTypeMacroCaseSimple(VTK_FLOAT, float, 0, call);                           \
    ttkTypeMacroCaseSimple(VTK_DOUBLE, double, 0, call);                         \
    ttkTypeMacroCaseSimple(VTK_INT, int, 0, call);                               \
    ttkTypeMacroCaseSimple(VTK_UNSIGNED_INT, unsigned int, 0, call);             \
    ttkTypeMacroCaseSimple(VTK_CHAR, char, 0, call);                             \
    ttkTypeMacroCaseSimple(VTK_SIGNED_CHAR, signed char, 0, call);               \
    ttkTypeMacroCaseSimple(VTK_UNSIGNED_CHAR, unsigned char, 0, call);           \
    ttkTypeMacroCaseSimple(VTK_LONG, long, 0, call);                             \
    ttkTypeMacroCaseSimple(VTK_LONG_LONG, long long, 0, call);                   \
    ttkTypeMacroCaseSimple(VTK_UNSIGNED_LONG, unsigned long, 0, call);           \
    ttkTypeMacroCaseSimple(VTK_UNSIGNED_LONG_LONG, unsigned long long, 0, call); \
    ttkTypeMacroCaseSimple(VTK_ID_TYPE, vtkIdType, 0, call);                     \
    ttkTypeMacroErrorCase(0, group);                                       \
  }

#ifdef TTK_REDUCE_TEMPLATE_INSTANTIATIONS
// reduced list of template instantiations by redefining ttkTypeMacroI
// & ttkTypeMacroA
#undef ttkTypeMacroI
#define ttkTypeMacroI(group, call)                       \
  switch(group) {                                        \
    ttkTypeMacroCaseSimple(VTK_INT, int, 0, call);             \
    ttkTypeMacroCaseSimple(VTK_LONG_LONG, long long, 0, call); \
    ttkTypeMacroErrorCase(0, group);                     \
  }
#undef ttkTypeMacroA
#define ttkTypeMacroA(group, call)                       \
  switch(group) {                                        \
    ttkTypeMacroCaseSimple(VTK_FLOAT, float, 0, call);         \
    ttkTypeMacroCaseSimple(VTK_DOUBLE, double, 0, call);       \
    ttkTypeMacroCaseSimple(VTK_INT, int, 0, call);             \
    ttkTypeMacroCaseSimple(VTK_LONG_LONG, long long, 0, call); \
    ttkTypeMacroErrorCase(0, group);                     \
  }
#endif // TTK_REDUCE_TEMPLATE_INSTANTIATIONS

#define ttkTypeMacroAT(group0, group1, card, call)                    \
  switch(group1) {                                              \
    ttkTypeMacroCase(ttk::Triangulation::Type::EXPLICIT,        \
                     ttk::ExplicitTriangulation, card, 1,             \
                     ttkTypeMacroA(group0, call));              \
    ttkTypeMacroCase(ttk::Triangulation::Type::COMPACT,         \
                     ttk::CompactTriangulation, card, 1,              \
                     ttkTypeMacroA(group0, call));              \
    ttkTypeMacroCase(ttk::Triangulation::Type::IMPLICIT,        \
                     ttk::ImplicitNoPreconditions, card, 1,           \
                     ttkTypeMacroA(group0, call));              \
    ttkTypeMacroCase(ttk::Triangulation::Type::HYBRID_IMPLICIT, \
                     ttk::ImplicitWithPreconditions, card, 1,         \
                     ttkTypeMacroA(group0, call));              \
    ttkTypeMacroCase(ttk::Triangulation::Type::PERIODIC,        \
                     ttk::PeriodicNoPreconditions, card, 1,           \
                     ttkTypeMacroA(group0, call));              \
    ttkTypeMacroCase(ttk::Triangulation::Type::HYBRID_PERIODIC, \
                     ttk::PeriodicWithPreconditions, card, 1,         \
                     ttkTypeMacroA(group0, call));              \
    ttkTypeMacroErrorCase(1, group1);                           \
  }

#define ttkTypeMacroRT(group0, group1, card, call)                    \
  switch(group1) {                                              \
    ttkTypeMacroCase(ttk::Triangulation::Type::EXPLICIT,        \
                     ttk::ExplicitTriangulation, card, 1,             \
                     ttkTypeMacroR(group0, call));              \
    ttkTypeMacroCase(ttk::Triangulation::Type::COMPACT,         \
                     ttk::CompactTriangulation, card, 1,              \
                     ttkTypeMacroR(group0, call));              \
    ttkTypeMacroCase(ttk::Triangulation::Type::IMPLICIT,        \
                     ttk::ImplicitNoPreconditions, card, 1,           \
                     ttkTypeMacroR(group0, call));              \
    ttkTypeMacroCase(ttk::Triangulation::Type::HYBRID_IMPLICIT, \
                     ttk::ImplicitWithPreconditions, card, 1,         \
                     ttkTypeMacroR(group0, call));              \
    ttkTypeMacroCase(ttk::Triangulation::Type::PERIODIC,        \
                     ttk::PeriodicNoPreconditions, card, 1,           \
                     ttkTypeMacroR(group0, call));              \
    ttkTypeMacroCase(ttk::Triangulation::Type::HYBRID_PERIODIC, \
                     ttk::PeriodicWithPreconditions, card, 1,         \
                     ttkTypeMacroR(group0, call));              \
    ttkTypeMacroErrorCase(1, group1);                           \
  }

#define ttkTypeMacroIT(group0, group1, card, call)                    \
  switch(group1) {                                              \
    ttkTypeMacroCase(ttk::Triangulation::Type::EXPLICIT,        \
                     ttk::ExplicitTriangulation, card, 1,             \
                     ttkTypeMacroI(group0, call));              \
    ttkTypeMacroCase(ttk::Triangulation::Type::COMPACT,         \
                     ttk::CompactTriangulation, card, 1,              \
                     ttkTypeMacroI(group0, call));              \
    ttkTypeMacroCase(ttk::Triangulation::Type::IMPLICIT,        \
                     ttk::ImplicitNoPreconditions, card, 1,           \
                     ttkTypeMacroI(group0, call));              \
    ttkTypeMacroCase(ttk::Triangulation::Type::HYBRID_IMPLICIT, \
                     ttk::ImplicitWithPreconditions, card, 1,         \
                     ttkTypeMacroI(group0, call));              \
    ttkTypeMacroCase(ttk::Triangulation::Type::PERIODIC,        \
                     ttk::PeriodicNoPreconditions, card, 1,           \
                     ttkTypeMacroI(group0, call));              \
    ttkTypeMacroCase(ttk::Triangulation::Type::HYBRID_PERIODIC, \
                     ttk::PeriodicWithPreconditions, card, 1,         \
                     ttkTypeMacroI(group0, call));              \
    ttkTypeMacroErrorCase(1, group1);                           \
  }



#define ttkTypeMacroAI(group0, group1, call)                                  \
  switch(group1) {                                                            \
    ttkTypeMacroCaseSimple(VTK_INT, int, 1, ttkTypeMacroA(group0, call));           \
    ttkTypeMacroCaseSimple(                                                         \
      VTK_LONG_LONG, long long, 1, ttkTypeMacroA(group0, call));              \
    ttkTypeMacroCaseSimple(VTK_ID_TYPE, vtkIdType, 1, ttkTypeMacroA(group0, call)); \
    ttkTypeMacroErrorCase(1, group1);                                         \
  }

#define ttkTypeMacroRR(group0, group1, call)                              \
  switch(group1) {                                                        \
    ttkTypeMacroCaseSimple(VTK_FLOAT, float, 1, ttkTypeMacroR(group0, call));   \
    ttkTypeMacroCaseSimple(VTK_DOUBLE, double, 1, ttkTypeMacroR(group0, call)); \
    ttkTypeMacroErrorCase(1, group1);                                     \
  }

#define ttkTypeMacroAA(group0, group1, call)                                  \
  switch(group1) {                                                            \
    ttkTypeMacroCaseSimple(VTK_FLOAT, float, 1, ttkTypeMacroA(group0, call));       \
    ttkTypeMacroCaseSimple(VTK_DOUBLE, double, 1, ttkTypeMacroA(group0, call));     \
    ttkTypeMacroCaseSimple(VTK_INT, int, 1, ttkTypeMacroA(group0, call));           \
    ttkTypeMacroCaseSimple(                                                         \
      VTK_UNSIGNED_INT, unsigned int, 1, ttkTypeMacroA(group0, call));        \
    ttkTypeMacroCaseSimple(VTK_CHAR, char, 1, ttkTypeMacroA(group0, call));         \
    ttkTypeMacroCaseSimple(                                                         \
      VTK_SIGNED_CHAR, signed char, 1, ttkTypeMacroA(group0, call));          \
    ttkTypeMacroCaseSimple(                                                         \
      VTK_UNSIGNED_CHAR, unsigned char, 1, ttkTypeMacroA(group0, call));      \
    ttkTypeMacroCaseSimple(VTK_LONG, long, 1, ttkTypeMacroA(group0, call));         \
    ttkTypeMacroCaseSimple(                                                         \
      VTK_LONG_LONG, long long, 1, ttkTypeMacroA(group0, call));              \
    ttkTypeMacroCaseSimple(                                                         \
      VTK_UNSIGNED_LONG, unsigned long, 1, ttkTypeMacroA(group0, call));      \
    ttkTypeMacroCaseSimple(VTK_UNSIGNED_LONG_LONG, unsigned long long, 1,           \
                     ttkTypeMacroA(group0, call));                            \
    ttkTypeMacroCaseSimple(VTK_ID_TYPE, vtkIdType, 1, ttkTypeMacroA(group0, call)); \
    ttkTypeMacroErrorCase(1, group1);                                         \
  }

#define ttkTypeMacroAAA(group0, group1, group2, call)                          \
  switch(group2) {                                                             \
    ttkTypeMacroCaseSimple(                                                          \
      VTK_FLOAT, float, 2, ttkTypeMacroAA(group0, group1, call));              \
    ttkTypeMacroCaseSimple(                                                          \
      VTK_DOUBLE, double, 2, ttkTypeMacroAA(group0, group1, call));            \
    ttkTypeMacroCaseSimple(VTK_INT, int, 2, ttkTypeMacroAA(group0, group1, call));   \
    ttkTypeMacroCaseSimple(VTK_UNSIGNED_INT, unsigned int, 2,                        \
                     ttkTypeMacroAA(group0, group1, call));                    \
    ttkTypeMacroCaseSimple(VTK_CHAR, char, 2, ttkTypeMacroAA(group0, group1, call)); \
    ttkTypeMacroCaseSimple(                                                          \
      VTK_SIGNED_CHAR, signed char, 2, ttkTypeMacroAA(group0, group1, call));  \
    ttkTypeMacroCaseSimple(VTK_UNSIGNED_CHAR, unsigned char, 2,                      \
                     ttkTypeMacroAA(group0, group1, call));                    \
    ttkTypeMacroCaseSimple(VTK_LONG, long, 2, ttkTypeMacroAA(group0, group1, call)); \
    ttkTypeMacroCaseSimple(                                                          \
      VTK_LONG_LONG, long long, 2, ttkTypeMacroAA(group0, group1, call));      \
    ttkTypeMacroCaseSimple(VTK_UNSIGNED_LONG, unsigned long, 2,                      \
                     ttkTypeMacroAA(group0, group1, call));                    \
    ttkTypeMacroCaseSimple(VTK_UNSIGNED_LONG_LONG, unsigned long long, 2,            \
                     ttkTypeMacroAA(group0, group1, call));                    \
    ttkTypeMacroCaseSimple(                                                          \
      VTK_ID_TYPE, vtkIdType, 2, ttkTypeMacroAA(group0, group1, call));        \
    ttkTypeMacroErrorCase(2, group2);                                          \
  }

#define ttkTypeMacroAII(group0, group1, group2, call)                        \
  switch(group2) {                                                           \
    ttkTypeMacroCaseSimple(VTK_INT, int, 2, ttkTypeMacroAI(group0, group1, call)); \
    ttkTypeMacroCaseSimple(                                                        \
      VTK_LONG_LONG, long long, 2, ttkTypeMacroAI(group0, group1, call));    \
    ttkTypeMacroCaseSimple(                                                        \
      VTK_ID_TYPE, vtkIdType, 2, ttkTypeMacroAI(group0, group1, call));      \
    ttkTypeMacroErrorCase(2, group2);                                        \
  }

#define ttkTypeMacroRRR(group0, group1, group2, call)               \
  switch(group2) {                                                  \
    ttkTypeMacroCaseSimple(                                               \
      VTK_FLOAT, float, 2, ttkTypeMacroRR(group0, group1, call));   \
    ttkTypeMacroCaseSimple(                                               \
      VTK_DOUBLE, double, 2, ttkTypeMacroRR(group0, group1, call)); \
    ttkTypeMacroErrorCase(2, group2);                               \
  }

#define ttkTypeMacroRRI(group0, group1, group2, call)                        \
  switch(group2) {                                                           \
    ttkTypeMacroCaseSimple(VTK_INT, int, 2, ttkTypeMacroRR(group0, group1, call)); \
    ttkTypeMacroCaseSimple(                                                        \
      VTK_LONG_LONG, long long, 2, ttkTypeMacroRR(group0, group1, call));    \
    ttkTypeMacroCaseSimple(                                                        \
      VTK_ID_TYPE, vtkIdType, 2, ttkTypeMacroRR(group0, group1, call));      \
    ttkTypeMacroErrorCase(2, group2);                                        \
  }
