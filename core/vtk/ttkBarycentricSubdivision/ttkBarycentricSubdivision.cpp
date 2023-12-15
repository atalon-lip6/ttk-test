#include <ttkBarycentricSubdivision.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkBarycentricSubdivision);

ttkBarycentricSubdivision::ttkBarycentricSubdivision() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

int ttkBarycentricSubdivision::FillInputPortInformation(int port,
                                                        vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkBarycentricSubdivision::FillOutputPortInformation(int port,
                                                         vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

vtkSmartPointer<vtkDataArray> ttkBarycentricSubdivision::AllocateScalarField(
  vtkDataArray *const inputScalarField, int ntuples) const {

  vtkSmartPointer<vtkDataArray> res;

  // allocate the memory for the output scalar field
  switch(inputScalarField->GetDataType()) {
    case VTK_CHAR:
    case VTK_DOUBLE:
    case VTK_FLOAT:
    case VTK_INT:
    case VTK_ID_TYPE:
    case VTK_LONG:
      res = inputScalarField->NewInstance();
      break;
    default:
      this->printErr("Unsupported data array type");
      break;
  }
  res->SetNumberOfComponents(1);
  res->SetNumberOfTuples(ntuples);
  res->SetName(inputScalarField->GetName());
  return res;
}

int ttkBarycentricSubdivision::InterpolateScalarFields(
  vtkDataSet *const input,
  vtkUnstructuredGrid *const output,
  ttk::Triangulation &inputTriangulation) const {

  const size_t npointdata = input->GetPointData()->GetNumberOfArrays();
  const size_t ncelldata = input->GetCellData()->GetNumberOfArrays();

  const auto outPointsNumber = this->getNumberOfVertices();

  for(size_t i = 0; i < npointdata; ++i) {
    auto inputScalarField = input->GetPointData()->GetArray(i);
    if(inputScalarField == nullptr) {
      return -2;
    }

#define DISPATCH_INTERPOLATE_DIS(CASE, TYPE)                             \
  case CASE:                                                             \
    this->interpolateDiscreteScalarField<TYPE>(                          \
      static_cast<TYPE *>(ttkUtils::GetVoidPointer(inputScalarField)),   \
      static_cast<TYPE *>(ttkUtils::GetVoidPointer(outputScalarField))); \
    break
#define DISPATCH_INTERPOLATE_CONT(CASE, TYPE)                                 \
  case CASE:                                                                  \
    switch(inputTriangulation.getType()) {                                    \
      BARYSUBD_TRIANGL_CALLS(                                                 \
        TYPE, ttk::Triangulation::Type::EXPLICIT, ttk::ExplicitTriangulation<0>) \
      BARYSUBD_TRIANGL_CALLS(TYPE, ttk::Triangulation::Type::IMPLICIT,        \
                             ttk::ImplicitNoPreconditions<0>)                    \
      BARYSUBD_TRIANGL_CALLS(TYPE, ttk::Triangulation::Type::HYBRID_IMPLICIT, \
                             ttk::ImplicitWithPreconditions<0>)                  \
      BARYSUBD_TRIANGL_CALLS(                                                 \
        TYPE, ttk::Triangulation::Type::COMPACT, ttk::CompactTriangulation<0>)   \
      BARYSUBD_TRIANGL_CALLS(TYPE, ttk::Triangulation::Type::PERIODIC,        \
                             ttk::PeriodicNoPreconditions<0>)                    \
      BARYSUBD_TRIANGL_CALLS(TYPE, ttk::Triangulation::Type::HYBRID_PERIODIC, \
                             ttk::PeriodicWithPreconditions<0>)                  \
    }                                                                         \
    break;
#define BARYSUBD_TRIANGL_CALLS(DATATYPE, TRIANGL_CASE, TRIANGL_TYPE)          \
  case TRIANGL_CASE: {                                                        \
    const auto inpTri                                                         \
      = static_cast<TRIANGL_TYPE *>(inputTriangulation.getData());            \
    if(inpTri != nullptr) {                                                   \
      this->interpolateContinuousScalarField<DATATYPE, TRIANGL_TYPE>(         \
        static_cast<DATATYPE *>(ttkUtils::GetVoidPointer(inputScalarField)),  \
        static_cast<DATATYPE *>(ttkUtils::GetVoidPointer(outputScalarField)), \
        *inpTri);                                                             \
    }                                                                         \
    break;                                                                    \
  }

    auto outputScalarField
      = AllocateScalarField(inputScalarField, outPointsNumber);
    if(outputScalarField == nullptr) {
      return -3;
    }

    // only for scalar fields
    switch(inputScalarField->GetDataType()) {
      DISPATCH_INTERPOLATE_DIS(VTK_CHAR, char);
      DISPATCH_INTERPOLATE_DIS(VTK_INT, int);
      DISPATCH_INTERPOLATE_DIS(VTK_LONG, long);
      DISPATCH_INTERPOLATE_DIS(VTK_ID_TYPE, vtkIdType);
      DISPATCH_INTERPOLATE_CONT(VTK_FLOAT, float);
      DISPATCH_INTERPOLATE_CONT(VTK_DOUBLE, double);
    }
    output->GetPointData()->AddArray(outputScalarField);
  }

  const auto outCellsNumber = this->getNumberOfTriangles();

  for(size_t i = 0; i < ncelldata; ++i) {
    auto inputScalarField = input->GetCellData()->GetArray(i);
    if(inputScalarField == nullptr) {
      return -2;
    }

    auto outputScalarField
      = AllocateScalarField(inputScalarField, outCellsNumber);
    if(outputScalarField == nullptr) {
      return -3;
    }

    // only for scalar fields
    switch(inputScalarField->GetDataType()) {
      vtkTemplateMacro(this->interpolateCellDataField<VTK_TT>(
        static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalarField)),
        static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputScalarField))));
    }
    output->GetCellData()->AddArray(outputScalarField);
  }

  return 0;
}




int ttkBarycentricSubdivision::RequestData(vtkInformation *ttkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector) {

  ttk::Timer tm;

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkUnstructuredGrid::GetData(outputVector);

  auto triangulation = ttkAlgorithm::GetTriangulation(input);
  ttk::ExplicitTriangulation<0> triangulationSubdivision0{};
  ttk::ExplicitTriangulation<1> triangulationSubdivision1{};
  ttk::ExplicitTriangulation<2> triangulationSubdivision2{};
  ttk::ExplicitTriangulation<3> triangulationSubdivision3{};


  int ret;// = this->execute(*triangulation, trigOutput);
  switch(triangulation->getDimensionality()) {
    case 0:
      ret = this->RequestDataTemplated<0>(inputVector, outputVector, triangulationSubdivision0);
      break;
    case 1:
      ret = this->RequestDataTemplated<1>(inputVector, outputVector, triangulationSubdivision1);
      break;
    case 2:
      ret = this->RequestDataTemplated<2>(inputVector, outputVector, triangulationSubdivision2);
      break;
    case 3:
      ret = this->RequestDataTemplated<3>(inputVector, outputVector, triangulationSubdivision3);
      break;
    default:
      ret = 0;
      break;
  }

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), this->threadNumber_);
  return ret;
}


