/// \ingroup vtk
/// \class ttkBarycentricSubdivision
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date July 2019.
///
/// \brief TTK VTK-filter that wraps the ttk::BarycentricSubdivision
/// processing package.
///
/// VTK wrapping code for the ttk::BarycentricSubdivision package.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::BarycentricSubdivision
#pragma once

// VTK Module
#include <ttkBarycentricSubdivisionModule.h>

// TTK code includes
#include <BarycentricSubdivision.h>
#include <ttkAlgorithm.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

class vtkUnstructuredGrid;

template <typename T>
class vtkSmartPointer;

class TTKBARYCENTRICSUBDIVISION_EXPORT ttkBarycentricSubdivision
  : public ttkAlgorithm,
    protected ttk::BarycentricSubdivision {

public:
  static ttkBarycentricSubdivision *New();
  vtkGetMacro(SubdivisionLevel, unsigned int);
  vtkSetMacro(SubdivisionLevel, unsigned int);

  vtkTypeMacro(ttkBarycentricSubdivision, ttkAlgorithm);

protected:
  ttkBarycentricSubdivision();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
template<size_t card>
int RequestDataTemplated(
                                           vtkInformationVector **inputVector,
                                           vtkInformationVector *outputVector,
                                           ttk::ExplicitTriangulation<card> &subdivTrig);



  /**
   * @brief Allocate an output array of same type that input array
   */
  vtkSmartPointer<vtkDataArray>
    AllocateScalarField(vtkDataArray *const inputScalarField,
                        int ntuples) const;

  int InterpolateScalarFields(vtkDataSet *const input,
                              vtkUnstructuredGrid *const output,
                              ttk::Triangulation &inputTriangulation) const;

private:
  // number of subdivisions
  unsigned int SubdivisionLevel{1};
};

template<size_t card>
int ttkBarycentricSubdivision::RequestDataTemplated(
                                           vtkInformationVector **inputVector,
                                           vtkInformationVector *outputVector,
                                           ttk::ExplicitTriangulation<card> &triangulationSubdivision) {


  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkUnstructuredGrid::GetData(outputVector);
  auto triangulation = ttkAlgorithm::GetTriangulation(input);
  this->preconditionTriangulation(triangulation);

  if(triangulation == nullptr) {
    printMsg("Error, internal triangulation is empty.");
    return 0;
  }

  // early return: copy input if no subdivision
  if(SubdivisionLevel == 0) {
    output->ShallowCopy(input);
    return 0;
  }

  // first iteration: generate the new triangulation
  int ret = this->execute(*triangulation, triangulationSubdivision);

  if(ret != 0) {
    this->printErr("Could not subdivide input mesh");
    return 0;
  }

  // first iteration: interpolate input scalar fields
  ret = InterpolateScalarFields(input, output, *triangulation);
  if(ret != 0) {
    this->printErr("Error interpolating input data array(s)");
    return 0;
  }

  for(unsigned int i = 1; i < SubdivisionLevel; ++i) {
    // move previous points to temp vector
    decltype(points_) tmpPoints{};
    std::swap(points_, tmpPoints);

    // move previous triangulation cells to temp vector
    decltype(cells_connectivity_) tmpCellsCo{};
    std::swap(cells_connectivity_, tmpCellsCo);
    decltype(cells_offsets_) tmpCellsOff{};
    std::swap(cells_offsets_, tmpCellsOff);

    // move previous triangulation to temp triangulation
    //decltype(triangulationSubdivision) tmpTr{};
    ttk::ExplicitTriangulation<card> tmpTr{};
    std::swap(triangulationSubdivision, tmpTr);

#ifdef TTK_CELL_ARRAY_NEW
    tmpTr.setInputCells(
        tmpCellsOff.size() - 1, tmpCellsCo.data(), tmpCellsOff.data());
#else
    ttk::LongSimplexId *tmpCells = nullptr;
    ttk::CellArray::TranslateToFlatLayout(tmpCellsCo, tmpCellsOff, tmpCells);
    tmpTr.setInputCells(tmpCellsCo.size() - 1, tmpCells);
#endif
    tmpTr.setInputPoints(tmpPoints.size() / 3, tmpPoints.data());
    this->preconditionTriangulation(&tmpTr);

    // generate the new triangulation
    this->execute(*triangulation, triangulationSubdivision);

    // interpolate output scalar fields
    InterpolateScalarFields(output, output, *triangulation);
  }

  // generated 3D coordinates
  auto points = vtkSmartPointer<vtkPoints>::New();
  for(size_t i = 0; i < points_.size() / 3; i++) {
    points->InsertNextPoint(&points_[3 * i]);
  }
  output->SetPoints(points);

  // generated triangles
  auto cells = vtkSmartPointer<vtkCellArray>::New();
  for(size_t i = 0; i < cells_offsets_.size() - 1; i++) {
    cells->InsertNextCell(3, &cells_connectivity_[cells_offsets_[i]]);
  }
  output->SetCells(VTK_TRIANGLE, cells);

  // cell id
  auto cellId = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  cellId->SetName("CellId");
  ttkUtils::SetVoidArray(cellId, pointId_.data(), pointId_.size(), 1);
  output->GetPointData()->AddArray(cellId);

  // cell dimension
  auto cellDim = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  cellDim->SetName("CellDimension");
  ttkUtils::SetVoidArray(cellDim, pointDim_.data(), pointDim_.size(), 1);
  output->GetPointData()->AddArray(cellDim);

  // shallow copy input field data
  output->GetFieldData()->ShallowCopy(input->GetFieldData());


  return 1;
}
