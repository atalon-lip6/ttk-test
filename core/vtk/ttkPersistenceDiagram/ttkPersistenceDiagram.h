/// \ingroup vtk
/// \class ttkPersistenceDiagram
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date September 2016.
///
/// \brief TTK VTK-filter for the computation of persistence diagrams.
///
/// This filter computes the persistence diagram of the extremum-saddle pairs
/// of an input scalar field. The X-coordinate of each pair corresponds to its
/// birth, while its smallest and highest Y-coordinates correspond to its birth
/// and death respectively.
///
/// In practice, the diagram is represented by a vtkUnstructuredGrid. Each
/// vertex of this mesh represent a critical point of the input data. It is
/// associated with point data (vertexId, critical type). Each vertical edge
/// of this mesh represent a persistence pair. It is associated with cell data
/// (persistence of the pair, critical index of the extremum of the pair).
///
/// Persistence diagrams are useful and stable concise representations of the
/// topological features of a data-set. It is useful to fine-tune persistence
/// thresholds for topological simplification or for fast similarity
/// estimations for instance.
///
/// \param Input Input scalar field, either 2D or 3D, regular grid or
/// triangulation (vtkDataSet)
/// \param Output Output persistence diagram (vtkUnstructuredGrid)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \b Related \b publication \n
/// "Computational Topology: An Introduction" \n
/// Herbert Edelsbrunner and John Harer \n
/// American Mathematical Society, 2010
///
/// \sa ttkFTMTreePP
/// \sa ttkPersistenceCurve
/// \sa ttkScalarFieldCriticalPoints
/// \sa ttkTopologicalSimplification
/// \sa ttk::PersistenceDiagram
#ifndef _TTK_PERSISTENCEDIAGRAM_H
#define _TTK_PERSISTENCEDIAGRAM_H

// VTK includes
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkFloatArray.h>
#include <vtkTable.h>

// VTK Module
#include <ttkPersistenceDiagramModule.h>

// ttk code includes
#include <PersistenceDiagram.h>
#include <ttkAlgorithm.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

class TTKPERSISTENCEDIAGRAM_EXPORT ttkPersistenceDiagram
  : public ttkAlgorithm,
    protected ttk::PersistenceDiagram  {

public:
  static ttkPersistenceDiagram *New();

  vtkTypeMacro(ttkPersistenceDiagram, ttkAlgorithm);

  void SetForceInputOffsetScalarField(int data) {
    ForceInputOffsetScalarField = data;
    Modified();
    computeDiagram_ = true;
  }
  vtkGetMacro(ForceInputOffsetScalarField, int);

  void SetComputeSaddleConnectors(int data) {
    ComputeSaddleConnectors = data;
    Modified();
    computeDiagram_ = true;
  }
  vtkGetMacro(ComputeSaddleConnectors, int);

  void SetShowInsideDomain(int onOff) {
    ShowInsideDomain = onOff;
    Modified();
  }
  vtkGetMacro(ShowInsideDomain, int);

  template <typename scalarType, class triangulationType>
  int setPersistenceDiagramInfo(
    ttk::SimplexId id,
    vtkSmartPointer<ttkSimplexIdTypeArray> vertexIdentifierScalars,
    vtkSmartPointer<vtkIntArray> nodeTypeScalars,
    vtkSmartPointer<vtkFloatArray> coordsScalars,
    const std::vector<std::tuple<ttk::SimplexId,
                                 ttk::CriticalType,
                                 ttk::SimplexId,
                                 ttk::CriticalType,
                                 scalarType,
                                 ttk::SimplexId>> &diagram,
    vtkSmartPointer<vtkPoints> points,
    vtkIdType ids[3], vtkDataArray *inputScalars,
  const triangulationType *triangulation);

  template <typename scalarType, class triangulationType>
  int getPersistenceDiagram(
    vtkUnstructuredGrid *outputCTPersistenceDiagram,
    ttk::ftm::TreeType treeType,
    const std::vector<std::tuple<ttk::SimplexId,
                                 ttk::CriticalType,
                                 ttk::SimplexId,
                                 ttk::CriticalType,
                                 scalarType,
                                 ttk::SimplexId>> &diagram, vtkDataArray *inputScalars,
  const triangulationType *triangulation);

  template <typename scalarType, class triangulationType>
  int setPersistenceDiagramInfoInsideDomain(
    ttk::SimplexId id,
    vtkSmartPointer<ttkSimplexIdTypeArray> vertexIdentifierScalars,
    vtkSmartPointer<vtkIntArray> nodeTypeScalars,
    vtkDataArray *birthScalars,
    vtkDataArray *deathScalars,
    const std::vector<std::tuple<ttk::SimplexId,
                                 ttk::CriticalType,
                                 ttk::SimplexId,
                                 ttk::CriticalType,
                                 scalarType,
                                 ttk::SimplexId>> &diagram,
    vtkSmartPointer<vtkPoints> points,
    vtkIdType ids[3], vtkDataArray *inputScalars,
  const triangulationType *triangulation);

  template <typename scalarType, class triangulationType>
  int getPersistenceDiagramInsideDomain(
    vtkUnstructuredGrid *outputCTPersistenceDiagram,
    ttk::ftm::TreeType treeType,
    const std::vector<std::tuple<ttk::SimplexId,
                                 ttk::CriticalType,
                                 ttk::SimplexId,
                                 ttk::CriticalType,
                                 scalarType,
                                 ttk::SimplexId>> &diagram, vtkDataArray *inputScalars,
  const triangulationType *triangulation);

  template <typename VTK_TT, typename TTK_TT>
  int dispatch(vtkUnstructuredGrid *outputCTPersistenceDiagram,
                                  vtkDataArray *inputScalarDataArray,
                                  const VTK_TT *inputScalars,
                                  int inputOffsetsDataType,
                                  const void *inputOffsets,
                                  const TTK_TT *triangulation);

  template <typename VTK_TT>
  int deleteDiagram();
  
protected:
  ttkPersistenceDiagram();
  ~ttkPersistenceDiagram() override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

private:
  bool ForceInputOffsetScalarField{false};
  int ShowInsideDomain{false};
  
  bool computeDiagram_{true};
  void *CTDiagram_{nullptr};
  int scalarDataType{0};
};

template <typename scalarType, class triangulationType>
int ttkPersistenceDiagram::setPersistenceDiagramInfo(
  ttk::SimplexId id,
  vtkSmartPointer<ttkSimplexIdTypeArray> vertexIdentifierScalars,
  vtkSmartPointer<vtkIntArray> nodeTypeScalars,
  vtkSmartPointer<vtkFloatArray> coordsScalars,
  const std::vector<std::tuple<ttk::SimplexId,
                               ttk::CriticalType,
                               ttk::SimplexId,
                               ttk::CriticalType,
                               scalarType,
                               ttk::SimplexId>> &diagram,
  vtkSmartPointer<vtkPoints> points,
  vtkIdType ids[3], const scalarType *inputScalars,
  const triangulationType *triangulation) {
  double p[3] = {0, 0, 0};
  const ttk::SimplexId a = std::get<0>(diagram[id]);
  const ttk::SimplexId na
    = static_cast<ttk::SimplexId>(std::get<1>(diagram[id]));
  const ttk::SimplexId b = std::get<2>(diagram[id]);
  const ttk::SimplexId nb
    = static_cast<ttk::SimplexId>(std::get<3>(diagram[id]));

  nodeTypeScalars->InsertTuple1(2 * id, na);
  nodeTypeScalars->InsertTuple1(2 * id + 1, nb);

  vertexIdentifierScalars->InsertTuple1(2 * id, a);
  vertexIdentifierScalars->InsertTuple1(2 * id + 1, b);

  float coords[3];
  triangulation->getVertexPoint(a, coords[0], coords[1], coords[2]);
  coordsScalars->InsertTuple3(2 * id, coords[0], coords[1], coords[2]);

  triangulation->getVertexPoint(b, coords[0], coords[1], coords[2]);
  coordsScalars->InsertTuple3(2 * id + 1, coords[0], coords[1], coords[2]);

  p[0] = inputScalars->GetTuple1(a);
  p[1] = inputScalars->GetTuple1(a);
  ids[0] = points->InsertNextPoint(p);

  p[0] = inputScalars->GetTuple1(a);
  p[1] = inputScalars->GetTuple1(b);
  ids[1] = points->InsertNextPoint(p);

  return 0;
}

template <typename scalarType, class triangulationType>
int ttkPersistenceDiagram::getPersistenceDiagram(
  vtkUnstructuredGrid *outputCTPersistenceDiagram,
  ttk::ftm::TreeType treeType,
  const std::vector<std::tuple<ttk::SimplexId,
                               ttk::CriticalType,
                               ttk::SimplexId,
                               ttk::CriticalType,
                               scalarType,
                               ttk::SimplexId>> &diagram,
                               vtkDataArray *inputScalars,
  const triangulationType *triangulation) {
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram
    = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<ttkSimplexIdTypeArray> vertexIdentifierScalars
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  vertexIdentifierScalars->SetNumberOfComponents(1);
  vertexIdentifierScalars->SetName(ttk::VertexScalarFieldName);

  vtkSmartPointer<vtkIntArray> nodeTypeScalars
    = vtkSmartPointer<vtkIntArray>::New();
  nodeTypeScalars->SetNumberOfComponents(1);
  nodeTypeScalars->SetName("CriticalType");

  vtkSmartPointer<ttkSimplexIdTypeArray> pairIdentifierScalars
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  pairIdentifierScalars->SetNumberOfComponents(1);
  pairIdentifierScalars->SetName("PairIdentifier");

  vtkSmartPointer<vtkDoubleArray> persistenceScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  persistenceScalars->SetNumberOfComponents(1);
  persistenceScalars->SetName("Persistence");

  vtkSmartPointer<vtkIntArray> extremumIndexScalars
    = vtkSmartPointer<vtkIntArray>::New();
  extremumIndexScalars->SetNumberOfComponents(1);
  extremumIndexScalars->SetName("PairType");

  vtkSmartPointer<vtkFloatArray> coordsScalars
    = vtkSmartPointer<vtkFloatArray>::New();
  coordsScalars->SetNumberOfComponents(3);
  coordsScalars->SetName("Coordinates");

  const ttk::SimplexId minIndex = 0;
  const ttk::SimplexId saddleSaddleIndex = 1;
  const ttk::SimplexId maxIndex = triangulation->getCellVertexNumber(0) - 2;

  const ttk::SimplexId diagramSize = diagram.size();
  if(diagramSize) {
    vtkIdType ids[2];
    vtkIdType oldIds[2];

    scalarType maxPersistenceValue = std::numeric_limits<scalarType>::min();
    oldIds[0] = 0;
    for(ttk::SimplexId i = 0; i < diagramSize; ++i) {
      const scalarType persistenceValue = std::get<4>(diagram[i]);
      const ttk::SimplexId type = std::get<5>(diagram[i]);
      maxPersistenceValue = std::max(persistenceValue, maxPersistenceValue);

      setPersistenceDiagramInfo(i, vertexIdentifierScalars, nodeTypeScalars,
                                coordsScalars, diagram, points, ids, inputScalars, triangulation);

      // add cell data
      persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);
      pairIdentifierScalars->InsertTuple1(i, i);
      if(!i)
        extremumIndexScalars->InsertTuple1(i, -1);
      else {
        switch(type) {
          case 0:
            extremumIndexScalars->InsertTuple1(i, minIndex);
            break;

          case 1:
            extremumIndexScalars->InsertTuple1(i, saddleSaddleIndex);
            break;

          case 2:
            extremumIndexScalars->InsertTuple1(i, maxIndex);
            break;
        }
      }
      persistenceScalars->InsertTuple1(i, persistenceValue);
    }
    oldIds[1] = ids[0];

    // add diag
    persistenceDiagram->InsertNextCell(VTK_LINE, 2, oldIds);
    pairIdentifierScalars->InsertTuple1(diagramSize, -1);
    extremumIndexScalars->InsertTuple1(diagramSize, -1);
    persistenceScalars->InsertTuple1(diagramSize, 2 * maxPersistenceValue);
  }

  persistenceDiagram->SetPoints(points);
  persistenceDiagram->GetPointData()->AddArray(vertexIdentifierScalars);
  persistenceDiagram->GetPointData()->AddArray(nodeTypeScalars);
  persistenceDiagram->GetPointData()->AddArray(coordsScalars);
  persistenceDiagram->GetCellData()->AddArray(pairIdentifierScalars);
  persistenceDiagram->GetCellData()->AddArray(extremumIndexScalars);
  persistenceDiagram->GetCellData()->AddArray(persistenceScalars);

  outputCTPersistenceDiagram->ShallowCopy(persistenceDiagram);

  return 0;
}

template <typename scalarType, class triangulationType>
int ttkPersistenceDiagram::setPersistenceDiagramInfoInsideDomain(
  ttk::SimplexId id,
  vtkSmartPointer<ttkSimplexIdTypeArray> vertexIdentifierScalars,
  vtkSmartPointer<vtkIntArray> nodeTypeScalars,
  vtkDataArray *birthScalars,
  vtkDataArray *deathScalars,
  const std::vector<std::tuple<ttk::SimplexId,
                               ttk::CriticalType,
                               ttk::SimplexId,
                               ttk::CriticalType,
                               scalarType,
                               ttk::SimplexId>> &diagram,
  vtkSmartPointer<vtkPoints> points,
  vtkIdType ids[3], vtkDataArray *inputScalars,
  const triangulationType *triangulation) {
  float p[3];
  const ttk::SimplexId a = std::get<0>(diagram[id]);
  const ttk::SimplexId na
    = static_cast<ttk::SimplexId>(std::get<1>(diagram[id]));
  const ttk::SimplexId b = std::get<2>(diagram[id]);
  const ttk::SimplexId nb
    = static_cast<ttk::SimplexId>(std::get<3>(diagram[id]));
  const double sa = inputScalars->GetTuple1(a);
  const double sb = inputScalars->GetTuple1(b);

  nodeTypeScalars->InsertTuple1(2 * id, na);
  nodeTypeScalars->InsertTuple1(2 * id + 1, nb);
  vertexIdentifierScalars->InsertTuple1(2 * id, a);
  vertexIdentifierScalars->InsertTuple1(2 * id + 1, b);
  birthScalars->InsertTuple1(2 * id, sa);
  birthScalars->InsertTuple1(2 * id + 1, sa);
  deathScalars->InsertTuple1(2 * id, sa);
  deathScalars->InsertTuple1(2 * id + 1, sb);

  triangulation->getVertexPoint(a, p[0], p[1], p[2]);
  ids[0] = points->InsertNextPoint(p);

  triangulation->getVertexPoint(b, p[0], p[1], p[2]);
  ids[1] = points->InsertNextPoint(p);

  return 0;
}

template <typename scalarType, class triangulationType>
int ttkPersistenceDiagram::getPersistenceDiagramInsideDomain(
  vtkUnstructuredGrid *outputCTPersistenceDiagram,
  ttk::ftm::TreeType treeType,
  const std::vector<std::tuple<ttk::SimplexId,
                               ttk::CriticalType,
                               ttk::SimplexId,
                               ttk::CriticalType,
                               scalarType,
                               ttk::SimplexId>> &diagram, vtkDataArray *inputScalars,
  const triangulationType *triangulation) {
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram
    = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<ttkSimplexIdTypeArray> vertexIdentifierScalars
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  vertexIdentifierScalars->SetNumberOfComponents(1);
  vertexIdentifierScalars->SetName(ttk::VertexScalarFieldName);

  vtkSmartPointer<vtkIntArray> nodeTypeScalars
    = vtkSmartPointer<vtkIntArray>::New();
  nodeTypeScalars->SetNumberOfComponents(1);
  nodeTypeScalars->SetName("CriticalType");

  vtkSmartPointer<ttkSimplexIdTypeArray> pairIdentifierScalars
    = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
  pairIdentifierScalars->SetNumberOfComponents(1);
  pairIdentifierScalars->SetName("PairIdentifier");

  vtkSmartPointer<vtkDoubleArray> persistenceScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  persistenceScalars->SetNumberOfComponents(1);
  persistenceScalars->SetName("Persistence");

  vtkSmartPointer<vtkIntArray> extremumIndexScalars
    = vtkSmartPointer<vtkIntArray>::New();
  extremumIndexScalars->SetNumberOfComponents(1);
  extremumIndexScalars->SetName("PairType");

  vtkDataArray *birthScalars = inputScalars->NewInstance();
  birthScalars->SetNumberOfComponents(1);
  birthScalars->SetName("Birth");

  vtkDataArray *deathScalars = inputScalars->NewInstance();
  deathScalars->SetNumberOfComponents(1);
  deathScalars->SetName("Death");

  const ttk::SimplexId minIndex = 0;
  const ttk::SimplexId saddleSaddleIndex = 1;
  const ttk::SimplexId maxIndex = triangulation->getCellVertexNumber(0) - 2;

  const ttk::SimplexId diagramSize = diagram.size();
  if(diagramSize) {
    vtkIdType ids[2];

    scalarType maxPersistenceValue = std::numeric_limits<scalarType>::min();
    for(ttk::SimplexId i = 0; i < diagramSize; ++i) {
      const scalarType persistenceValue = std::get<4>(diagram[i]);
      const ttk::SimplexId type = std::get<5>(diagram[i]);
      maxPersistenceValue = std::max(persistenceValue, maxPersistenceValue);

      setPersistenceDiagramInfoInsideDomain(i, vertexIdentifierScalars,
                                            nodeTypeScalars, birthScalars,
                                            deathScalars, diagram, points, ids, inputScalars, triangulation);

      // add cell data
      persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);
      pairIdentifierScalars->InsertTuple1(i, i);
      if(!i)
        extremumIndexScalars->InsertTuple1(i, -1);
      else {
        switch(type) {
          case 0:
            extremumIndexScalars->InsertTuple1(i, minIndex);
            break;

          case 1:
            extremumIndexScalars->InsertTuple1(i, saddleSaddleIndex);
            break;

          case 2:
            extremumIndexScalars->InsertTuple1(i, maxIndex);
            break;
        }
      }
      persistenceScalars->InsertTuple1(i, persistenceValue);
    }
  }

  persistenceDiagram->SetPoints(points);
  persistenceDiagram->GetPointData()->AddArray(vertexIdentifierScalars);
  persistenceDiagram->GetPointData()->AddArray(nodeTypeScalars);
  persistenceDiagram->GetPointData()->AddArray(birthScalars);
  persistenceDiagram->GetPointData()->AddArray(deathScalars);
  persistenceDiagram->GetCellData()->AddArray(pairIdentifierScalars);
  persistenceDiagram->GetCellData()->AddArray(extremumIndexScalars);
  persistenceDiagram->GetCellData()->AddArray(persistenceScalars);

  outputCTPersistenceDiagram->ShallowCopy(persistenceDiagram);
  return 0;
}

#endif // _TTK_PERSISTENCEDIAGRAM_H
