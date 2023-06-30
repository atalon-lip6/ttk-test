#include <ttkTopologicalMapper.h>
#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAttributes.h>
#include <vtkTable.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkNew.h>

#include <ttkMacros.h>
#include <ttkUtils.h>
#include <ttkMapper.h>

#include <regex>
#include <string>

vtkStandardNewMacro(ttkTopologicalMapper);


static void extractInputMatrix(ttk::Mapper::Matrix &inputMatrix,
                               std::vector<std::vector<float>> &distMatrix,
                               std::vector<std::string> &arrayNames,
                               vtkDataSetAttributes *const pd,
                               const std::string &regexp,
                               const int nThreads,
                               const bool useRegexp) {
  if(pd == nullptr) {
    return;
  }

  if(useRegexp) {
    // select all input columns whose name is matching the regexp
    arrayNames.clear();
    const auto n = pd->GetNumberOfArrays();
    arrayNames.reserve(n);
    for(int i = 0; i < n; ++i) {
      const auto name = pd->GetArrayName(i);

      if(name == nullptr) {
        continue;
      }

      // check regexp
      if(!std::regex_match(name, std::regex(regexp))) {
        continue;
      }

      // sanity check
      const auto arr{pd->GetArray(name)};
      if(arr == nullptr || arr->GetNumberOfComponents() != 1) {
        continue;
      }

      // store matching array names
      arrayNames.emplace_back(name);
    }
  }

  std::vector<vtkDataArray *> arrays(arrayNames.size());
  for(size_t i = 0; i < arrayNames.size(); ++i) {
    arrays[i] = pd->GetArray(arrayNames[i].data());
  }

  inputMatrix.alloc(pd->GetNumberOfTuples(), arrayNames.size());
//#ifdef TTK_ENABLE_OPENMP
//#pragma omp parallel for num_threads(nThreads)
//#endif // TTK_ENABLE_OPENMP
  //for(size_t i = 0; i < inputMatrix.nRows(); ++i) {
    //for(size_t j = 0; j < inputMatrix.nCols(); ++j) {
  for(size_t i = 0; i < distMatrix.size(); ++i) {
    for(size_t j = 0; j < distMatrix.size(); ++j) {
      const auto array{arrays[j]};
      distMatrix[i][j] = array->GetVariantValue(i).ToDouble();
      //inputMatrix.get(i, j) = array->GetVariantValue(i).ToDouble();
    }
  }
}

ttkTopologicalMapper::ttkTopologicalMapper() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}


int ttkTopologicalMapper::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable"); // ou table
    return 1;
  }
  return 0;
}


int ttkTopologicalMapper::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set( vtkDataObject::DATA_TYPE_NAME(), "vtkTable" );//table ?
    return 1;
  }
  return 0;
}
/*
static void setNodes(vtkUnstructuredGrid *const nodes,
                     const std::vector<std::array<float, 3>> &compBaryCoords,
                     const std::vector<int> &compBucketId) {

  vtkNew<vtkIntArray> compId{};
  compId->SetName("ComponentId");
  compId->SetNumberOfComponents(1);
  compId->SetNumberOfTuples(compBaryCoords.size());

  vtkNew<vtkIntArray> bucket{};
  bucket->SetName("BucketId");
  bucket->SetNumberOfComponents(1);
  bucket->SetNumberOfTuples(compBucketId.size());

  vtkNew<vtkPoints> points{};
  points->SetNumberOfPoints(compBaryCoords.size());
  for(size_t i = 0; i < compBaryCoords.size(); ++i) {
    points->SetPoint(i, compBaryCoords[i].data());
    bucket->SetTuple1(i, compBucketId[i]);
    compId->SetTuple1(i, i);
  }

  ttkUtils::CellVertexFromPoints(nodes, points);
  nodes->GetPointData()->AddArray(bucket);
  nodes->GetPointData()->SetScalars(compId);
}
*/





int ttkTopologicalMapper::RequestData(vtkInformation *ttkNotUsed(request),
                           vtkInformationVector **inputVector,
                           vtkInformationVector *outputVector) {


  auto input = vtkTable::GetData(inputVector[0]);
  auto outputNodes = vtkTable::GetData(outputVector, 0);


  /*
  Matrix inputDistMat{};
  if(this->ReEmbedMapper) {
    ttk::Timer tm{};
    extractInputMatrix(inputDistMat, this->DistanceMatrix,
                       input->GetPointData(), this->DistanceMatrixRegexp,
                       this->threadNumber_, this->SelectMatrixWithRegexp);
    if(inputDistMat.nRows() != inputDistMat.nCols() || inputDistMat.nCols() == 0
       || inputDistMat.nRows() == 0) {
      this->printErr("Invalid input distance matrix");
      return 0;
    }
    */
  size_t nbPoint = input->GetNumberOfRows();
  vtkNew<vtkPoints> outputPoints{};
  outputPoints->SetNumberOfPoints(nbPoint);
  std::cout << "IL Y A " << nbPoint << " POINTS :D\n";
  outputPoints->GetData()->Fill(0.0);
  //outputPoints->SetData(input->GetPointData);

  std::vector<std::array<float, 3>> compBaryCoords{};

  std::vector<std::vector<float>> distMatrix(nbPoint);
  for (int i = 0; i < nbPoint; i++)
    distMatrix[i].resize(nbPoint);
  int status;


  /*auto *pts = input->GetRowData();
  vtkDataArray* ScalarFields = pts->GetScalars();
  size_t dim = ScalarFields.size();
  vector<double> inputCoords(dim*nbPoint);
  std::vector<vtkAbstractArray *> arrays{};
  for(const auto &s : ScalarFields) {
    arrays.push_back(input->GetColumnByName(s.data()));
  }
  */
  char colNames[2][2] = {"x", "y"};
  vtkDoubleArray *colX = vtkDoubleArray::SafeDownCast(input->GetColumnByName("x")),
                   *colY = vtkDoubleArray::SafeDownCast(input->GetColumnByName("y"));
      //inputMatrix[i*dim] = (arrays[j]->GetVariantValue(i).ToDouble());
  std::vector<double> inputCoords(nbPoint*2);
  for (size_t i = 0; i < nbPoint; i++)
  {
    inputCoords[2*i] = colX->GetValue(i);
    inputCoords[2*i+1] = colY->GetValue(i);
  }
  ttk::Mapper::Matrix inputDistMat{};
  ttk::Timer tm{};
  extractInputMatrix(inputDistMat,
                     distMatrix, this->DistanceMatrix,
      input->GetRowData(), this->DistanceMatrixRegexp,
      this->threadNumber_, this->SelectMatrixWithRegexp);
  //if(inputDistMat.nRows() != inputDistMat.nCols() || inputDistMat.nCols() == 0
  //    || inputDistMat.nRows() == 0) {
  if (false) { //TODO
    this->printErr("Invalid input distance matrix");
    return 0;
  }
  //if (false)
  this->execute(inputCoords, ttkUtils::GetPointer<float>(outputPoints->GetData()), distMatrix);
  std::string nameCoords[3] = {"x", "y", "z"};
  for (int i = 0; i < 2; i++)
  {
    vtkNew<vtkDoubleArray> col{};
    col->SetNumberOfTuples(nbPoint);
    col->SetName(nameCoords[i].c_str());
    for (int j = 0; j < nbPoint; j++)
    {
      //std::cout << j  << " => " << outputPoints->GetData()[3*j+i];
      col->SetTuple1(j, ttkUtils::GetPointer<float>(outputPoints->GetData())[2*j+i]);
    }

    outputNodes->AddColumn(col);
  }
  //outputNodes->GetPointData()->SetNumberOfPoints(nbPoint);
  //outputNodes->SetPoints(outputPoints);
  //setNodes(outputNodes, compBaryCoords, compBucketId);
  //setArcs(outputArcs, compArcs, outputNodes->GetPoints());

  return 1;
}
