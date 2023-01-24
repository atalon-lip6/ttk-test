#include <ttkDistanceMatrixDistorsion.h>

#include <vtkInformation.h>

#include <vtkTable.h>
#include <vtkDoubleArray.h>
#include <vtkVariantArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <regex>

#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkDistanceMatrixDistorsion);

ttkDistanceMatrixDistorsion::ttkDistanceMatrixDistorsion() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

int ttkDistanceMatrixDistorsion::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  }
  else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkDistanceMatrixDistorsion::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkDistanceMatrixDistorsion::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  // Get input objects from input vector
  vtkTable *inputHigh = vtkTable::GetData(inputVector[0]);
  //tester taille de l'inputVector ? => request dit taille inputvector
  vtkTable *inputLow = vtkTable::GetData(inputVector[1]);

  vtkTable *output = vtkTable::GetData(outputVector);
  //tester si sortie vide

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputLow || !inputHigh || !output)
    return 0;
#endif

  if(SelectFieldsWithRegexpHigh) {
    // select all input columns whose name is matching the regexp
    ScalarFieldsHigh.clear();
    const auto n = inputHigh->GetNumberOfColumns();
    for(int i = 0; i < n; ++i) {
      const auto &name = inputHigh->GetColumnName(i);
      if(std::regex_match(name, std::regex(RegexpStringHigh))) {
        ScalarFieldsHigh.emplace_back(name);
      }
    }
  }

  if(SelectFieldsWithRegexpLow) {
    // select all input columns whose name is matching the regexp
    ScalarFieldsLow.clear();
    const size_t n = inputLow->GetNumberOfColumns();
    for(size_t i = 0; i < n; ++i) {
      const auto &name = inputLow->GetColumnName(i);
      if(std::regex_match(name, std::regex(RegexpStringLow))) {
        ScalarFieldsLow.emplace_back(name);
      }
    }
  }


  const size_t nRowsHigh = inputHigh->GetNumberOfRows();
  const size_t nRowsLow = inputLow->GetNumberOfRows();

  const size_t nColsHigh = ScalarFieldsHigh.size();
  const size_t nColsLow = ScalarFieldsLow.size();

#ifndef TTK_ENABLE_KAMIKAZE
  if (nRowsHigh <= 0 || nRowsHigh != nColsHigh) {
    this->printErr("High input matrix cannot be a valid square matrix (rows: "
                   + std::to_string(nRowsHigh)
                   + ", columns: " + std::to_string(nColsHigh) + ")");
    return 0;
  }

  if (nRowsLow <= 0 || nRowsLow != nColsLow) {
    this->printErr("Low input matrix is not a valid square matrix (rows: "
                   + std::to_string(nRowsLow)
                   + ", columns: " + std::to_string(nColsLow) + ")");
    return 0;
  }
  if (nRowsHigh != nRowsLow) {
    this->printErr("High and low input matrices must have same size (rows(high): "
                   + std::to_string(nRowsHigh)
                   + ", rows(low): " + std::to_string(nRowsLow) + ")");
    return 0;
  }
#endif

  int n = nRowsHigh;
  std::vector<std::vector<double>> vectMatHigh(n), vectMatLow(n);
  //TODO séparer le resize du remplissage pour utiliser pragma openmp
/*
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
*/

  std::vector<vtkAbstractArray *> arraysHigh{}, arraysLow{};
  for (const auto &s : ScalarFieldsHigh)
  {
    arraysHigh.push_back(inputHigh->GetColumnByName(s.data()));
    //this->printMsg(" coucou high: " + std::string(s.data()) + "\n");
  }
  for (const auto &s : ScalarFieldsLow)
  {
    arraysLow.push_back(inputLow->GetColumnByName(s.data()));
    //this->printMsg(" coucou low: " + std::string(s.data()) + "\n");
  }

  for (int i = 0; i < n; i++)
  {
    vectMatHigh[i].resize(n);
    vectMatLow[i].resize(n);
    for (int j = 0; j < n; j++)
    {
      vectMatHigh[i][j] = arraysHigh[j]->GetVariantValue(i).ToDouble();
      vectMatLow[i][j] = arraysLow[j]->GetVariantValue(i).ToDouble();
    }
  }

  std::vector<double> vectOutput;
  double distorsionValue = 0;
  this->printMsg("Starting computation of sim distorsion value...");
  this->test(10000, vectOutput);
  //this->execute(vectMatHigh, vectMatLow, distorsionValue, vectOutput);
  //this->printMsg(std::to_string(n) + " VS " + std::to_string(vectOutput.size()) + "\n");


  output->SetNumberOfRows(n);
  vtkNew<vtkDoubleArray>distorsionValArray{}, tmpCol{};
  tmpCol->SetNumberOfTuples(n);
  tmpCol->SetName("SimValue");
  for (int i = 0; i < n; i++)
    tmpCol->SetTuple1(i, vectOutput[i]);
  output->AddColumn(tmpCol);

  distorsionValArray->SetName("DistorsionValue");
  distorsionValArray->SetNumberOfTuples(1);
  distorsionValArray->SetTuple1(0, distorsionValue);
  output->GetFieldData()->AddArray(distorsionValArray);

  // return success
  return 1;
}
