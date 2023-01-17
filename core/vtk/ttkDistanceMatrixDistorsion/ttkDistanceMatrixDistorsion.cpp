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

/**
 * TODO 7: Implement the filter constructor and destructor in the cpp file.
 *
 * The constructor has to specify the number of input and output ports
 * with the functions SetNumberOfInputPorts and SetNumberOfOutputPorts,
 * respectively. It should also set default values for all filter
 * parameters.
 *
 * The destructor is usually empty unless you want to manage memory
 * explicitly, by for example allocating memory on the heap that needs
 * to be freed when the filter is destroyed.
 */
ttkDistanceMatrixDistorsion::ttkDistanceMatrixDistorsion() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

/**
 * TODO 8: Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
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

/**
 * TODO 9: Specify the data object type of each output port
 *
 * This method specifies in the port information object the data type of the
 * corresponding output objects. It is possible to either explicitly
 * specify a type by adding a vtkDataObject::DATA_TYPE_NAME() key:
 *
 *      info->Set( vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
 *
 * or to pass a type of an input port to an output port by adding the
 * ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key (see below).
 *
 * Note: prior to the execution of the RequestData method the pipeline will
 * initialize empty output data objects based on this information.
 */
int ttkDistanceMatrixDistorsion::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

/**
 * TODO 10: Pass VTK data to the base code and convert base code output to VTK
 *
 * This method is called during the pipeline execution to update the
 * already initialized output data objects based on the given input
 * data objects and filter parameters.
 *
 * Note:
 *     1) The passed input data objects are validated based on the information
 *        provided by the FillInputPortInformation method.
 *     2) The output objects are already initialized based on the information
 *        provided by the FillOutputPortInformation method.
 */
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
    const auto n = inputLow->GetNumberOfColumns();
    for(int i = 0; i < n; ++i) {
      const auto &name = inputLow->GetColumnName(i);
      if(std::regex_match(name, std::regex(RegexpStringLow))) {
        ScalarFieldsLow.emplace_back(name);
      }
    }
  }


  const auto nRowsHigh = inputHigh->GetNumberOfRows();
  const auto nRowsLow = inputLow->GetNumberOfRows();

  const auto nColsHigh = ScalarFieldsHigh.size()
  const auto nColsLow = ScalarFieldsLow.size()

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
    this->printMsg(" coucou high: " + std::string(s.data()) + "\n");
  }
  for (const auto &s : ScalarFieldsLow)
  {
    arraysLow.push_back(inputLow->GetColumnByName(s.data()));
    this->printMsg(" coucou low: " + std::string(s.data()) + "\n");
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
  double distorsionValue;
  this->printMsg("Starting computation of sim distorsion value...");
  this->execute(vectMatHigh, vectMatLow, distorsionValue, vectOutput);
  this->printMsg(std::to_string(n) + " VS " + std::to_string(vectOutput.size()) + "\n");


  output->SetNumberOfRows(n);
  vtkNew<vtkDoubleArray>distorsionValArray{}, tmpCol{};
  this->printMsg("totot1\n");
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
