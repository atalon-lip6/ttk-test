#include <ttkClusteringMetrics.h>

#include <vtkInformation.h>

#include <vtkIntArray.h>
#include <vtkTable.h>
#include <vtkDoubleArray.h>
#include <vtkVariantArray.h> //TODO pourquoi et voir dans sim truc
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkClusteringMetrics);

ttkClusteringMetrics::ttkClusteringMetrics() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}


int ttkClusteringMetrics::FillInputPortInformation(int port, vtkInformation *info) {
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

int ttkClusteringMetrics::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}


int ttkClusteringMetrics::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  // Get input object from input vector
  vtkTable *input1 = vtkTable::GetData(inputVector[0]);
  //tester taille de l'inputVector ? => request dit taille inputvector
  vtkTable *input2 = vtkTable::GetData(inputVector[1]);
  vtkTable *output = vtkTable::GetData(outputVector);
  //tester si sortie vide
  if(!input1 || ! input2 || !output)
    return 0;
 //       During the RequestData execution one can then retrieve an actual
  //       array with the method "GetInputArrayToProcess".
  vtkAbstractArray *inputClustering1 = this->GetInputAbstractArrayToProcess(0, input1);
  vtkAbstractArray *inputClustering2 = this->GetInputAbstractArrayToProcess(1, input2);
  vtkDataArray *inputArray = this->GetInputArrayToProcess(0, input1);
  if(!inputArray) {
    this->printErr("Unable to retrieve input array.");
    return 1;
  }

  const auto intArray1 = vtkIntArray::SafeDownCast(inputClustering1);//    this->GetInputAbstractArrayToProcess(0, inputVector));
  const auto intArray2 = vtkIntArray::SafeDownCast(inputClustering2);//    this->GetInputAbstractArrayToProcess(0, inputVector));


/*
  // To make sure that the selected array can be processed by this filter,
  // one should also check that the array association and format is correct.
  if(this->GetInputArrayAssociation(0, inputVector) != 0) {
    this->printErr("Input array needs to be a point data array.");
    return 0;
  }
  if(inputArray->GetNumberOfComponents() != 1) {
    this->printErr("Input array needs to be a scalar array.");
    return 0;
  }
*/
  // If all checks pass then log which array is going to be processed.
  this->printMsg("Starting computation...");
  this->printMsg("  Scalar Array: " + std::string(inputArray->GetName()));
  // TODO ^

  size_t nbVal1 = inputClustering1->GetNumberOfTuples();
  size_t nbVal2 = inputClustering2->GetNumberOfTuples();
  this->printMsg("toto1 = " + std::to_string(nbVal1) + "\n");
  this->printMsg("toto2 = " + std::to_string(nbVal2) + "\n");
  if (nbVal1 != nbVal2)
  {
    this->printMsg("Error : the two columns must have the same size\n");
    return 0;
  }

  size_t nbVal = nbVal1;
  //TODO assert same size
  std::vector<int> values1(nbVal), values2(nbVal);
  for (size_t i = 0; i < nbVal; i++)
  {
    values1[i] = intArray1->GetValue(i);
    values2[i] = intArray2->GetValue(i);
    this->printMsg(std::to_string(values1[i]) + " " + std::to_string(values2[i]) + "\n");
  }


  double nmiValue = 0, ariValue = 0;
  //TODO calculer

  output->SetNumberOfRows(0);
  /*vtkNew<vtkDoubleArray>distorsionValArray{}, tmpCol{};
  tmpCol->SetNumberOfTuples(n);
  tmpCol->SetName("SimValue");
  for (int i = 0; i < n; i++)
    tmpCol->SetTuple1(i, vectOutput[i]);
  output->AddColumn(tmpCol);
  */

  vtkNew<vtkDoubleArray>nmiValArray{}, ariValArray{};

  nmiValArray->SetName("NMIValue");
  nmiValArray->SetNumberOfTuples(1);
  nmiValArray->SetTuple1(0, nmiValue);
  output->GetFieldData()->AddArray(nmiValArray);

  ariValArray->SetName("ARIValue");
  ariValArray->SetNumberOfTuples(1);
  ariValArray->SetTuple1(0, ariValue);
  output->GetFieldData()->AddArray(ariValArray);


  // return success
  return 1;
}
