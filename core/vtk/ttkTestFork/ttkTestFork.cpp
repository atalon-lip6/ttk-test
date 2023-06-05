#include <ttkTestFork.h>

#include <vtkInformation.h>
#include <vtkNew.h>

#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkDoubleArray.h>
#include <vtkTable.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <sys/wait.h> // for wait()
#include <unistd.h> // for fork()

vtkStandardNewMacro(ttkTestFork);

ttkTestFork::ttkTestFork() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkTestFork::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkTestFork::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
    return 1;
  }
  return 0;
}

int ttkTestFork::RequestData(vtkInformation *ttkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  if(!input)
    return 0;


  if (fork() == 0)
  {
    sleep(3);
    std::cout << "EXITING" << std::endl;
    exit(0);
  }

  wait(nullptr);
  std::cout << "YEAH" << std::endl;

  this->printMsg("Starting computation...");
  vtkNew<vtkDoubleArray> tmpCol{};
  tmpCol->SetName("TestCol");
  tmpCol->SetNumberOfTuples(7); // Why not seven? :p

  vtkTable *output = vtkTable::GetData(outputVector);
  output->AddColumn(tmpCol);
  return 1;
}
