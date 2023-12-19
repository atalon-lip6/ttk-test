#include <ttkTemplatingTest.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

using namespace ttk;

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkTemplatingTest);

ttkTemplatingTest::ttkTemplatingTest() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

/**
 * TODO 8: Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkTemplatingTest::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}


int ttkTemplatingTest::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkTemplatingTest::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  auto input = vtkDataSet::GetData(inputVector[0]);
  //auto output = vtkDataSet::GetData(outputVector);

  Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);

/*
  vtkDataArray *inputArray = this->GetInputArrayToProcess(0, inputVector);
  if(!inputArray) {
    this->printErr("Unable to retrieve input array.");
    return 0;
  }

  this->preconditionTriangulation(triangulation);

  vtkDataArray *inputScalarField = this->GetInputArrayToProcess(0, inputVector);
  if(!inputScalarField)
    return 0;

  if(inputScalarField->GetNumberOfComponents() != 1) {
    printErr("Invalid scalar field ("
             + std::to_string(inputScalarField->GetNumberOfComponents())
             + " components)");
    return 0;
  }
*/
  // If all checks pass then log which array is going to be processed.
  this->printMsg("Starting computation...");
  //this->printMsg("  Scalar Array: " + std::string(inputArray->GetName()));

  // Create an output array that has the same data type as the input array
  // Note: vtkSmartPointers are well documented
  //       (https://vtk.org/Wiki/VTK/Tutorials/SmartPointers)
  vtkSmartPointer<vtkDataArray> const outputArray;
    /*= vtkSmartPointer<vtkDataArray>::Take(inputArray->NewInstance());
  outputArray->SetName(this->OutputArrayName.data()); // set array name
  outputArray->SetNumberOfComponents(1); // only one component per tuple
  outputArray->SetNumberOfTuples(inputArray->GetNumberOfTuples());
*/
  // Get ttk::triangulation of the input vtkDataSet (will create one if one does
  // not exist already).
/*
  ttk::Triangulation *triangulation
    = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation)
    return 0;
*/
  // Precondition the triangulation (e.g., enable fetching of vertex neighbors)
  this->preconditionTriangulation(triangulation); // implemented in base class


  // On error cancel filter execution

  // Get output vtkDataSet (which was already instantiated based on the
  // information provided by FillOutputPortInformation)
  vtkDataSet *outputDataSet = vtkDataSet::GetData(outputVector, 0);

  // make a SHALLOW copy of the input
  //outputDataSet->ShallowCopy(inputDataSet);

  // add to the output point data the computed output array
  outputDataSet->GetPointData()->AddArray(outputArray);

  // calling the smoothing package
  ttkTypeMacroAT(
    VTK_INT, triangulation->getType(), triangulation->getDimensionality(),
    (execute<T0, T1>(
      static_cast<const T1 *>(triangulation->getData()), NumberOfIterations)));

  return 1;
}
