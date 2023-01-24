#include <ttkClusteringMetrics.h>

#include <vtkInformation.h>

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

  // Get input array that will be processed
  //
  // Note: VTK provides abstract functionality to handle array selections, but
  //       this essential functionality is unfortunately not well documented.
  //       Before you read further, please keep in mind the the TTK developer
  //       team is not responsible for the existing VTK Api ;-)
  //
  //       In a nutshell, prior to the RequestData execution one has to call
  //
  //           SetInputArrayToProcess (
  //               int idx,
  //               int port,
  //               int connection,
  //               int fieldAssociation,
  //               const char *name
  //            )
  //
  //       The parameter 'idx' is often misunderstood: lets say the filter
  //       requires n arrays, then idx enumerates them from 0 to n-1.
  //
  //       The 'port' is the input port index at which the object is connected
  //       from which we want to get the array.
  //
  //       The 'connection' is the connection index at that port (we have to
  //       specify this because VTK allows multiple connections at the same
  //       input port).
  //
  //       The 'fieldAssociation' integer specifies if the array should be taken
  //       from 0: point data, 1: cell data, or 2: field data.
  //
  //       The final parameter is the 'name' of the array.
  //
  //       Example: SetInputArrayToProcess(3,1,0,1,"EdgeLength") will store that
  //                for the 3rd array the filter needs the cell data array named
  //                "EdgeLength" that it will retrieve from the vtkDataObject
  //                at input port 1 (first connection). During the RequestData
  //                method one can then actually retrieve the 3rd array it
  //                requires for its computation by calling
  //                GetInputArrayToProcess(3, inputVector)
  //
  //       If this filter is run within ParaView, then the UI will automatically
  //       call SetInputArrayToProcess (see ClusteringMetrics.xml file).
  //
  //       During the RequestData execution one can then retrieve an actual
  //       array with the method "GetInputArrayToProcess".
  vtkDataArray *inputArray = this->GetInputArrayToProcess(0, inputVector);
  if(!inputArray) {
    this->printErr("Unable to retrieve input array.");
    return 0;
  }

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

  // If all checks pass then log which array is going to be processed.
  this->printMsg("Starting computation...");
  this->printMsg("  Scalar Array: " + std::string(inputArray->GetName()));

  // Create an output array that has the same data type as the input array
  // Note: vtkSmartPointers are well documented
  //       (https://vtk.org/Wiki/VTK/Tutorials/SmartPointers)
  vtkSmartPointer<vtkDataArray> outputArray
    = vtkSmartPointer<vtkDataArray>::Take(inputArray->NewInstance());
  outputArray->SetName(this->OutputArrayName.data()); // set array name
  outputArray->SetNumberOfComponents(1); // only one component per tuple
  outputArray->SetNumberOfTuples(inputArray->GetNumberOfTuples());

  // Get ttk::triangulation of the input vtkDataSet (will create one if one does
  // not exist already).
  /*ttk::Triangulation *triangulation
    = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation)
    return 0;
  */


  // return success
  return 1;
}
