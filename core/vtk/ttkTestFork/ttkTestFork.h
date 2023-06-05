#pragma once

// VTK Module
#include <ttkTestForkModule.h>

// VTK Includes
#include <ttkAlgorithm.h>


class TTKTESTFORK_EXPORT ttkTestFork
: public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
{
  private:
  public:
    static ttkTestFork *New();
    vtkTypeMacro(ttkTestFork, ttkAlgorithm);

  protected:
    ttkTestFork();
    ~ttkTestFork() override = default;

    int FillInputPortInformation(int port, vtkInformation *info) override;

    int FillOutputPortInformation(int port, vtkInformation *info) override;

    int RequestData(vtkInformation *request,
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector) override;
};
