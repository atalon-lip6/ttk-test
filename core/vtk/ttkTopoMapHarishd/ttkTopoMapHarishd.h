/// TODO 4: Provide your information and **update** the documentation (in
/// particular regarding the order convention if input arrays need to be
/// specified with the standard VTK call SetInputArrayToProcess()).
///
/// \ingroup vtk
/// \class ttkTopoMapHarishd
/// \author Alexandre Talon <alexandre.talon@lip6.fr>
/// \date June 2023
///
/// \brief TTK VTK-filter that wraps the ttk::TopoMapHarishd module.
///
/// This class generates a mapper from a data-set. The mapper is computed using some minimum spanning tree
/// and convex hull arguments.
///
/// \param Input vtkDataSet.
/// \param Output vtkDataSet.
///

///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).

#pragma once

// VTK Module
#include <ttkTopoMapHarishdModule.h>

// ttk code includes
#include <TopoMapHarishd.h>
#include <ttkAlgorithm.h>
#include <ttkMacros.h>
#include <string>

  //vtkTypeMacro(ttkTopoMapHarishd, ttkAlgorithm);

class TTKTOPOMAPHARISHD_EXPORT ttkTopoMapHarishd
  : public ttkAlgorithm , protected ttk::TopoMapHarishd {
private:
  /**
   * TODO 5: Add all filter parameters only as private member variables and
   *         initialize them here.
   */
  std::string OutputArrayName{"AveragedScalarField"};

public:

  vtkSetMacro(SelectMatrixWithRegexp, bool);
  vtkGetMacro(SelectMatrixWithRegexp, bool);

  vtkSetMacro(DistanceMatrixRegexp, const std::string &);
  vtkGetMacro(DistanceMatrixRegexp, std::string);

  void SetDistMat(const std::string &s) {
    this->DistanceMatrixNames.emplace_back(s);
    this->Modified();
  }
  void ClearDistMat() {
    this->DistanceMatrixNames.clear();
    this->Modified();
  }

  vtkSetMacro(AngleSamplingFreq, int);
  vtkGetMacro(AngleSamplingFreq, int);
  
  vtkSetMacro(CheckMST, bool);
  vtkGetMacro(CheckMST, bool);

  static ttkTopoMapHarishd *New();
  vtkTypeMacro(ttkTopoMapHarishd, ttkAlgorithm);

protected:

  ttkTopoMapHarishd();
  ~ttkTopoMapHarishd() override = default;


  int FillInputPortInformation(int port, vtkInformation *info) override;


  int FillOutputPortInformation(int port, vtkInformation *info) override;


  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  // domain vertices coordinates in high dimension
  bool SelectMatrixWithRegexp{false};
  std::string DistanceMatrixRegexp{".*"};
  std::vector<std::string> DistanceMatrixNames{};
};
