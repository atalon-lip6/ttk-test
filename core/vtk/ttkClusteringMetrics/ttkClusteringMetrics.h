/// \ingroup vtk
/// \class ttkClusteringMetrics
/// \author Alexandre Talon <alexandre.talon@lip6.fr>
/// \date January 2022.
///
/// \brief TTK VTK-filter that wraps the ttk::ClusteringMetrics module.
///
/// This VTK filter uses the ttk::ClusteringMetrics module to compute two scores (NMI and ARI) to compare two clustering for
/// the same points.
///
/// \param Input0 vtkTable.
/// \param Input1 vtkTable.
/// \param Output vtkTable.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// See the corresponding standalone program for a usage example:
///   - standalone/ClusteringMetrics/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::ClusteringMetrics
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkClusteringMetricsModule.h>

// VTK Includes
#include <ClusteringMetrics.h>
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <ClusteringMetrics.h>

class TTKCLUSTERINGMETRICS_EXPORT ttkClusteringMetrics
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::ClusteringMetrics // and we inherit from the base class
{
private:
  /**
   * TODO 5: Add all filter parameters only as private member variables and
   *         initialize them here.
   */
  std::string OutputArrayName{"TODO"};
  // To store which column names we want to extract from the inputs.
  std::vector<std::string> ScalarFieldsHigh{}, ScalarFieldsLow{};
public:
  /**
   * TODO 6: Automatically generate getters and setters of filter
   *         parameters via vtkMacros.
   */
  vtkSetMacro(OutputArrayName, const std::string &);
  vtkGetMacro(OutputArrayName, std::string);

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkClusteringMetrics *New();
  vtkTypeMacro(ttkClusteringMetrics, ttkAlgorithm);

protected:
  /**
   * TODO 7: Implement the filter constructor and destructor
   *         (see cpp file)
   */
  ttkClusteringMetrics();
  ~ttkClusteringMetrics() override = default;


  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
