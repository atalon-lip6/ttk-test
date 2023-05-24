/// \ingroup base
/// \class ttk::TopologicalMapper
/// \author Alexandre Talon <alexandre.talon@lip6.fr>
/// \date March 2023
///
///
///TODO
/// This module defines the %TopologicalMapper class that computes a
/// score indicating how good the low dimension distance matrix represents the
/// high dimension one. The score is computed according to the SIM formula.
///
///
/// \b Related \b publication: \n
/// 'Principal Geodesic Analysis of Merge Trees (and Persistence Diagrams)'
/// Mathieu Pont, Jules Vidal and Julien Tierny.
/// IEEE Transactions on Visualization and Computer Graphics, 2023.
///
/// \sa TopologicalMapper

#pragma once

#include <vector>
// ttk common includes
#include <Debug.h>
#include <DimensionReduction.h>

namespace ttk {

  /**
   * The TopologicalMapper class TODO
   */
  class TopologicalMapper : virtual public Debug, public DimensionReduction {

  public:
    TopologicalMapper();
    ~TopologicalMapper();
    enum class LOWER_DIMENSION {
      LOWER_DIM_2D = 2,
      LOWER_DIM_3D = 3,
    };

    int execute(float* outputCoords, const std::vector<std::vector<float>> &distMatrix) const;
  protected:
    LOWER_DIMENSION LowerDimension{LOWER_DIMENSION::LOWER_DIM_2D}; //TODO enum de ttk...
  };
} // namespace ttk
