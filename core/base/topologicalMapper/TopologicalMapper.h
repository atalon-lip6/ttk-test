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

    int execute(const std::vector<std::vector<double>> distMatrix, unsigned int nDim, std::vector<std::vector<double>> &outputCoords);
  };
} // namespace ttk
