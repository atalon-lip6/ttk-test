/// \ingroup base
/// \class ttk::ClusteringMetrics
/// \author Alexandre Talon <alexandre.talon@lip6.fr>
/// \date January 2022
///
/// This module defines the %ClusteringMetrics class that computes both the NMI
/// (normalized mutual information) and the ARI (adjusted rand index) scores
/// between two clustering of the same points.
///
/// \b Related \b publication: \n
/// 'Principal Geodesic Analysis of Merge Trees (and Persistence Diagrams)'
/// Mathieu Pont, Jules Vidal and Julien Tierny.
/// TTK Publications.
/// 2022.
///
/// \sa ClusteringMetrics

#pragma once

#include <vector>
// ttk common includes
#include <Debug.h>

namespace ttk {

  /**
   * The ClusteringMetrics class provides methods to compute two scores
   * comparing two clustering of the same points. It can compute the NMI score
   * and the ARI score.
   */
  class ClusteringMetrics : virtual public Debug {

  public:
    ClusteringMetrics();

    int execute(const std::vector<int> &clustering1,
                const std::vector<int> &clustering2,
                double &nmiValue,
                double &ariValue) const;

    int computeContingencyTables(
      const std::vector<int> &clust1,
      const std::vector<int> &clust2,
      std::vector<std::vector<int>> &contingencyMatrix,
      std::vector<int> &sumLin,
      std::vector<int> &sumCol) const;

    // For single ARI or single NMI computing, the contingency matrix, as well
    // as the sums of each lines and each column must be provided as arguments.
    int computeARI(std::vector<std::vector<int>> &contingencyMatrix,
                    const std::vector<int> &sumLin,
                    const std::vector<int> &sumCol,
                    int nPoint,
                    double &ariValue) const;
    int computeNMI(std::vector<std::vector<int>> &contingencyMatrix,
                    const std::vector<int> &sumLin,
                    const std::vector<int> &sumCol,
                    int nPoint,
                    double &nmiValue) const;
  }; // ClusteringMetrics class

} // namespace ttk
