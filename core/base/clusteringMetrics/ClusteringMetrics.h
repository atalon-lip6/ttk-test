/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::ClusteringMetrics
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// This module defines the %ClusteringMetrics class that computes for each vertex of a
/// triangulation the average scalar value of itself and its direct neighbors.
///
/// \b Related \b publication: \n
/// 'ClusteringMetrics'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

#include <vector>

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The ClusteringMetrics class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class ClusteringMetrics : virtual public Debug {

  public:
    ClusteringMetrics();

  int execute(const std::vector<int> &clustering1, const std::vector<int> &clustering2, double &nmiValue, double &ariValue) const;

  int compute_contingency_tables(const std::vector<int> &clust1, const std::vector<int> &clust2, std::vector<std::vector<int>> &contingencyMatrix, std::vector<int> &sumLin, std::vector<int> &sumCol) const;


  int compute_ARI(std::vector<std::vector<int>> &contingencyMatrix, const std::vector<int> &sumLin, const std::vector<int> &sumCol, int nPoint, double &ariValue) const;
  int compute_NMI(std::vector<std::vector<int>> &contingencyMatrix, const std::vector<int> &sumLin, const std::vector<int> &sumCol, double &nmiValue) const;
  }; // ClusteringMetrics class

} // namespace ttk
