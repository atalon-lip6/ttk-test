/// \ingroup base
/// \class ttk::DistanceMatrixDistorsion
/// \author Alexandre Talon <alexandre.talon@lip6.fr>
/// \date January 2022
///
/// This module defines the %DistanceMatrixDistorsion class that computes a
/// score telling how good the low dimension distance matrix represents the high
/// dimension one. The score is computed according to the SIM formula.
///
///
/// TODO \b Related \b publication: \n
/// 'Principal Geodesic Analysis of Merge Trees (and Persistence Diagrams)'
/// Mathieu Pont, Jules Vidal and Julien Tierny.
/// TTK Publications.
/// 2022.
///
/// \sa DistanceMatrixDistorsion
///
/// TODO \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/clusteringKelvinHelmholtzInstabilities/">
///   Clustering Kelvin Helmholtz Instabilities example</a> \n

#pragma once

#include <vector>
// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The DistanceMatrixDistorsion class provides a method to compute the
   * distorsion score between two distance matrices representing the same
   * points.
   */
  class DistanceMatrixDistorsion : virtual public Debug {

  public:
    DistanceMatrixDistorsion();

    int test(int n, std::vector<double> &res) const;
    int execute(const std::vector<std::vector<double>> &highDistMatrix,
                const std::vector<std::vector<double>> &lowDistMatrix,
                double &distorsionValue,
                std::vector<double> &distorsionVerticesValues) const;
  };

} // namespace ttk
