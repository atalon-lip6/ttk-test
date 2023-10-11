/// \ingroup base
/// \class ttk::DistanceMatrixDistortion
/// \author Alexandre Talon <alexandre.talon@lip6.fr>
/// \date January 2022
///
/// This module defines the %DistanceMatrixDistortion class that computes a
/// score indicating how good the low dimension distance matrix represents the
/// high dimension one. The score is computed according to the SIM formula.
///
///
/// \b Related \b publication: \n
/// 'Principal Geodesic Analysis of Merge Trees (and Persistence Diagrams)'
/// Mathieu Pont, Jules Vidal and Julien Tierny.
/// IEEE Transactions on Visualization and Computer Graphics, 2023.
///
/// \sa DistanceMatrixDistortion

#pragma once

#include <vector>
// ttk common includes
#include <Debug.h>
#include <Geometry.h> // To check wheter a double is zero.

namespace ttk {

  /**
   * The DistanceMatrixDistortion class provides a method to compute the
   * distorsion score between two distance matrices representing the same
   * points.
   */
  template<size_t toto>
  class DistanceMatrixDistortion : virtual public Debug {

  public:
    DistanceMatrixDistortion();

    int execute(const std::vector<double *> &highDistMatrix,
                const std::vector<double *> &lowDistMatrix,
                double &distorsionValue,
                double *distorsionVerticesValues) const;
  };





} // namespace ttk
  //
  //
template<size_t toto>
ttk::DistanceMatrixDistortion<toto>::DistanceMatrixDistortion() {
  this->setDebugMsgPrefix("DistanceMatrixDistortion");
}

template<size_t toto>
int ttk::DistanceMatrixDistortion<toto>::execute(
  const std::vector<double *> &highDistMatrix,
  const std::vector<double *> &lowDistMatrix,
  double &distorsionValue,
  double *distorsionVerticesValues) const {
  ttk::Timer timer;
  auto n = highDistMatrix.size();

  if(lowDistMatrix.size() != n) {
    this->printErr(" Sizes mismatch: the high distance matrix has "
                   + std::to_string(n)
                   + " rows and the low distance matrix has "
                   + std::to_string(lowDistMatrix.size()) + " rows.");
    return 0;
  }

  /* The computation, which is optimised for performance here, can be decomposed
   * as follows: compute for each (x,y) delta(x,y) =
   * (dist_low(x,y)-dist_high(x,y))^2. Compute maxi = maximum (delta(x,y)) over
   * all (x,y). Compute delta2(x,y) = 1-delta(x,y) for each (x,y). The sim value
   * is the mean of the n^2 values of delta2(x,y). The distorsion for a vertex
   * x0 is the mean of delta2(x0,y) over all y's.
   */

  double maxi = 0;
  if(distorsionVerticesValues == nullptr) {
    this->printErr(
      " The output pointer to the distorsionValues must be non NULL. "
      "It must point to an allocated array of the right size.");
    return 1;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_) reduction(max     \
                                                                    : maxi) \
  schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < n; i++) {
    for(size_t j = i + 1; j < n; j++) {
      const double diff = lowDistMatrix[i][j] - highDistMatrix[i][j];
      maxi = std::max(maxi, diff * diff);
    }
  }
  const double EPS = ttk::Geometry::powIntTen(-DBL_DIG);
  if(maxi <= EPS) // We consider maxi is equal to zero.
  {
    this->printMsg(
      "The two distance matrices provided for SIM computation are equal.");
    maxi = 1;
  }

  double totalSum = 0;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_), reduction(+:totalSum)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < n; i++) {
    double sum = 0;
    for(size_t j = 0; j < n; j++) {
      const double diff = lowDistMatrix[i][j] - highDistMatrix[i][j];
      const double diff2 = diff * diff;
      sum += diff2;
    }
    const double sumHarmonized = sum / maxi;
    distorsionVerticesValues[i] = 1 - sumHarmonized / n;
    totalSum += 1 - sumHarmonized / n;
  }

  distorsionValue = totalSum / n;

  this->printMsg("Size of output in ttk/base = " + std::to_string(n));

  this->printMsg("Computed distorsion value: "
                 + std::to_string(distorsionValue));
  this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
  this->printMsg("Complete", 1, timer.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator

  return 0;
}
