#include <DistanceMatrixDistorsion.h>
#include <random>
#include <vector>
// To check wheter a double is zero.
#include <Geometry.h>
#include <Os.h>

ttk::DistanceMatrixDistorsion::DistanceMatrixDistorsion() {
  this->setDebugMsgPrefix("DistanceMatrixDistorsion");
}

int ttk::DistanceMatrixDistorsion::execute(
  const std::vector<std::vector<double>> &highDistMatrix,
  const std::vector<std::vector<double>> &lowDistMatrix,
  double &distorsionValue,
  std::vector<double> &distorsionVerticesValues) const {
  ttk::Timer timer;
  auto n = highDistMatrix.size();

#ifndef TTK_ENABLE_KAMIKAZE
  // print horizontal separator
  this->printMsg(
    ttk::debug::Separator::L1); // L1 is the '=' separator
                                // print input parameters in table format
  this->printMsg({
    {"#Threads", std::to_string(this->threadNumber_)},
    {"#Vertices", std::to_string(n)},
  });
  this->printMsg(ttk::debug::Separator::L1);

  if(lowDistMatrix.size() != n) {
    this->printErr(" Sizes mismatch: the high distance matrix has "
                   + std::to_string(n)
                   + " rows and the low distance matrix has "
                   + std::to_string(lowDistMatrix.size()) + " rows\n.");
    return 0;
  }

  for(size_t i = 0; i < n; i++) {
    if(highDistMatrix[i].size() != n) {
      this->printErr(
        " Sizes mismatch: high distance matrix is not a square matrix: it has "
        + std::to_string(n) + " rows and  row " + std::to_string(i) + " has "
        + std::to_string(highDistMatrix[i].size()) + " elements.\n");
      return 1;
    }
    if(lowDistMatrix[i].size() != n) {
      this->printErr(
        " Sizes mismatch: low distance matrix is not a square matrix: it has "
        + std::to_string(n) + " rows and  row " + std::to_string(i) + " has "
        + std::to_string(lowDistMatrix[i].size()) + "elements .\n");
      return 1;
    }
  }
#endif

  /* The computation, which is optimised for performance here, can be decomposed
   * as follows: compute for each (x,y) delta(x,y) =
   * (dist_low(x,y)-dist_high(x,y))^2. Compute maxi = maximum (delta(x,y)) over
   * all (x,y). Compute delta2(x,y) = 1-delta(x,y) for each (x,y). The sim value
   * is the mean of the n^2 values of delta2(x,y). The distorsion for a vertex
   * x0 is the mean of delta2(x0,y) over all y's.
   */

  double maxi = 0;
  distorsionVerticesValues.resize(n);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_) reduction(max     \
                                                                    : maxi) \
  schedule(dynamic)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < n; i++) {
    for(size_t j = i + 1; j < n; j++) {
      double diff = lowDistMatrix[i][j] - highDistMatrix[i][j];
      maxi = std::max(maxi, diff * diff);
    }
  }
  const double EPS = ttk::Geometry::powIntTen(-DBL_DIG);
  if(maxi <= EPS) // We consider maxi is equal to zero.
  {
#ifndef TTK_ENABLE_KAMIKAZE
    this->printMsg(
      "The two distance matrices provided for SIM computation are equal.\n");
#endif
    maxi = 1;
  }

  double totalSum = 0;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_), reduction(+:totalSum)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < n; i++) {
    double sum = 0;
    for(size_t j = 0; j < n; j++) {
      double diff = lowDistMatrix[i][j] - highDistMatrix[i][j];
      double diff2 = diff * diff;
      sum += diff2;
    }
    double sumHarmonized = sum / maxi;
    distorsionVerticesValues[i] = 1 - sumHarmonized / n;
    totalSum += 1 - sumHarmonized / n;
  }

  distorsionValue = totalSum / n;

#ifndef TTK_ENABLE_KAMIKAZE
  this->printMsg("Size of output in ttk/base = "
                 + std::to_string(distorsionVerticesValues.size()) + "\n");

  this->printMsg("Computed distorsion value: "
                 + std::to_string(distorsionValue));
  this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
  this->printMsg("Complete", 1, timer.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
#endif

  return 0;
}
