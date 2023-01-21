#include <ClusteringMetrics.h>
#include <cmath> // For the log2 function
#include <map>
#include <vector>

ttk::ClusteringMetrics::ClusteringMetrics() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("ClusteringMetrics");
}

inline int nChoose2(int x) {
  return x * (x - 1) / 2;
}

inline int checkContingencyMatSize(const ttk::ClusteringMetrics *object,
                                   const std::vector<std::vector<int>> &matrix,
                                   int nPoint) {
#ifndef TTK_ENABLE_KAMIKAZE
  if(nPoint == 0) {
    object->printErr("Error: clustering on zero points.");
    return 0;
  }

  size_t nLin = matrix.size();
  if(nLin == 0) {
    object->printErr("The provided contingency matrix is empty.\n");
    return 0;
  }
  size_t nCol = matrix[0].size();

  for(size_t i = 0; i < nLin; i++) {
    size_t curNCol = matrix[i].size();
    if(curNCol == 0) {
      object->printErr("Line " + std::to_string(i)
                       + " of the contingency matrix is empty.\n");
      return 0;
    } else if(curNCol != nCol) {
      object->printErr(
        "Line " + std::to_string(i)
        + " of the contingency matrix has wrong number of columns : "
        + std::to_string(curNCol) + " instead of " + std::to_string(nCol)
        + ".\n");
      return 0;
    }
  }
#endif // TTK_ENABLE_KAMIKAZE
  return 1;
}

int ttk::ClusteringMetrics::compute_contingency_tables(
  const std::vector<int> &clust1,
  const std::vector<int> &clust2,
  std::vector<std::vector<int>> &contingencyMatrix,
  std::vector<int> &sumLin,
  std::vector<int> &sumCol) const {
  size_t nPoint = clust1.size();
#ifndef TTK_ENABLE_KAMIKAZE
  if(nPoint == 0) {
    this->printErr("Error: clustering on zero points.");
    return 0;
  }
#endif

  std::map<int, int> values1ToId, values2ToId;
  size_t nbVal1 = 0, nbVal2 = 0;
  for(int x : clust1) {
    auto found = values1ToId.find(x);
    if(found == values1ToId.end()) {
      values1ToId[x] = nbVal1;
      nbVal1++;
    }
  }
  for(int x : clust2) {
    auto found = values2ToId.find(x);
    if(found == values2ToId.end()) {
      values2ToId[x] = nbVal2;
      nbVal2++;
    }
  }

  size_t nCluster1 = nbVal1, nCluster2 = nbVal2;
  contingencyMatrix.resize(nCluster1);
  for(size_t i = 0; i < nCluster1; i++)
    contingencyMatrix[i].resize(nCluster2, 0);
  sumLin.resize(nCluster1);
  sumCol.resize(nCluster2, 0);

  for(size_t i = 0; i < nPoint; i++) {
    int x1 = values1ToId[clust1[i]], x2 = values2ToId[clust2[i]];
    contingencyMatrix[x1][x2]++;
  }

  for(size_t i1 = 0; i1 < nCluster1; i1++) {
    int sum = 0;
    for(size_t i2 = 0; i2 < nCluster2; i2++) {
      sumCol[i2] += contingencyMatrix[i1][i2];
      sum += contingencyMatrix[i1][i2];
    }

    sumLin[i1] = sum;
  }

  return 0;
}

int ttk::ClusteringMetrics::compute_ARI(
  std::vector<std::vector<int>> &contingencyMatrix,
  const std::vector<int> &sumLin,
  const std::vector<int> &sumCol,
  int nPoint,
  double &ariValue) const {
  if(!checkContingencyMatSize(this, contingencyMatrix, nPoint))
    return 0;
  size_t nCluster1 = contingencyMatrix.size();
  size_t nCluster2 = contingencyMatrix[0].size();

  double sumNChooseContingency = 0;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_) reduction(+:sumNChooseContingency)
#endif // TTK_ENABLE_OPENMP
  for(size_t i1 = 0; i1 < nCluster1; i1++) {
    for(size_t i2 = 0; i2 < nCluster2; i2++)
      sumNChooseContingency += nChoose2(contingencyMatrix[i1][i2]);
  }

  double sumNChoose2_1 = 0, sumNChoose2_2 = 0;
  for(size_t i = 0; i < nCluster1; i++) {
#ifndef TTK_ENABLE_KAMIKAZE
    if(sumLin[i] == 0) {
      this->printErr("Error: the sum of a line in the contingency matrix is "
                     "zero. This should not happen.");
    }
#endif
    sumNChoose2_1 += nChoose2(sumLin[i]);
  }
  for(size_t i = 0; i < nCluster2; i++) {
#ifndef TTK_ENABLE_KAMIKAZE
    if(sumCol[i] == 0) {
      this->printErr("Error: the sum of a column in the contingency matrix is "
                     "zero. This should not happen.");
    }
#endif
    sumNChoose2_2 += nChoose2(sumCol[i]);
  }

  double numerator = sumNChooseContingency
                     - (sumNChoose2_1 * sumNChoose2_2) / nChoose2(nPoint);
  double denominator = 0.5 * (sumNChoose2_1 + sumNChoose2_2)
                       - (sumNChoose2_1 * sumNChoose2_2) / nChoose2(nPoint);
  ariValue = numerator / denominator;

  return 0;
}

int ttk::ClusteringMetrics::compute_NMI(
  std::vector<std::vector<int>> &contingencyMatrix,
  const std::vector<int> &sumLin,
  const std::vector<int> &sumCol,
  int nPoint,
  double &nmiValue) const {
  if(!checkContingencyMatSize(this, contingencyMatrix, nPoint))
    return 0;
  size_t nCluster1 = contingencyMatrix.size();
  size_t nCluster2 = contingencyMatrix[0].size();

  double mutualInfo = 0;
  bool invalidCell = false;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_) reduction(+:mutualInfo)
#endif // TTK_ENABLE_OPENMP
  for(size_t i1 = 0; i1 < nCluster1; i1++) {
    for(size_t i2 = 0; i2 < nCluster2; i2++) {
      if(contingencyMatrix[i1][i2] == 0)
        continue;
#ifndef TTK_ENABLE_KAMIKAZE
      if(sumLin[i1] == 0 || sumCol[i2] == 0) {
        this->printErr("Error: a sum of a line or a column of the contingency "
                       "matrix is zero. This should not happen.");
        invalidCell = true;
        continue;
      }
#endif
      double logArg = (double)nPoint * contingencyMatrix[i1][i2]
                      / (sumLin[i1] * sumCol[i2]);
      double curAdd = contingencyMatrix[i1][i2] * log2(logArg) / (nPoint);
      mutualInfo += curAdd;
    }
  }
  if(invalidCell)
    return 0;

  double entropy1 = 0, entropy2 = 0;
  for(size_t i = 0; i < nCluster1; i++) {
    double eltLin = (double)sumLin[i] / nPoint;
    entropy1 -= eltLin * log2(eltLin);
  }
  for(size_t i = 0; i < nCluster2; i++) {
    double eltCol = (double)sumCol[i] / nPoint;
    entropy2 -= eltCol * log2(eltCol);
  }

  nmiValue = 2 * mutualInfo / (entropy1 + entropy2);

  return 0;
}

int ttk::ClusteringMetrics::execute(const std::vector<int> &clustering1,
                                    const std::vector<int> &clustering2,
                                    double &nmiValue,
                                    double &ariValue) const {
  ttk::Timer timer;
  size_t n = clustering1.size();

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

  if(clustering2.size() != n) {
    this->printMsg(" Sizes mismatch: clustering one represents "
                   + std::to_string(n)
                   + " points and clustering two represents "
                   + std::to_string(clustering2.size()) + " points\n");
    return 1;
  }
#endif

  std::vector<std::vector<int>> contingencyMatrix;
  std::vector<int> sumLines, sumColumns;
  compute_contingency_tables(
    clustering1, clustering2, contingencyMatrix, sumLines, sumColumns);

  compute_ARI(contingencyMatrix, sumLines, sumColumns, n, ariValue);
  compute_NMI(contingencyMatrix, sumLines, sumColumns, n, nmiValue);

#ifndef TTK_ENABLE_KAMIKAZE
  this->printMsg("Size of output in ttk/base = 0\n");

  this->printMsg("Computed NMI value: " + std::to_string(nmiValue) + "\n");
  this->printMsg("Computed ARI value: " + std::to_string(ariValue) + "\n");
  this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
  this->printMsg("Complete", 1, timer.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
#endif
  return 0;
}
