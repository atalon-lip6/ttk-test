#include <ClusteringMetrics.h>
#include <cmath> // For log2 function
#include <unordered_map>
#include <vector>

ttk::ClusteringMetrics::ClusteringMetrics() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("ClusteringMetrics");
}

int ttk::ClusteringMetrics::test(void) const {
  std::vector<int> t1 = {
    2,  31, 10, 50, 19, 30, 8,  36, 25, 14, 26, 30, 28, 1,  14, 44, 29, 49, 7,
    21, 44, 21, 25, 30, 40, 20, 32, 24, 32, 25, 8,  46, 37, 6,  34, 24, 16, 41,
    13, 50, 33, 15, 42, 2,  26, 14, 10, 34, 15, 21, 49, 34, 48, 25, 10, 19, 33,
    23, 50, 20, 43, 18, 34, 36, 9,  20, 44, 37, 15, 45, 5,  4,  50, 29, 27, 23,
    2,  5,  34, 15, 12, 27, 5,  0,  42, 27, 46, 50, 4,  2,  49, 22, 20, 16, 1,
    3,  18, 1,  20, 26, 24, 30, 18, 44, 49, 43, 31, 46, 7,  21, 4,  30, 35, 8,
    12, 14, 26, 36, 27, 32, 22, 26, 21, 8,  35, 30, 1,  3,  50, 32, 49, 3,  36,
    7,  8,  49, 3,  4,  17, 18, 14, 26, 15, 0,  13, 16, 14, 18, 42, 47, 23, 9,
    9,  28, 42, 48, 49, 12, 43, 9,  24, 10, 42, 19, 4,  18, 11, 13, 29, 15, 43,
    15, 13, 42, 50, 25, 20, 8,  18, 1,  48, 13, 9,  29, 29, 41, 18, 15, 23, 9,
    10, 30, 7,  21, 26, 7,  13, 6,  32, 26, 10, 41, 18, 29, 33, 27, 0,  42, 41,
    23, 19, 42, 21, 32, 42, 40, 9,  20, 7,  40, 22, 11, 1,  3,  45, 33, 44, 45,
    5,  25, 48, 18, 33, 43, 8,  35, 14, 22, 21, 23, 39, 7,  21, 7,  42, 0,  47,
    19, 50, 28, 28, 26, 7,  18, 9,  3,  12, 2,  9,  16, 20, 8,  27, 8,  29, 31,
    43, 26, 35, 32, 26, 7,  20, 22, 34, 7,  45, 0,  43, 30, 20, 17, 23, 0,  44,
    43, 1,  14, 12, 42, 2,  43, 43, 40, 44, 11, 22, 47, 8,  23, 19, 37, 44, 32,
    0,  30, 50, 43, 9,  13, 28, 44, 6,  9,  32, 22, 11, 46, 4,  8,  38, 26, 50,
    40, 44, 20, 0,  39, 19, 9,  36, 22, 28, 48, 6,  42, 42, 44, 25, 5,  0,  31,
    6,  42, 23, 31, 47, 8,  27, 14, 46, 29, 20, 21, 33, 9,  25, 46, 43, 43, 18,
    20, 15, 38, 39, 40, 9,  30, 34, 2,  7,  7,  17, 40, 36, 11, 0,  10, 21, 44,
    48, 1,  30, 7,  11, 13, 44, 4,  1,  6,  28, 25, 36, 20, 3,  35, 37, 40, 48,
    29, 37, 17, 50, 14, 28, 5,  16, 28, 30, 41, 40, 43, 42, 22, 23, 32, 44, 3,
    34, 35, 42, 13, 14, 22, 47, 4,  21, 37, 15, 29, 12, 12, 10, 0,  50, 38, 28,
    20, 0,  35, 5,  44, 36, 27, 14, 40, 36, 6,  26, 46, 31, 45, 40, 3,  1,  3,
    35, 1,  29, 42, 12, 24, 23, 41, 35, 22, 8,  10, 11, 4,  1,  9,  5,  31, 34,
    34, 5,  27, 41, 47, 17, 39, 8,  15, 42, 16, 6,  15, 26, 25, 26, 43, 37, 42,
    31, 9,  26, 16, 26, 50, 9,  40, 5,  20, 31, 24, 9,  34, 22, 6,  18, 49, 28,
    50, 0,  40, 31, 43, 20, 35, 33, 26, 36, 14, 49, 7,  25, 17, 26, 11, 17, 29,
    48, 26, 10, 4,  33, 5,  1,  16, 44, 29, 15, 8,  18, 48, 40, 6,  16, 36, 5,
    6,  27, 38, 38, 15, 11, 12, 28, 17, 8,  50, 20, 49, 45, 18, 48, 32, 14, 49,
    47, 20, 7,  33, 3,  31, 32, 30, 30, 20, 39, 19, 17, 47, 31, 14, 35, 19, 46,
    33, 8,  6,  11, 37, 45, 27, 31, 2,  0,  41, 44, 32, 36, 46, 33, 18, 13, 8,
    41, 25, 35, 42, 23, 40, 24, 39, 42, 36, 50, 48, 8,  47, 33, 13, 22, 8,  28,
    20, 2,  24, 42, 43, 50, 48, 5,  39, 24, 18, 25, 38, 39, 44, 36, 49, 46, 19,
    19, 9,  9,  39, 30, 13, 29, 0,  45, 3,  28, 48, 3,  0,  9,  3,  44, 26, 29,
    11, 19, 30, 0,  46, 40, 44, 0,  50, 8,  0,  14, 39, 20, 2,  14, 39, 13, 41,
    10, 46, 31, 40, 2,  9,  4,  14, 43, 2,  14, 9,  16, 36, 40, 48, 14, 38, 34,
    35, 13, 42, 43, 41, 23, 14, 15, 14, 7,  8,  27, 42, 8,  8,  24, 21, 47, 46,
    39, 12, 18, 44, 6,  42, 47, 4,  32, 6,  11, 34, 8,  27, 33, 27, 20, 44, 16,
    13, 48, 42, 15, 28, 26, 21, 48, 4,  33, 30, 33, 48, 39, 48, 10, 19, 32, 17,
    49, 14, 10, 27, 32, 4,  8,  12, 5,  12, 34, 50, 27, 11, 16, 38, 6,  13, 26,
    41, 40, 31, 24, 12, 19, 26, 39, 19, 11, 21, 2,  0,  38, 50, 10, 26, 0,  31,
    29, 37, 22, 41, 13, 48, 31, 26, 3,  4,  37, 24, 16, 9,  20, 12, 37, 38, 12,
    5,  2,  45, 44, 46, 39, 3,  8,  41, 47, 5,  7,  16, 11, 7,  0,  9,  31, 2,
    14, 30, 50, 43, 31, 36, 42, 8,  39, 36, 14, 39, 43, 23, 13, 20, 22, 29, 24,
    33, 40, 44, 25, 9,  36, 13, 2,  10, 49, 17, 46, 17, 31, 11, 29, 14, 18, 36,
    31, 26, 36, 4,  12, 8,  48, 17, 34, 39, 18, 7,  46, 40, 1,  8,  26, 46, 39,
    35, 21, 21, 41, 42, 50, 9,  3,  3,  6,  5,  40, 21, 48, 17, 33, 48, 27, 23,
    16, 36, 22, 48, 41, 16, 15, 39, 42, 14, 3,  34, 0,  24, 13, 45, 2,  46, 46,
    0,  27, 32, 44, 41, 50, 18, 38, 16, 20, 50, 35, 40, 10, 10, 3,  27, 45, 41,
    44, 27, 6,  36, 31, 14, 32, 48, 3,  14, 32, 26, 12, 6,  49, 32, 29, 29, 42,
    4,  14, 17, 43, 4,  32, 40, 50, 25, 12, 23, 15, 40, 3,  37, 50, 44, 22, 20,
    49, 16, 30, 25, 2,  45, 40, 20, 33, 43, 14, 34};

  std::vector<int> t2 = {
    37, 35, 1,  42, 42, 38, 14, 34, 34, 27, 39, 12, 6,  7,  50, 10, 23, 42, 48,
    36, 50, 40, 33, 41, 45, 7,  4,  33, 36, 13, 38, 14, 29, 2,  3,  49, 33, 20,
    45, 7,  25, 22, 38, 32, 15, 18, 1,  35, 0,  50, 3,  30, 49, 24, 13, 32, 0,
    2,  38, 9,  10, 36, 14, 35, 28, 4,  37, 5,  11, 19, 34, 37, 12, 39, 28, 11,
    11, 36, 21, 11, 15, 31, 27, 27, 3,  3,  10, 47, 50, 50, 1,  42, 23, 49, 41,
    12, 22, 4,  15, 7,  50, 0,  33, 1,  11, 2,  16, 18, 10, 30, 48, 43, 5,  49,
    39, 16, 10, 19, 47, 3,  28, 0,  7,  15, 26, 45, 40, 13, 41, 36, 6,  30, 35,
    46, 36, 50, 11, 28, 18, 3,  11, 43, 37, 36, 49, 43, 28, 47, 1,  33, 1,  19,
    47, 28, 44, 1,  8,  28, 39, 15, 29, 20, 19, 19, 31, 14, 5,  47, 10, 19, 14,
    27, 4,  21, 4,  48, 47, 20, 28, 0,  36, 3,  22, 30, 20, 22, 4,  20, 28, 7,
    14, 23, 11, 14, 35, 10, 34, 19, 21, 29, 20, 28, 2,  11, 30, 22, 3,  5,  41,
    5,  35, 15, 29, 20, 7,  18, 16, 10, 35, 30, 8,  44, 49, 36, 3,  32, 36, 22,
    35, 46, 28, 49, 18, 41, 34, 15, 22, 22, 43, 45, 27, 47, 5,  22, 1,  32, 42,
    42, 50, 28, 41, 34, 38, 23, 8,  43, 22, 35, 36, 35, 42, 39, 11, 14, 45, 50,
    37, 23, 39, 26, 20, 8,  24, 43, 19, 16, 13, 14, 47, 29, 17, 42, 24, 44, 28,
    27, 31, 0,  19, 24, 7,  42, 50, 26, 6,  35, 24, 17, 23, 7,  24, 25, 47, 16,
    32, 3,  37, 49, 7,  14, 17, 29, 34, 49, 10, 42, 13, 50, 48, 36, 10, 26, 10,
    50, 4,  17, 25, 46, 38, 36, 2,  37, 26, 26, 8,  18, 18, 46, 38, 14, 42, 43,
    0,  11, 12, 14, 29, 10, 45, 42, 0,  16, 34, 17, 0,  32, 31, 11, 21, 3,  50,
    7,  37, 38, 27, 24, 1,  31, 2,  20, 50, 35, 14, 4,  9,  48, 44, 32, 17, 13,
    35, 38, 6,  23, 28, 50, 30, 21, 18, 42, 42, 14, 25, 46, 37, 9,  38, 7,  48,
    4,  45, 21, 18, 14, 16, 32, 11, 48, 43, 10, 48, 46, 24, 40, 49, 43, 8,  27,
    14, 19, 22, 14, 40, 14, 33, 25, 7,  11, 22, 17, 16, 4,  30, 41, 27, 39, 8,
    3,  38, 16, 50, 42, 14, 8,  18, 21, 40, 28, 11, 50, 21, 42, 35, 9,  35, 50,
    10, 49, 11, 15, 15, 0,  50, 17, 48, 38, 32, 4,  18, 3,  24, 36, 50, 28, 18,
    36, 14, 39, 34, 12, 33, 28, 20, 34, 12, 3,  3,  23, 10, 45, 36, 43, 29, 11,
    11, 21, 45, 24, 50, 50, 35, 42, 29, 21, 20, 2,  17, 36, 31, 13, 14, 19, 10,
    6,  41, 30, 37, 1,  49, 20, 14, 9,  14, 0,  42, 29, 6,  23, 50, 28, 0,  48,
    25, 14, 3,  36, 24, 48, 13, 41, 28, 14, 27, 8,  1,  0,  30, 10, 29, 19, 34,
    37, 26, 17, 28, 0,  2,  48, 23, 37, 41, 48, 24, 15, 45, 45, 0,  38, 12, 34,
    44, 8,  10, 22, 10, 16, 13, 29, 44, 9,  26, 30, 0,  30, 15, 7,  32, 9,  10,
    22, 13, 10, 27, 41, 44, 29, 46, 42, 40, 37, 26, 36, 46, 30, 33, 17, 42, 37,
    43, 17, 33, 49, 10, 10, 34, 48, 12, 20, 10, 6,  20, 21, 7,  24, 45, 5,  16,
    47, 33, 8,  40, 18, 32, 33, 50, 39, 26, 49, 41, 21, 29, 31, 33, 49, 9,  32,
    18, 30, 21, 33, 47, 7,  43, 28, 44, 34, 24, 7,  0,  28, 3,  29, 31, 11, 41,
    14, 48, 39, 0,  15, 16, 16, 26, 2,  45, 8,  10, 41, 37, 7,  22, 33, 18, 28,
    5,  27, 22, 8,  38, 16, 33, 4,  47, 42, 11, 34, 35, 41, 32, 15, 8,  26, 26,
    15, 3,  15, 29, 5,  33, 0,  32, 22, 24, 28, 47, 2,  8,  20, 12, 28, 7,  6,
    46, 44, 48, 19, 33, 1,  50, 25, 41, 17, 9,  23, 22, 15, 37, 38, 21, 37, 19,
    46, 11, 44, 10, 19, 13, 19, 30, 32, 28, 9,  43, 26, 46, 18, 34, 33, 38, 48,
    39, 27, 13, 49, 21, 34, 45, 35, 14, 30, 44, 40, 39, 35, 39, 41, 26, 32, 47,
    44, 29, 18, 16, 36, 10, 23, 1,  31, 42, 14, 6,  48, 28, 42, 11, 17, 18, 12,
    28, 41, 14, 21, 32, 5,  20, 30, 9,  33, 17, 23, 39, 9,  36, 11, 44, 3,  34,
    41, 23, 16, 20, 40, 10, 12, 16, 0,  5,  34, 42, 8,  26, 27, 24, 37, 13, 46,
    33, 26, 46, 12, 22, 5,  15, 10, 32, 3,  42, 9,  25, 50, 25, 36, 24, 41, 23,
    31, 34, 22, 26, 6,  25, 20, 50, 42, 48, 13, 23, 20, 11, 38, 1,  33, 46, 19,
    10, 33, 1,  24, 29, 32, 37, 0,  26, 8,  16, 39, 21, 29, 4,  4,  2,  32, 25,
    19, 31, 5,  6,  9,  6,  16, 15, 36, 36, 39, 14, 32, 37, 31, 29, 47, 7,  0,
    12, 36, 14, 34, 36, 38, 40, 46, 3,  22, 49, 0,  4,  15, 33, 0,  3,  0,  1,
    3,  46, 19, 29, 8,  12, 26, 34, 38, 49, 3,  5,  2,  45, 8,  9,  35, 21, 2,
    18, 40, 27, 23, 46, 42, 41, 32, 38, 42, 25, 22, 48, 2,  35, 50, 12, 7,  10,
    31, 41, 16, 42, 29, 32, 39, 16, 13, 7,  35, 5,  35, 4,  9,  4,  12, 2,  13,
    8,  37, 14, 9,  33, 2,  35, 0,  41, 16, 5,  9};

  this->printMsg("TOOOOOOOOOOOOOOOOOOOOOOo\n");
  double nmiValue = 0, ariValue = 0;
  execute(t1, t2, nmiValue, ariValue);
  return 1;
}

inline int nChoose2(int x) {
  return x * (x - 1) / 2;
}

int ttk::ClusteringMetrics::compute_contingency_tables(
  const std::vector<int> &clust1,
  const std::vector<int> &clust2,
  std::vector<std::vector<int>> &contingencyMatrix,
  std::vector<int> &sumLin,
  std::vector<int> &sumCol) const {
  size_t n = clust1.size();

  this->printMsg("toto0: " + std::to_string(clust1.size()) + "\n");
  this->printMsg("toto0bis: " + std::to_string(clust2.size()) + "\n");

  std::unordered_map<int, int> values1ToId, values2ToId;
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
  this->printMsg("YYYYYEEEEAAAAAHHHHH\n");
  size_t nCluster1 = nbVal1, nCluster2 = nbVal2;
  contingencyMatrix.resize(nCluster1);
  for(size_t i = 0; i < nCluster1; i++)
    contingencyMatrix[i].resize(nCluster2, 0);
  sumLin.resize(nCluster1);
  sumCol.resize(nCluster2);

  this->printMsg("toto1: " + std::to_string(contingencyMatrix.size()) + "\n");

  for(size_t i = 0; i < n; i++) {
    int x1 = values1ToId[clust1[i]], x2 = values2ToId[clust2[i]];
    contingencyMatrix[x1][x2]++;
  }
  /*
  this->printMsg("toto1.5\n");
  for (size_t i = 0; i < nCluster1; i++)
  {
    for (size_t j = 0; j < nCluster2; j++)
      this->printMsg(std::to_string(contingencyMatrix[i][j])+"\t");
    this->printMsg("\n");
  }
  */

  for(size_t i1 = 0; i1 < nCluster1; i1++) {
    int sum = 0;
    for(size_t i2 = 0; i2 < nCluster2; i2++)
      sum += contingencyMatrix[i1][i2];

    sumLin[i1] = sum;
  }
  this->printMsg("toto2\n");

  // TODO changer cache miss...
  for(size_t i2 = 0; i2 < nCluster2; i2++) {
    int sum = 0;
    for(size_t i1 = 0; i1 < nCluster1; i1++)
      sum += contingencyMatrix[i1][i2];

    sumCol[i2] = sum;
  }
  this->printMsg("toto3\n");
  return 0;
}

int ttk::ClusteringMetrics::compute_ARI(
  std::vector<std::vector<int>> &contingencyMatrix,
  const std::vector<int> &sumLin,
  const std::vector<int> &sumCol,
  int nPoint,
  double &ariValue) const {
  size_t nCluster1 = contingencyMatrix.size();
  size_t nCluster2 = contingencyMatrix[0].size(); // TODO check matrices vides
  this->printMsg(std::to_string(nCluster1) + " was nb of cluster values\n");

  double sumNChooseContingency = 0;
  for(size_t i1 = 0; i1 < nCluster1; i1++) {
    for(size_t i2 = 0; i2 < nCluster2; i2++)
      sumNChooseContingency += nChoose2(contingencyMatrix[i1][i2]);
  }
  this->printMsg(std::to_string(sumNChooseContingency)
                 + " was contingency choose sum\n");

  double sumNChoose2_1 = 0;
  double sumNChoose2_2 = 0;
  for(size_t i = 0; i < nCluster1; i++) {
    sumNChoose2_1 += nChoose2(sumLin[i]);
    this->printMsg(std::to_string(i) + "->" + std::to_string(sumLin[i]) + "\n");
  }
  for(size_t i = 0; i < nCluster2; i++) {
    sumNChoose2_2 += nChoose2(sumCol[i]);
    this->printMsg(std::to_string(i) + "->" + std::to_string(sumCol[i]) + "\n");
  }
  this->printMsg(std::to_string(sumNChoose2_1) + " was col choose sum\n");
  this->printMsg(std::to_string(sumNChoose2_2) + " was lin choose sum\n");

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
  size_t nCluster1 = contingencyMatrix.size();
  size_t nCluster2 = contingencyMatrix[0].size();

  double mutualInfo = 0;

  for(size_t i1 = 0; i1 < nCluster1; i1++) {
    for(size_t i2 = 0; i2 < nCluster2; i2++) {
      if(contingencyMatrix[i1][i2] == 0)
        continue;
      if(sumLin[i1] == 0 || sumCol[i2] == 0) // TODO ne dois plus arriver quand
                                             // généraliser étiquettes clusters
        continue;
      double logArg = (double)nPoint * contingencyMatrix[i1][i2]
                      / (sumLin[i1] * sumCol[i2]);
      double curAdd = contingencyMatrix[i1][i2] * log2(logArg) / (nPoint);
      mutualInfo += curAdd;
      this->printMsg("toto : " + std::to_string(curAdd) + "( "
                     + std::to_string(logArg));
      this->printMsg("\t : " + std::to_string(contingencyMatrix[i1][i2]) + "( "
                     + std::to_string(sumLin[i1]));
      this->printMsg("\t : " + std::to_string(contingencyMatrix[i1][i2]) + "( "
                     + std::to_string(sumCol[i2]));
    }
  }

  double entropy1 = 0, entropy2 = 0;
  for(size_t i = 0; i < nCluster1; i++) {
    double eltLin = (double)sumLin[i] / nPoint;
    // TODO plutôt assert plus tôt, je devrait pas arriver
    if(sumLin[i] != 0)
      entropy1 -= eltLin * log2(eltLin);
  }
  for(size_t i = 0; i < nCluster2; i++) {
    double eltCol = (double)sumCol[i] / nPoint;
    if(sumCol[i] != 0)
      entropy2 -= eltCol * log2(eltCol);
  }

  nmiValue = 2 * mutualInfo / (entropy1 + entropy2);
  this->printMsg("entropy1 = " + std::to_string(entropy1));
  this->printMsg("entropy2 = " + std::to_string(entropy2));
  this->printMsg("ou bien nmi = "
                 + std::to_string(2 * mutualInfo / (entropy1 + entropy2)));
  this->printMsg("ou bien mi = " + std::to_string(mutualInfo));

  // nmiValue = 1-mutualInfo/(entropy1+entropy2);
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

  /*
  this->printMsg("\n La high matrice :\n\n");
  for (size_t i = 0; i < n; i++)
  {
    this->printMsg("\t"+std::to_string(clustering1[i]));
  }
  this->printMsg("\n Et la low matrice :\n\n");
  for (size_t i = 0; i < n; i++)
  {
    this->printMsg("\t"+std::to_string(clustering2[i]));
  }
  */

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
