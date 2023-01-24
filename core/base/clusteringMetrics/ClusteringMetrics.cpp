#include <ClusteringMetrics.h>
#include <vector>
#include <cmath> // For log2 function

ttk::ClusteringMetrics::ClusteringMetrics() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("ClusteringMetrics");
}


inline int nChoose2(int x)
{
  return x*(x-1)/2;
}

int ttk::ClusteringMetrics::compute_contingency_tables(const std::vector<double> &clust1, const std::vector<double> &clust2, std::vector<std::vector<double>> &contingencyMatrix, std::vector<double> &sumLin, std::vector<double> &sumCol) const
{
  size_t n = clust1.size();

  int max1 = *std::max_element(clust1.begin(), clust1.end());
  int max2 = *std::max_element(clust2.begin(), clust2.end());
  int nCluster = std::max(max1, max2);
  contingencyMatrix.resize(nCluster);
  for (auto i = 0; i < nCluster; i++)
    contingencyMatrix[i].resize(nCluster, 0);

  for (auto i = 0; i < n; i++)
  {
    int x1 = clust1[i], x2 = clust2[i];
    contingencyMatrix[x1][x2]++;
  }

  for (size_t i1 = 0; i1 < nCluster; i1++)
  {
    int sum = 0;
    for (size_t i2 = 0; i2 < nCluster; i2++)
      sum += contingencyMatrix[i1][i2];

    sumLin[i1] = sum;
  }

  //TODO changer cache miss...
  for (size_t i2 = 0; i2 < nCluster; i2++)
  {
    int sum = 0;
    for (size_t i1 = 0; i1 < nCluster; i1++)
      sum += contingencyMatrix[i1][i2];

    sumCol[i2] = sum;
  }
  return 0;
}


int ttk::ClusteringMetrics::compute_ARI(std::vector<std::vector<double>> &contingencyMatrix, const std::vector<double> &sumLin, const std::vector<double> &sumCol, double &ariValue) const
{
  size_t nCluster = contingencyMatrix.size();

  double sumNChooseContingency = 0;
  for (size_t i1 = 0; i1 < nCluster; i1++)
  {
    for (size_t i2 = 0; i2 < nCluster; i2++)
      sumNChooseContingency += nChoose2(contingencyMatrix[i1][i2]);
  }

  double sumNChoose2_1 = 0;
  double sumNChoose2_2 = 0;
  for (size_t i = 0; i < nCluster; i++)
  {
    sumNChoose2_1 += nChoose2(sumLin[i]);
    sumNChoose2_2 += nChoose2(sumCol[i]);
  }

  double numerator = sumNChooseContingency - (sumNChoose2_1*sumNChoose2_2)/nChoose2(nCluster);
  double denominator = 0.5*(sumNChoose2_1+sumNChoose2_2)-(sumNChoose2_1*sumNChoose2_2)/nChoose2(nCluster);

  ariValue = numerator/denominator;


  return 0;
}


int ttk::ClusteringMetrics::compute_NMI(std::vector<std::vector<double>> &contingencyMatrix, const std::vector<double> &sumLin, const std::vector<double> &sumCol, double &nmiValue) const
{
  size_t nCluster = contingencyMatrix.size();

  double mutualInfo = 0;

  for (size_t i1 = 0; i1 < nCluster; i1++)
  {
    for (size_t i2 = 0; i2 < nCluster; i2++)
    {
      double logArg = contingencyMatrix[i1][i2]/(sumLin[i1]*sumCol[i2]);
      double curAdd = contingencyMatrix[i1][i2]*log2(logArg);
      mutualInfo += curAdd;
    }
  }

  double entropy1 = 0, entropy2 = 0;
  for (size_t i = 0; i < nCluster; i++)
  {
    entropy1 -= sumLin[i]*log2(sumLin[i]);
    entropy2 -= sumCol[i]*log2(sumCol[i]);
  }

  nmiValue = 2*mutualInfo/(entropy1+entropy2);

  return 0;
}
