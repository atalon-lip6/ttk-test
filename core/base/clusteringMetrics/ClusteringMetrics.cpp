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

int ttk::ClusteringMetrics::compute_contingency_tables(const std::vector<int> &clust1, const std::vector<int> &clust2, std::vector<std::vector<int>> &contingencyMatrix, std::vector<int> &sumLin, std::vector<int> &sumCol) const
{
  size_t n = clust1.size();

  this->printMsg("toto0: "+std::to_string(clust1.size())+"\n");
  this->printMsg("toto0bis: "+std::to_string(clust2.size())+"\n");

  size_t max1 = *std::max_element(clust1.begin(), clust1.end());
  size_t max2 = *std::max_element(clust2.begin(), clust2.end());
  size_t nCluster = std::max(max1, max2)+1;
  contingencyMatrix.resize(nCluster);
  for (size_t i = 0; i < nCluster; i++)
    contingencyMatrix[i].resize(nCluster, 0);
  sumLin.resize(nCluster);
  sumCol.resize(nCluster);

  this->printMsg("toto1: "+std::to_string(contingencyMatrix.size())+"\n");

  for (size_t i = 0; i < n; i++)
  {
    int x1 = clust1[i], x2 = clust2[i];
    contingencyMatrix[x1][x2]++;
  }
  this->printMsg("toto1.5\n");
  for (int i = 0; i < nCluster; i++)
  {
    for (int j = 0; j < nCluster; j++)
      this->printMsg(std::to_string(contingencyMatrix[i][j])+"\t");
    this->printMsg("\n");
  }

  for (size_t i1 = 0; i1 < nCluster; i1++)
  {
    int sum = 0;
    for (size_t i2 = 0; i2 < nCluster; i2++)
      sum += contingencyMatrix[i1][i2];

    sumLin[i1] = sum;
  }
  this->printMsg("toto2\n");

  //TODO changer cache miss...
  for (size_t i2 = 0; i2 < nCluster; i2++)
  {
    int sum = 0;
    for (size_t i1 = 0; i1 < nCluster; i1++)
      sum += contingencyMatrix[i1][i2];

    sumCol[i2] = sum;
  }
  this->printMsg("toto2\n");
  return 0;
}


int ttk::ClusteringMetrics::compute_ARI(std::vector<std::vector<int>> &contingencyMatrix, const std::vector<int> &sumLin, const std::vector<int> &sumCol, int nPoint, double &ariValue) const
{
  size_t nCluster = contingencyMatrix.size();
  this->printMsg(std::to_string(nCluster) + " was nb of cluster values\n");

  double sumNChooseContingency = 0;
  for (size_t i1 = 0; i1 < nCluster; i1++)
  {
    for (size_t i2 = 0; i2 < nCluster; i2++)
      sumNChooseContingency += nChoose2(contingencyMatrix[i1][i2]);
  }
  this->printMsg(std::to_string(sumNChooseContingency) + " was contingency choose sum\n");

  double sumNChoose2_1 = 0;
  double sumNChoose2_2 = 0;
  for (size_t i = 0; i < nCluster; i++)
  {
    this->printMsg(std::to_string(i) + "->"+std::to_string(sumCol[i])+"\n");
    sumNChoose2_1 += nChoose2(sumCol[i]);
    sumNChoose2_2 += nChoose2(sumLin[i]);
  }
  this->printMsg(std::to_string(sumNChoose2_1) + " was col choose sum\n");
  this->printMsg(std::to_string(sumNChoose2_2) + " was lin choose sum\n");

  double numerator = sumNChooseContingency - (sumNChoose2_1*sumNChoose2_2)/nChoose2(nPoint);
  double denominator = 0.5*(sumNChoose2_1+sumNChoose2_2)-(sumNChoose2_1*sumNChoose2_2)/nChoose2(nPoint);

  ariValue = numerator/denominator;


  return 0;
}


int ttk::ClusteringMetrics::compute_NMI(std::vector<std::vector<int>> &contingencyMatrix, const std::vector<int> &sumLin, const std::vector<int> &sumCol, double &nmiValue) const
{
  size_t nCluster = contingencyMatrix.size();

  double mutualInfo = 0;

  for (size_t i1 = 0; i1 < nCluster; i1++)
  {
    for (size_t i2 = 0; i2 < nCluster; i2++)
    {
      if (contingencyMatrix[i1][i2] == 0)
        continue;
      double logArg = contingencyMatrix[i1][i2]/(sumLin[i1]*sumCol[i2]);
      double curAdd = contingencyMatrix[i1][i2]*log2(logArg);
      mutualInfo += curAdd;
    }
  }

  double entropy1 = 0, entropy2 = 0;
  for (size_t i = 0; i < nCluster; i++)
  {
    //TODO plutôt assert plus tôt, je devrait pas arriver
    if (sumLin[i] != 0)
      entropy1 -= sumLin[i]*log2(sumLin[i]);
    if (sumCol[i] != 0)
      entropy2 -= sumCol[i]*log2(sumCol[i]);
  }

  nmiValue = 2*mutualInfo/(entropy1+entropy2);

  return 0;
}


  int ttk::ClusteringMetrics::execute(const std::vector<int> &clustering1, const std::vector<int> &clustering2, double &nmiValue, double &ariValue) const
{
  ttk::Timer timer;
  size_t n = clustering1.size();

#ifndef TTK_ENABLE_KAMIKAZE
  // print horizontal separator
  this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator
                                             // print input parameters in table format
  this->printMsg({
      {"#Threads", std::to_string(this->threadNumber_)},
      {"#Vertices", std::to_string(n)},
      });
  this->printMsg(ttk::debug::Separator::L1);

  if (clustering2.size() != n)
  {
    this->printMsg(" Sizes mismatch: clustering one represents " + std::to_string(n) + " points and clustering two represents " + std::to_string(clustering2.size()) + " points\n");
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
  compute_contingency_tables(clustering1, clustering2, contingencyMatrix, sumLines, sumColumns);

  compute_ARI(contingencyMatrix, sumLines, sumColumns, n, ariValue);
  compute_NMI(contingencyMatrix, sumLines, sumColumns, nmiValue);

#ifndef TTK_ENABLE_KAMIKAZE
  this->printMsg("Size of output in ttk/base = 2\n");

  this->printMsg("Computed NMI value: " + std::to_string(nmiValue)+"\n");
  this->printMsg("Computed ARI value: " + std::to_string(ariValue)+"\n");
  this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
  this->printMsg("Complete", 1, timer.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
#endif
  return 0;
}
