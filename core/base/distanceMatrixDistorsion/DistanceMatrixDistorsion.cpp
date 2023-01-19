#include <vector>
#include <random>
#include <DistanceMatrixDistorsion.h>

ttk::DistanceMatrixDistorsion::DistanceMatrixDistorsion() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("DistanceMatrixDistorsion");
}


int ttk::DistanceMatrixDistorsion::test(int n, std::vector<double> &res) const
{
  std::vector<std::vector<double>> mat1(n), mat2(n);
  std::uniform_real_distribution<double> unif(-1000000,1000000);
  std::default_random_engine re;

  for (int i = 0; i < n; i++)
  {
    mat1[i].resize(n);
    mat2[i].resize(n);

    for (int j = 0; j < n; j++)
    {
      mat1[i][j] = unif(re);
      mat2[i][j] = unif(re);
    }
  }

  ttk::Timer timer;
  double sim;
  for (int i = 0; i < 20; i++)
    execute(mat1, mat2, sim, res);

  this->printMsg("CompleteAll", 1, timer.getElapsedTime());
  return 1;
}

int ttk::DistanceMatrixDistorsion::execute(const std::vector<std::vector<double>> &highDistMatrix, const std::vector<std::vector<double>> &lowDistMatrix, double &distorsionValue, std::vector<double> &distorsionVerticesValues) const
{
  ttk::Timer timer;
  auto n = highDistMatrix.size();

  // print horizontal separator
  this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator
                                             // print input parameters in table format
  this->printMsg({
      {"#Threads", std::to_string(this->threadNumber_)},
      {"#Vertices", std::to_string(n)},
      });
  this->printMsg(ttk::debug::Separator::L1);

#ifndef TTK_ENABLE_KAMIKAZE
  if (lowDistMatrix.size() != n)
  {
    this->printMsg(" Sizes mismatch: the high distance matrix has " + std::to_string(n) + " rows and the low distance matrix has " + std::to_string(lowDistMatrix.size()) + " rows\n.");
    return 0;
  }
#endif
  distorsionVerticesValues.resize(n);

  // Actually we first compute delta(x,y) = (d_l(x,y)-d_h(x,y))^2
  double maxi = 0;
  std::vector<std::vector<double>> deltaBis(lowDistMatrix.size());

  /*
  this->printMsg("\n La high matrice :\n\n");
  for (long unsigned int i = 0; i < n; i++)
  {
    for (long unsigned int j = 0; j < n; j++)
      this->printMsg("\t"+std::to_string(highDistMatrix[i][j]));
    this->printMsg("\n");
  }
  this->printMsg("\n Et la low matrice :\n\n");
  for (long unsigned int i = 0; i < n; i++)
  {
    for (long unsigned int j = 0; j < n; j++)
      this->printMsg("\t"+std::to_string(lowDistMatrix[i][j]));
    this->printMsg("\n");
  }
  */

#ifndef TTK_ENABLE_KAMIKAZE
  for (long unsigned int i = 0; i < n; i++)
  {
    if (highDistMatrix[i].size() != n)
    {
      this->printMsg(" Sizes mismatch: high distance matrix is not a square matrix: it has " + std::to_string(n) + " rows and  row " + std::to_string(i) + " has " + std::to_string(highDistMatrix[i].size()) + " elements.\n");
      return 0;
    }
    if (lowDistMatrix[i].size() != n)
    {
      this->printMsg(" Sizes mismatch: low distance matrix is not a square matrix: it has " + std::to_string(n) + " rows and  row " + std::to_string(i) + " has " + std::to_string(lowDistMatrix[i].size()) + "elements .\n");
      return 0;
    }
  }
#endif


  for (long unsigned int i = 0; i < n; i++)
    deltaBis[i].resize(n);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_), reduction(max:maxi)
#endif // TTK_ENABLE_OPENMP
  for (long unsigned int i = 0; i < n; i++)
  {
    for (long unsigned int j = 0; j < n; j++) // plus lent si seulement calcul triangle supérieur
    {
      double diff = lowDistMatrix[i][j] - highDistMatrix[i][j];
      deltaBis[i][j] = diff*diff;
      maxi = std::max(maxi, deltaBis[i][j]);
    }
  }
  if (maxi < 1e-8)
    maxi = 1;
  //TODO et si maxi ~= 0 ?
  // Then the actual delta'(x,y) = 1-(delta(x,y)/max{delta}).
  double totalSum = 0;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_), reduction(+:totalSum)
#endif // TTK_ENABLE_OPENMP
  for (long unsigned int i = 0; i < n; i++)
  {
    double sum = 0;
    for (long unsigned int j = 0; j < n; j++)
    {
      deltaBis[i][j] = 1-(deltaBis[i][j]/maxi);
      sum += deltaBis[i][j];
    }
    distorsionVerticesValues[i] = sum/n;
    totalSum += sum;
  }

  this->printMsg("Size of output in ttk/base = " + std::to_string(distorsionVerticesValues.size()) + "\n");

  this->printMsg("Computed distorsion value: " + std::to_string(totalSum/(n*n)));
  distorsionValue = totalSum/(n*n);
  this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
  this->printMsg("Complete", 1, timer.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator

  return 1;
}
