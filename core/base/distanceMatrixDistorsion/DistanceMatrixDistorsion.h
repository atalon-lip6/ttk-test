/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::DistanceMatrixDistorsion
/// \author Alexandre Talon <alexandre.talon@lip6.fr>
/// \date January 2022
///
/// This module defines the %DistanceMatrixDistorsion class that computes a score telling
/// how good the low dimension distance matrix represents the high dimension one. The score
/// is computed according to the distorsion formula.
///
///
/// TODO \b Related \b publication: \n
/// 'Principal Geodesic Analysis of Merge Trees (and Persistence Diagrams)'
/// Mathieu Pont, Jules Vidal and Julien Tierny.
/// TTK Publications.
/// 2022.
///

#pragma once

#include <vector>
// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The DistanceMatrixDistorsion class provides a method to compute the distorsion score between
   * two distance matrix representing the same points.
   */
  class DistanceMatrixDistorsion : virtual public Debug {

  public:
    DistanceMatrixDistorsion();

template <class dataType,
              class triangulationType = ttk::AbstractTriangulation>
    int computeAverages(dataType *outputData,
                        const dataType *inputData,
                        const triangulationType *triangulation) const {
      // start global timer
      ttk::Timer globalTimer;

      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator

      // print input parameters in table format
      this->printMsg({
        {"#Threads", std::to_string(this->threadNumber_)},
        {"#Vertices", std::to_string(triangulation->getNumberOfVertices())},
      });
      this->printMsg(ttk::debug::Separator::L1);

      // -----------------------------------------------------------------------
      // Compute Vertex Averages
      // -----------------------------------------------------------------------
      {
        // start a local timer for this subprocedure
        ttk::Timer localTimer;

        // print the progress of the current subprocedure (currently 0%)
        this->printMsg("Computing Averages",
                       0, // progress form 0-1
                       0, // elapsed time so far
                       this->threadNumber_, ttk::debug::LineMode::REPLACE);

        // compute the average of each vertex in parallel
        size_t nVertices = triangulation->getNumberOfVertices();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
        for(size_t i = 0; i < nVertices; i++) {
          // initialize average
          outputData[i] = inputData[i];

          // add neighbor values to average
          size_t nNeighbors = triangulation->getVertexNeighborNumber(i);
          ttk::SimplexId neighborId{-1};
          for(size_t j = 0; j < nNeighbors; j++) {
            triangulation->getVertexNeighbor(i, j, neighborId);
            outputData[i] += inputData[neighborId];
          }

          // divide by neighbor number
          outputData[i] /= (nNeighbors + 1);
        }

        // print the progress of the current subprocedure with elapsed time
        this->printMsg("Computing Averages",
                       1, // progress
                       localTimer.getElapsedTime(), this->threadNumber_);
      }

      // ---------------------------------------------------------------------
      // print global performance
      // ---------------------------------------------------------------------
      {
        this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
        this->printMsg(
          "Complete", 1, globalTimer.getElapsedTime() // global progress, time
        );
        this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
      }

      return 1; // return success
    }

    int test(int n, std::vector<double> &res) const;
    int execute(const std::vector<std::vector<double>> &highDistMatrix, const std::vector<std::vector<double>> &lowDistMatrix, double &distorsionValue, std::vector<double> &distorsionVerticesValues) const;


  }; // DistanceMatrixDistorsion class

} // namespace ttk
