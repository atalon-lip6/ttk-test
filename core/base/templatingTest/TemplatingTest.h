/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::TemplatingTest
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// This module defines the %TemplatingTest class that computes for each vertex of a
/// triangulation the average scalar value of itself and its direct neighbors.
///
/// \b Related \b publication: \n
/// 'TemplatingTest'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The TemplatingTest class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class TemplatingTest : virtual public Debug {

  public:
    TemplatingTest();

    /**
     * TODO 2: This method preconditions the triangulation for all operations
     *         the algorithm of this module requires. For instance,
     *         preconditionVertexNeighbors, preconditionBoundaryEdges, ...
     *
     *         Note: If the algorithm does not require a triangulation then
     *               this method can be deleted.
     */
    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    }

    /**
     * TODO 3: Implementation of the algorithm.
     *
     *         Note: If the algorithm requires a triangulation then this
     *               method must be called after the triangulation has been
     *               preconditioned for the upcoming operations.
     */
// template functions
template <class dataType, class triangulationType>
int execute(const triangulationType *triangulation,
                                     const int &numberOfIterations) const {
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
        // start a local timer for this subprocedure
        ttk::Timer localTimer;

        // print the progress of the current subprocedure (currently 0%)
        this->printMsg("Computing Averages",
                       0, // progress form 0-1
                       0, // elapsed time so far
                       this->threadNumber_, ttk::debug::LineMode::REPLACE);

        
        int nbCell = triangulation->getNumberOfCells();
        std::cout << "Il y a " << nbCell << " edges.\n";
        //auto cellEdges = triangulation->getCellEdges();
        //std::cout << "Il y a " << cellEdges.size() << "cell edges.\n";
        long long  sum = 0;
        for (int i = 0; i < numberOfIterations; i++)
        {
          for (int x = 0; x < nbCell; x++)
            sum += triangulation->getCellEdgeNumber(x);
        }

        std::cout << "sum =  " << sum << ".\n";

        // compute the average of each vertex in parallel
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

  }; // TemplatingTest class

} // namespace ttk
