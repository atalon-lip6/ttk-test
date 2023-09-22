/// \ingroup base
/// \class ttk::TopoMapHarishd
/// \author Alexandre Talon <alexandre.talon@lip6.fr>
/// \date March 2023
///
///
///TODO
/// This module defines the %TopoMapHarishd class that computes a
/// score indicating how good the low dimension distance matrix represents the
/// high dimension one. The score is computed according to the SIM formula.
///
///
/// \b Related \b publication: \n
/// 'Principal Geodesic Analysis of Merge Trees (and Persistence Diagrams)'
/// Mathieu Pont, Jules Vidal and Julien Tierny.
/// IEEE Transactions on Visualization and Computer Graphics, 2023.
///
/// \sa TopoMapHarishd

#pragma once

#define VERB 0
#define CHECK 1

// ttk common includes
#include <Debug.h>
#include <DimensionReduction.h>
#include <Geometry.h>
#include <UnionFind.h>
#include <geomutils.h>
#include <DisjointSets.h>

// STL includes
#include <map>
#include <set>
#include <vector>
#include <string>
// Boost includes
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>




namespace ttk {

  /**
   * The TopoMapHarishd class TODO
   */
  class TopoMapHarishd : virtual public Debug, public DimensionReduction {

    public:
      TopoMapHarishd();
      ~TopoMapHarishd();
      std::vector<Point> project(std::vector<double> &data, int dimension);

      /**
       * @brief Computes the topological mapper
       *
       * @param[out] outputCoords the final coordinates of the points, computed by the topological mapper
       *
       * @param[in] distMatrix the high-dimension input distance matrix
       * 
       * @return 0 in case of success.
       */
      template<typename T>
        int execute(T* outputCoords, std::vector<double> &inputCoords);
    protected:
      size_t AngleSamplingFreq{20};
      bool CheckMST{false};
      DisjointSets<int> comps;
      std::vector<Component> compMap;
      std::vector<Vertex> verts;

      size_t leafSize{1};
      bool verbose;

      void emst(std::vector<double> &data, int dimension, std::vector<std::pair<int, int>> &edges, std::vector<double> &weights);
      std::vector<Point> placePoints(const std::vector<std::pair<int, int>> &edges, const std::vector<double> &weights, int n);

      Component mergeComponents(Component &c1, Component &c2, int v1, int v2, double length);
      void transformComponent(const Component &c, const Transformation &t, double yOffset);
      Transformation alignHull(const Polygon &hull, const Point &p, bool topEdge);
      void log(std::string str);

  };





  //Sketch of the algorithm:
  //1. Sort the edges
  //2. For each edge uv with cost c_uv
  //     a. Get connected components comp(u), comp(v)
  //     b. Compute convex hull of comp(u) and of (comp(v))
  //     c. Find in hull(u) (resp. hull(v)) the point pU (resp. pV) in the convex hull which is closest to u (resp. v) in high dimension, u if possible.
  //     Let angleU (resp. angleV) the angle of hull(u) (resp. hull(v)) at vertex pU (resp. pV).
  //     d. Translate comp(u) such that pU lies on the line of the bissector of angleV, at distance c_uv of pV, and on the opposite direction of hull(v).
  //     e. Rotate comp(u) such that the bissector of angleU prolongates the one of angleU, and such that pU and pV are closest than any pair of points in comp(u),comp(v): the two components are facing each other at pU and pV.
  //     f. Try several rotations of comp(u) and comp(v) such that pU and pV still are the closest pair of vertices from the two components. We keep the angle of rotation which minimizes the difference between the high dimension distance matrix and the new distance matrix between the two components.
  //     g. comp(u) and comp(v) were merged, we now consider them as a single component.


  //TODO inputPoints remove
  template<typename T>
    int TopoMapHarishd::execute(T* outputCoords, std::vector<double> &inputCoords)
    {
      int n = inputCoords.size()/3;
      auto vect = project(inputCoords, 3);
      for (int i = 0; i < n; i++)
      {
        outputCoords[2*i] = vect[i].x;
        outputCoords[2*i+1] = vect[i].y;
      }
      return 0;
    }

} //ttk namespace
