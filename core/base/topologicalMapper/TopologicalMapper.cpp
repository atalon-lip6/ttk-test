#include <TopologicalMapper.h>
#include <UnionFind.h>
#include <vector>
#include <set>
#include <queue>
#include <map>
#include <algorithm>
#include "libqhullcpp/Qhull.h"


ttk::TopologicalMapper::TopologicalMapper() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("TopologicalMapper");
}

  //1. Sort the edges
  //2. For each edge uv with cost c_uv

  //TODO 1 'n' ou 2 'n' ? :p
    // Get conected components comp(u), comp(v)

    // Compute convex hull of comp(u) and of (comp(v))


    // Find longest edge (between two points) e_u of comp(u), same for e_v

    // Rotate the two conected components so that they are paraller to the x axis

    // Put c_uv distance between the two components



int ttk::TopologicalMapper::execute(const std::vector<std::vector<double>> distMatrix, unsigned int nDim, std::vector<std::vector<double>> &outputCoords)
{
#ifndef TTK_ENABLE_QHULL
  printErr("Error, qhull is not enabled. Please see the cmake configuration and enable it, and check that the package is installed on your system.");
  return 1;
#else
  printErr("Trobi1");
#endif

  size_t n = distMatrix.size();
  outputCoords.resize(n);
  for (size_t i = 0; i < n; i++)
    outputCoords[i].resize(nDim, 0);

  std::priority_queue<std::pair<double, std::pair<size_t, size_t>>> edgeHeap; // We store elements as (c_uv,(u,v))
                                                                //TODO reverse order

  std::map<UnionFind*, std::set<size_t>> ufToSets;

  std::vector<UnionFind> ufVector(n);


  while (!edgeHeap.empty())
  {
    const auto &elt = edgeHeap.top();
    edgeHeap.pop();
    double edgeCost = elt.first;
    size_t u = elt.second.first;
    size_t v = elt.second.second;

    UnionFind *reprU = ufVector[u].find();
    UnionFind *reprV = ufVector[v].find();
    if (reprU == reprV) // Already in the same component
      continue;

    std::set<size_t> &compU = ufToSets[reprU], &compV = ufToSets[reprV];

    if (compU.size() == 1 || compV.size() == 1) // We just put the point at the right of the convex hull
    {
      size_t idSingle = compU.size() == 1 ? u : v;
      std::set<size_t> bigComp = (idSingle == u) ? compV : compU;
      double xMax = -DBL_MAX; //TODO include geometry?
      for (size_t eltBig : bigComp)
        xMax = std::max(xMax, outputCoords[eltBig][0]);

      outputCoords[idSingle][0] = xMax+edgeCost;
    }

    else // Each component has at least two points
    {
      size_t nU = compU.size(), nV = compV.size();
      std::vector<double> pointsU, pointsV;
      pointsU.reserve(nU*nDim);
      pointsV.reserve(nV*nDim);
      char qHullFooStr[1] = "";
        for (size_t idU : compU)
          for (double& coord : outputCoords[idU])
            pointsU.reserve(coord);
        for (size_t idV : compV)
          for (double& coord : outputCoords[idV])
            pointsV.reserve(coord);

      try
      {
        orgQhull::Qhull qhullU, qhullV;
        qhullU.runQhull(qHullFooStr, nDim, nU, pointsU.data(), qHullFooStr);
        qhullV.runQhull(qHullFooStr, nDim, nU, pointsU.data(), qHullFooStr);

        //for (orgQhull::Qhull &qhull : {qhullU, qhullV})
        {
          std::cout << "qhullU.hullDimension(): " << qhullU.hullDimension() << "\n";
          std::cout << "qhullU.volume(): " << qhullU.volume() << "\n";
          std::cout << "qhullU.area(): " << qhullU.area() << "\n";
        }


      }
      catch (orgQhull::QhullError &e)
      {
        printErr("Error with qHull module: " + std::string(e.what()));
        return -1;
      }

    }

    UnionFind* unionRepr = UnionFind::makeUnion(reprU, reprV);
    UnionFind* otherRepr = (unionRepr == reprU) ? reprV : reprU;
    std::set<size_t> &unionSet = ufToSets[unionRepr];
    std::set<size_t> &otherSet = ufToSets[otherRepr];

    unionSet.insert(otherSet.begin(), otherSet.end());
    ufToSets.erase(otherRepr); // TODO vérifier que ça clear le set;
    //TODO faire des trucs;

  }

  return 0;
}
