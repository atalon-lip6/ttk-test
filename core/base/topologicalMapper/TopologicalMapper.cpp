#include <TopologicalMapper.h>
#include <UnionFind.h>
#include <vector>
#include <set>
#include <queue>
#include <map>
#include <utility>
#include <algorithm>
#include "libqhullcpp/Qhull.h"


ttk::TopologicalMapper::TopologicalMapper() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("TopologicalMapper");
}
ttk::TopologicalMapper::~TopologicalMapper() = default;
  //1. Sort the edges
  //2. For each edge uv with cost c_uv

    // Get connected components comp(u), comp(v)

    // Compute convex hull of comp(u) and of (comp(v))


    // Find longest edge (between two points) e_u of comp(u), same for e_v

    // Rotate the two connected components so that they are paraller to the x axis

    // Put c_uv distance between the two components


//int ttk::TopologicalMapper::execute(const std::vector<std::vector<float>> &distMatrix, unsigned int nDim, const triangulationType &triangulation, std::vector<std::vector<float>> &outputCoords) const
int ttk::TopologicalMapper::execute(float* outputCoords, const std::vector<std::vector<float>> &distMatrix) const
{
#ifndef TTK_ENABLE_QHULL
  printErr("Error, qhull is not enabled. Please see the cmake configuration and enable it, and check that the package is installed on your system.");
  return 1;
#else
  printErr("Trobi1");
#endif

  size_t n = distMatrix.size();
  size_t nDim = (LowerDimension == ttk::TopologicalMapper::LOWER_DIMENSION::LOWER_DIM_2D ? 2 : 3);
  //outputCoords.resize(n);
  //for (size_t i = 0; i < n; i++)
  //  outputCoords[i].resize(nDim, 0);

  std::priority_queue<std::pair<float, std::pair<size_t, size_t>>> edgeHeap; // We store elements as (c_uv,(u,v))
                                                                //TODO reverse order
  for (int u1 = 0; u1 < n; u1++)
  {
    for (int u2 = u1+1; u2 < n; u2++)
    {
      edgeHeap.push({distMatrix[u1][u2], {u1,u2}});

    }
  }


  std::map<UnionFind*, std::set<size_t>> ufToSets;

  std::vector<UnionFind> ufVector(n);
  for (int i = 0; i < n; i++)
    ufToSets[&ufVector[i]].insert(i);


  while (!edgeHeap.empty())
  {
    const auto &elt = edgeHeap.top();
    edgeHeap.pop();
    float edgeCost = elt.first;
    size_t u = elt.second.first;
    size_t v = elt.second.second;

    UnionFind *reprU = ufVector[u].find();
    UnionFind *reprV = ufVector[v].find();
    if (reprU == reprV) // Already in the same component
      continue;

    std::set<size_t> &compU = ufToSets[reprU], &compV = ufToSets[reprV];

    if (compU.size() <= 2 || compV.size() <= 2) // We just put the point at the right of the convex hull
    {
      size_t idSingle = compU.size() <= 2 ? u : v;
      size_t idOther = idSingle ^ u ^ v;
      std::set<size_t> &compBig = (idSingle == u ? compV : compU);
      size_t nBig = std::max(compU.size(), compV.size());
      std::set<size_t> bigComp = (idSingle == u) ? compV : compU;
      float xMax = -DBL_MAX; //TODO include geometry?
      std::vector<double> pointsBig((std::max(3,nBig)*nDim));
      std::cout << "nb pts = " << nBig << endl;
      int indexInSetOtherVert = -1;
      int cpt = 0;
      for (int vertId : compBig)
      {
        if (vertId == idOther)
          indexInSetOtherVert = cpt;
        for (int k = 0; k < 3; k++) //TODO seulement la dimension ?
          pointsBig[3*cpt+k] = outputCoords[3*vertId+k];
        cpt++;
      }
      if (compBig.size() <= 2)
      {
        pointsBig[3] = 0.34;
        pointsBig[4] = 0;
        pointsBig[5] = 0;

        if (compBig.size() == 1)
        {
          pointsBig[6] = 0;
          pointsBig[7] = 0.34;
          pointsBig[8] = 0;
        }
      }
      char qHullFooStr[1] = "";

      orgQhull::Qhull qhullBig;
      try
      {
        qhullBig.runQhull(qHullFooStr, nDim, nBig, pointsBig.data(), qHullFooStr);

        {
          std::cout << "qhullU.hullDimension(): " << qhullBig.hullDimension() << "\n";
          std::cout << "qhullU.volume(): " << qhullBig.volume() << "\n";
          std::cout << "qhullU.area(): " << qhullBig.area() << "\n";
        }
      }
      catch (orgQhull::QhullError &e)
      {
        printErr("Error with qHull module: " + std::string(e.what()));
        return -1;
      }
      size_t chosenVert = *compBig.begin();

      for(const orgQhull::QhullVertex &vert: qhullBig.vertexList())
      {
        if (vert.id() == indexInSetOtherVert)
        {
          chosenVert = idOther;
          break;
        }
      }


      for (size_t eltBig : bigComp) //TODO tester chaque pt enveloppe
        xMax = std::max(xMax, outputCoords[3*eltBig]);

      outputCoords[3*idSingle] = xMax+edgeCost;
    }

    else // Each component has at least two points
    {
      size_t nU = compU.size(), nV = compV.size();
      std::vector<double> pointsU, pointsV;
      pointsU.reserve(nU*nDim); //TODO reserve before looop... scope
      pointsV.reserve(nV*nDim);
      char qHullFooStr[1] = "";

      size_t cptU = 0;
      for (int vertId : compU)
      {
        for (int k = 0; k < 3; k++) //TODO seulement la dimension ?
          pointsU[3*cptU+k] = outputCoords[3*vertId+k];
        cptU++;

      }
      size_t cptV = 0;
      for (int vertId : compV)
      {
        for (int k = 0; k < 3; k++) //TODO seulement la dimension ?
          pointsV[3*cptV+k] = outputCoords[3*vertId+k];
        cptV++;

      }

      /*for (size_t idU : compU)
        for (float& coord : outputCoords[idU])
          pointsU.reserve(coord);
      for (size_t idV : compV)
        for (float& coord : outputCoords[idV])
          pointsV.reserve(coord);
      */
      try
      {
        orgQhull::Qhull qhullU, qhullV;
      std::cout << "nb ptsU, nbPtsV = " << nU << "," << nV << endl;
        qhullU.runQhull(qHullFooStr, nDim, nU, pointsU.data(), qHullFooStr);
        qhullV.runQhull(qHullFooStr, nDim, nV, pointsV.data(), qHullFooStr);

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


    // We change the coordinates
    float minXUnion = 1e50, maxXOther = -1e50;
    for (size_t id : unionSet)
      minXUnion = std::min(minXUnion, outputCoords[3*id]);

    for (size_t id : otherSet)
      minXUnion = std::max(maxXOther, outputCoords[3*id]);

    float shift = maxXOther + edgeCost - minXUnion;
    for (size_t id : unionSet)
      outputCoords[3*id] += shift;
  }

  return 0;
}
