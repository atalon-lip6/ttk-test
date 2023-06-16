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
  std::vector<UnionFind*> ufPtrVector(n);
  for (int i = 0; i < n; i++)
  {
    ufPtrVector[i] = &ufVector[i];
    ufToSets[ufPtrVector[i]].insert(i);
  }


  while (!edgeHeap.empty())
  {
    const auto &elt = edgeHeap.top();
    edgeHeap.pop();
    float edgeCost = elt.first;
    size_t u = elt.second.first;
    size_t v = elt.second.second;

    UnionFind *reprU = ufPtrVector[u]->find();
    UnionFind *reprV = ufPtrVector[v]->find();
    if (reprU == reprV) // Already in the same component
      continue;
    std::set<size_t> &compU = ufToSets[reprU], &compV = ufToSets[reprV];
    size_t idSmall = compU.size() < compV.size() ? u:v;
    size_t idBig = idSmall == u ? v:u;
    std::set<size_t> &compSmall = idSmall == u ? compU:compV;
    std::set<size_t> &compBig = idSmall == u ? compV:compU;
    size_t nBig = compBig.size();
    size_t nSmall = compSmall.size();


    if (nBig <= 2) // We just put the point at the right of the convex hull
    {
      float xMax = -DBL_MAX; //TODO include geometry?
      std::vector<double> pointsBig(nDim*(std::max((size_t)3,nBig)));
      std::cout << "nb pts = " << nBig << endl;
      int vertIdInSmallSet = -1;
      int cpt = 0;
      for (int vertId : compBig)
      {
        if (vertId == idBig)
          vertIdInSmallSet = cpt; //TODO in small or in big, utilisé ?
        for (int k = 0; k < nDim; k++) //TODO seulement la dimension ?
          pointsBig[nDim*cpt+k] = outputCoords[nDim*vertId+k];
        cpt++;
      }
      if (compBig.size() <= 2)
      {
        pointsBig[3] = 0.34;
        pointsBig[4] = 0.13;
        pointsBig[5] = 0.07;

        if (compBig.size() == 1)
        {
          pointsBig[6] = 0.5;
          pointsBig[7] = 0.34;
          pointsBig[8] = 0.23;
        }
      }
      char qHullFooStr[1] = "";

      orgQhull::Qhull qhullBig;
      try
      {
        qhullBig.runQhull(qHullFooStr, nDim, std::max((size_t)3, nBig), pointsBig.data(), qHullFooStr);

        {
          std::cout << "qhullBig.hullDimension(): " << qhullBig.hullDimension() << "\n";
          std::cout << "qhullBig.volume(): " << qhullBig.volume() << "\n";
          std::cout << "qhullBig.area(): " << qhullBig.area() << "\n";
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
        if (vert.id() == vertIdInSmallSet)
        {
          chosenVert = idBig;
          break;
        }
      }


      for (size_t eltBig : compBig) //TODO tester chaque pt enveloppe
        xMax = std::max(xMax, outputCoords[nDim*eltBig]);

      outputCoords[nDim*idSmall] = xMax+edgeCost;
    }

    else // Each component has at least two points
    {
      std::vector<double> pointsBig, pointsSmall;
      pointsBig.reserve(nBig*nDim); //TODO reserve before looop... scope
      pointsSmall.reserve(nSmall*nDim);
      char qHullFooStr[1] = "";

      size_t cptU = 0;
      for (int vertId : compU)
      {
        for (int k = 0; k < nDim; k++) //TODO seulement la dimension ?
          pointsBig[nDim*cptU+k] = outputCoords[nDim*vertId+k];
        cptU++;

      }
      size_t cptV = 0;
      for (int vertId : compV)
      {
        for (int k = 0; k < nDim; k++) //TODO seulement la dimension ?
          pointsSmall[nDim*cptV+k] = outputCoords[nDim*vertId+k];
        cptV++;

      }

      /*for (size_t idU : compU)
        for (float& coord : outputCoords[idU])
          pointsBig.reserve(coord);
      for (size_t idV : compV)
        for (float& coord : outputCoords[idV])
          pointsSmall.reserve(coord);
      */
      try
      {
        orgQhull::Qhull qhullBig, qhullSmall;
        std::cout << "nb ptsU, nbPtsV = " << nBig << "," << nSmall << endl;
        qhullBig.runQhull(qHullFooStr, nDim, nBig, pointsBig.data(), qHullFooStr);
        qhullSmall.runQhull(qHullFooStr, nDim, nSmall, pointsSmall.data(), qHullFooStr);

        //for (orgQhull::Qhull &qhull : {qhullBig, qhullSmall})
        {
          std::cout << "qhullBig.hullDimension(): " << qhullBig.hullDimension() << "\n";
          std::cout << "qhullBig.volume(): " << qhullBig.volume() << "\n";
          std::cout << "qhullBig.area(): " << qhullBig.area() << "\n";
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

    if (otherRepr == reprU)
      ufPtrVector[u] = unionRepr;
    else
      ufPtrVector[v] = unionRepr;
    unionSet.insert(otherSet.begin(), otherSet.end());
    //TODO faire des trucs;

    ufPtrVector[u] = unionRepr; //TODO plutôt inutile
    ufPtrVector[v] = unionRepr;

    // We change the coordinates
    float minXUnion = 1e50, maxXOther = -1e50;
    for (size_t id : unionSet)
      minXUnion = std::min(minXUnion, outputCoords[nDim*id]);

    for (size_t id : otherSet)
      minXUnion = std::max(maxXOther, outputCoords[nDim*id]);

    float shift = maxXOther + edgeCost - minXUnion;
    for (size_t id : unionSet)
      outputCoords[nDim*id] += shift;
    ufToSets.erase(otherRepr); // TODO vérifier que ça clear le set;
    //ufToSets[otherRepr] = unionSet;
  }

  return 0;
}
