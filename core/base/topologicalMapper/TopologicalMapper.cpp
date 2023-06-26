#include <TopologicalMapper.h>
#include <UnionFind.h>
#include <vector>
#include <set>
#include <queue>
#include <map>
#include <utility>
#include <algorithm>
#include "libqhullcpp/Qhull.h"

using namespace std;

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
    //std::cerr << "item " << i << " ==> " << ufPtrVector[i] << std::endl;
  }

  for (size_t i = 0; i < nDim*n; i++)
    outputCoords[i] = 0;


  while (!edgeHeap.empty())
  {
    const auto &elt = edgeHeap.top();
    edgeHeap.pop();
    float edgeCost = elt.first;
    size_t u = elt.second.first;
    size_t v = elt.second.second;

    UnionFind *reprU = ufPtrVector[u]->find();
    UnionFind *reprV = ufPtrVector[v]->find();
    //ufPtrVector[u] = reprU;
    //ufPtrVector[v] = reprV;
    //std::cerr << "u,reprU = " << u << "," <<reprU << " AND v,reprV = " << v << "," << reprV << std::endl;
    if (reprU == reprV) // Already in the same component
    {
      //std::cerr << "same component\n";
      continue;
    }
    std::set<size_t> &compU = ufToSets[reprU], &compV = ufToSets[reprV];
    size_t idSmall = compU.size() < compV.size() ? u:v;
    size_t idBig = idSmall == u ? v:u;
    std::set<size_t> &compSmall = idSmall == u ? compU:compV;
    std::set<size_t> &compBig = idSmall == u ? compV:compU;
    size_t nBig = compBig.size();
    size_t nSmall = compSmall.size();

    /*
    if (nBig <= 2) // We just put the point at the right of the convex hull
    {
      cout << "nBIg = " << nBig << "    and     nSMall = " << nSmall << endl;
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
        {
          cout << outputCoords[nDim*vertId+k] << " ";
          pointsBig[nDim*cpt+k] = outputCoords[nDim*vertId+k];
        }
        cout << "      WAS " << vertId << endl;
        cpt++;
      }
      if (compBig.size() <= 2)
      {
        cout << "DOING STUFF\n";
        pointsBig[nDim] = 0.34;
        pointsBig[nDim+1] = 0.13;
        if (nDim == 3)
          pointsBig[nDim+2] = 0.07;

        if (compBig.size() == 1)
        {
          pointsBig[2*nDim] = 0.5;
          pointsBig[2*nDim+1] = 0.34;
          if (nDim == 3)
            pointsBig[2*nDim+2] = 0.23;
        }
      }

      char qHullFooStr[3] = "QJ";

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
        cout << endl;
        for (int i = 0; i < std::max((size_t)3, nBig); i++)
        {
          for (int k = 0; k < nDim; k++)
            cout << pointsBig[3*i+k] << " ";
          cout << endl;
        }
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
    */

    std::vector<double> pointsSets[2];
    std::vector<double> &pointsBig = pointsSets[1], &pointsSmall = pointsSets[0];
    std::vector<int> idsInHullSets[2];
    std::vector<int> &idsInHullSmall = idsInHullSets[0], idsInHullBig = idsInHullSets[1];


    for (int idSet = 0; idSet < 2; idSet++)
    {
      std::vector<double> &curPointsSet = pointsSets[idSet];
      size_t idCur = idSet == 0 ? idSmall:idBig;
      std::set<size_t> &curComp = idCur == u ? compU:compV;
      size_t nCur = curComp.size();
      curPointsSet.resize(nCur*nDim);

      size_t cptCur = 0;
      for (int vertId : curComp)
      {
        for (int k = 0; k < nDim; k++) //TODO seulement la dimension ?
          curPointsSet[nDim*cptCur+k] = outputCoords[nDim*vertId+k];
        cptCur++;
      }

      char qHullFooStr[1] = "";
      try
      {
        orgQhull::Qhull qhullCur;
        qhullCur.runQhull(qHullFooStr, nDim, nCur, curPointsSet.data(), qHullFooStr);

        //for (orgQhull::Qhull &qhull : {qhullBig, qhullSmall})
        {
          std::cout << "qhullBig.hullDimension(): " << qhullCur.hullDimension() << "\n";
          std::cout << "qhullBig.volume(): " << qhullCur.volume() << "\n";
          std::cout << "qhullBig.area(): " << qhullCur.area() << "\n";
        }


      }
      catch (orgQhull::QhullError &e)
      {
        printErr("Error with qHull module: " + std::string(e.what()));
        return -1;
      }



    }


    //else // Each component has at least two points
    {
      std::vector<double> pointsBig, pointsSmall;
      pointsBig.resize(nBig*nDim); //TODO reserve before looop... scope
      pointsSmall.resize(nSmall*nDim);
      char qHullFooStr[3] = "QJ";

      size_t cptBig = 0;
      for (int vertId : compBig)
      {
        for (int k = 0; k < nDim; k++) //TODO seulement la dimension ?
          pointsBig[nDim*cptBig+k] = outputCoords[nDim*vertId+k];
        cptBig++;

      }
      size_t cptSMall = 0;
      for (int vertId : compSmall)
      {
        for (int k = 0; k < nDim; k++) //TODO seulement la dimension ?
          pointsSmall[nDim*cptSMall+k] = outputCoords[nDim*vertId+k];
        cptSMall++;

      }

      /* TODO voir ça !!!!
      for (size_t idBig : compBig)
        for (float& coord : outputCoords[idBig])
          pointsBig.reserve(coord);
      for (size_t idSmall : compSmall)
        for (float& coord : outputCoords[idSmall])
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
    if (unionRepr != reprU && unionRepr != reprV)
      std::cerr << "NOOOOOOOOOOOOOOO\n";
    std::set<size_t> &unionSet = ufToSets[unionRepr];
    std::set<size_t> &otherSet = ufToSets[otherRepr];

    unionSet.insert(otherSet.begin(), otherSet.end());
    //TODO faire des trucs;

    ufPtrVector[u] = unionRepr;
    ufPtrVector[v] = unionRepr;

    // We change the coordinates
    float minXUnion = 1e50, maxXOther = -1e50;
    for (size_t id : unionSet)
      minXUnion = std::min(minXUnion, outputCoords[nDim*id]);

    //auto it = otherSet.begin(;
    for (size_t id : otherSet)
      maxXOther = std::max(maxXOther, outputCoords[nDim*id]);

    float shift = maxXOther + edgeCost - minXUnion;
    for (size_t id : unionSet)
      outputCoords[nDim*id] += shift;
    std::cerr << "\tDELETING " << otherRepr << "( and was " <<unionRepr << ")"<< std::endl;
    ufToSets.erase(otherRepr); // TODO vérifier que ça clear le set;
    //ufToSets[otherRepr] = unionSet;
  }

  return 0;
}


void getConvexHull(const std::vector<float>& coords, size_t dim, std::vector<size_t> &idsInHull)
{
  //Test if all points the same
  //Test if all points aligned
  //Test if all points coplanar

  bool runQhull = true;
  size_t idMins[3], idMaxs[3];
  size_t n = coords.size()/dim;
  for (size_t i = 1; i < n; i++)
  {
    size_t iPt = nDim*i;
    for (size_t k = 0; k < dim; k++)
    {
      if (coords[iPt+k] < coords[dim*idMins[k]+k])
        idMins[k] = iPt;
      if (coords[iPt+k] > coords[dim*idMaxs[k]+k])
        idMaxs[k] = iPt;

    }
    if (coords[iPt] < coords[dim*idMinX])
      idMinX = iPt;
    if (coords[iPt] > coords[



  }
  if (dim == 3)
  {


  }



}
