#include <TopologicalMapper.h>
#include <UnionFind.h>
#include <vector>
#include <set>
#include <queue>
#include <map>
#include <utility>
#include <algorithm>
#include <Geometry.h>
#include "libqhullcpp/Qhull.h"

#include <cmath>

using namespace std;

float computeAngle(float const ptA[], float const ptB[], float const ptC[], bool is3D)
{
  float angle;
  if (is3D)
  {
    angle = ttk::Geometry::angle(ptA, ptB, ptA, ptC);
  }
  else
  {
    float ptA3D[3] = {ptA[0], ptA[1], 0}, ptB3D[3] = {ptB[0], ptB[1], 0}, ptC3D[3] = {ptC[0], ptC[1], 0};
    angle = ttk::Geometry::angle(ptA3D, ptB3D, ptA3D, ptC3D);
  }

  return angle;
}

void rotate(float ptToRotate[], float const centre[], float angle)
{
  //(cos(theta)+isin(theta))*(x+iy)
  //eitheta*(z-z0)+z0;
  const float &xCtr = centre[0], &yCtr = centre[1];
  float &xPt = ptToRotate[0], &yPt = ptToRotate[1];
  float dx = xPt-xCtr, dy = yPt-yCtr;
  xPt = dx*cos(angle)-dy*sin(angle)+xCtr;
  yPt = dx*sin(angle)+dy*cos(angle)+yCtr;
}


void printCoords(const char prefix[], float* coords)
{
  std::cout << prefix <<  coords[0] << "," << coords[1] << "\n";
}

inline float compute_dist2(const float ptA[], const float ptB[])
{
  float dx = ptB[0]-ptA[0], dy = ptB[1]-ptA[1];

  return dx*dx + dy*dy;
}

inline float compute_dist(const float ptA[], const float ptB[])
{
  float dx = ptB[0]-ptA[0], dy = ptB[1]-ptA[1];

  return sqrt(dx*dx + dy*dy);
}

inline void computeBarycenter(const std::vector<double> &coords, size_t dim, float baryCoords[])
{
  size_t nbPoint = coords.size()/dim;
  std::cout << " Il y a " << nbPoint << " pour le calcul de barycentre\n";
  for (size_t k = 0; k < dim; k++)
    baryCoords[k] = 0;
  for (size_t iPt = 0; iPt < nbPoint; iPt++)
  {
    for (size_t k = 0; k < dim; k++)
      baryCoords[k]+=coords[iPt*dim+k];
  }
  for (size_t k = 0; k < dim; k++)
    baryCoords[k] /= nbPoint;
}

inline void computeUnitVector(float* const coordOrig, float* const coordDest, float* const coordVect)
{
  float tmp[2] = {coordDest[0] - coordOrig[0], coordDest[1] - coordOrig[1]};
  float dist = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1]);
  if (dist < 1e-7)
  {
    coordVect[0] = 1;
    coordVect[1] = 0;
  }
  else
  {
    coordVect[0] = tmp[0]/dist;
    coordVect[1] = tmp[1]/dist;
  }
}

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


//int ttk::TopologicalMapper::execute(const std::vector<std::vector<float>> &distMatrix, unsigned int dim, const triangulationType &triangulation, std::vector<std::vector<float>> &outputCoords) const
int ttk::TopologicalMapper::execute(std::vector<double> &inputPoints, float* outputCoords, const std::vector<std::vector<float>> &distMatrix) const
{
#ifndef TTK_ENABLE_QHULL
  printErr("Error, qhull is not enabled. Please see the cmake configuration and enable it, and check that the package is installed on your system.");
  return 1;
#else
  printErr("Trobi1");
#endif

  float aa[2] = {0,0}, bb[2] = {1,0}, cc[2] = {0,1};
  float angleTest = computeAngle(aa, bb, cc, false);
  std::cout << angleTest << " was test angle :D\n";
  size_t n = distMatrix.size();
  size_t dim = (LowerDimension == ttk::TopologicalMapper::LOWER_DIMENSION::LOWER_DIM_2D ? 2 : 3);
    //TODO reverse order


  std::vector<std::pair<float, std::pair<size_t, size_t>>> edgeHeapVect;
  for (int u1 = 0; u1 < n; u1++)
  {
    for (int u2 = u1+1; u2 < n; u2++)
    {
      edgeHeapVect.push_back({distMatrix[u1][u2], {u1, u2}});

    }
  }
  sort(edgeHeapVect.begin(), edgeHeapVect.end());


  std::map<UnionFind*, std::set<size_t>> ufToSets;

  std::vector<UnionFind> ufVector(n);
  std::vector<UnionFind*> ufPtrVector(n);
  for (int i = 0; i < n; i++)
  {
    ufPtrVector[i] = &ufVector[i];
    ufToSets[ufPtrVector[i]].insert(i);
    //std::cerr << "item " << i << " ==> " << ufPtrVector[i] << std::endl;
  }
  for (auto& x : ufToSets)
  {
    for (auto &y : x.second)
      cout << y;
    cout << endl;
  }
  cout << " FINI LISTE \n";

  for (size_t i = 0; i < dim*n; i++)
    outputCoords[i] = 0;



  for (const auto &elt : edgeHeapVect)
  {
    cout << endl << endl;
    for (size_t iPt = 0; iPt < n; iPt++)
    {
      std::cout << "Coords of " << iPt << ": ";
      for (size_t k = 0; k < dim; k++)
      {
        std::cout << outputCoords[iPt*dim+k];
        if ((int)k < (int)dim-1)
          std::cout << ",";
      }
      std::cout << "\n";
    }
    std::cout << "....................................\n\n";

    float edgeCost = elt.first;
    size_t u = elt.second.first;
    size_t v = elt.second.second;

    std::cout << "considering edge " << u << "-" << v << " : " << edgeCost << endl;
    UnionFind *reprU = ufPtrVector[u]->find();
    UnionFind *reprV = ufPtrVector[v]->find();
    //std::cerr << "u,reprU = " << u << "," <<reprU << " AND v,reprV = " << v << "," << reprV << std::endl;
    if (reprU == reprV) // Already in the same component
    {
      continue;
    }
    std::set<size_t> &compU = ufToSets[reprU], &compV = ufToSets[reprV];
    std::cout << "Component1 = ";
    for (auto x : compU)
      std::cout << x << " ";

    std::cout << "\nComponent2 = ";
    for (auto x : compV)
      std::cout << x << " ";
    cout << "\n FIN DES COMPOSANTES\n";

    size_t idSmall = compU.size() < compV.size() ? u:v;
    size_t idBig = idSmall == u ? v:u;
    std::set<size_t> &compSmall = idSmall == u ? compU:compV;
    std::set<size_t> &compBig = idSmall == u ? compV:compU;
    size_t nBig = compBig.size();
    size_t nSmall = compSmall.size();

    size_t idEdgeVert[2] = {idSmall == u ? u:v, idSmall == u ? v:u};
    std::vector<double> pointsSets[2];
    std::vector<size_t> idsInHullSets[2];
    size_t idChosenVerts[2];

    for (int idSet = 0; idSet < 2; idSet++)
    {
      std::vector<double> &curPointsSet = pointsSets[idSet];// stores coordinates
      std::vector<size_t> &curHullVerts = idsInHullSets[idSet];
      size_t &idChosenVert = idChosenVerts[idSet];
      size_t idCur = idSet == 0 ? idSmall:idBig;
      std::set<size_t> &curComp = idCur == u ? compU:compV;
      size_t nCur = curComp.size();
      curPointsSet.resize(nCur*dim);
      std::vector<size_t> curCompVect(nCur);

      size_t cptCur = 0;
      for (int vertId : curComp)
      {
        for (int k = 0; k < dim; k++) //TODO seulement la dimension ?
        {
          curPointsSet[dim*cptCur+k] = outputCoords[dim*vertId+k];
          //curPointsSet[dim*cptCur+k] = inputPoints[dim*vertId+k];
        }
        curCompVect[cptCur] = vertId;
        cptCur++;
      }
      for (auto x : curPointsSet)
        cout << x << " ";
      cout << "\nFIN DES COORDS\n";

      getConvexHull(curPointsSet, dim, curHullVerts);
      // The ids in curHullVerts are the index of the vertices in the component list, not
      // the real ids of the vertices. The loop just below solves this.
      for (size_t &v : curHullVerts)
        v = curCompVect[v];
      idChosenVert = curHullVerts[0];
      float* coordEdgeVert = &outputCoords[idChosenVert];
      // We want to select, among all vertices in the convex hull, the one which is
      // closest to the vertex of the edge we work on.
      for (size_t vert : curHullVerts)
      {
        float* coordVert = &outputCoords[dim*vert];
        float dist2 = compute_dist2(coordVert, coordEdgeVert);
        if (dist2 < compute_dist2(&outputCoords[idChosenVert], coordEdgeVert))
          idChosenVert = vert;
      }

      cout << endl;
      std::cout << "\n\nConvex hull : ";
      for (size_t &x : idsInHullSets[idSet])
      {
        std::cout << x << " ";
      }
      std::cout << std::endl;

    } //TODO split into function ?

    std::vector<size_t> &idsInHullSmall = idsInHullSets[0], idsInHullBig = idsInHullSets[1];
    std::vector<double> &pointsBig = pointsSets[1], &pointsSmall = pointsSets[0];
    size_t sizeBigHull = idsInHullSets[1].size(), sizeSmallHull = idsInHullSets[0].size();//idsInHullBig.size();;
    std::vector<double> coordsBigHull(sizeBigHull*2), coordsSmallHull(sizeSmallHull*2);
    for (size_t iHull = 0; iHull < sizeBigHull; iHull++)
    {
      size_t vert = idsInHullBig[iHull];
      coordsBigHull[iHull*2] = outputCoords[vert*dim];
      coordsBigHull[iHull*2+1] = outputCoords[vert*dim+1];
    }
    for (size_t iHull = 0; iHull < sizeSmallHull; iHull++)
    {
      size_t vert = idsInHullSmall[iHull];
      coordsSmallHull[iHull*2] = outputCoords[vert*dim];
      coordsSmallHull[iHull*2+1] = outputCoords[vert*dim+1];
    }

    size_t idChosenBig = idChosenVerts[1], idChosenSmall = idChosenVerts[0];
    float coordsCentreBig[2], coordsCentreSmall[2];
    computeBarycenter(coordsBigHull, 2, coordsCentreBig);
    computeBarycenter(coordsSmallHull, 2, coordsCentreSmall);

    float unitCentreBigVect[2];
    computeUnitVector(&outputCoords[idChosenBig*dim], coordsCentreBig, unitCentreBigVect);
    printCoords("Destination for unit is: ", &outputCoords[idChosenBig*dim]);
    printCoords("Barycentre for unit is: ", coordsCentreBig);
        //{outputCoords[idChosenBig*dim] - coordsCentreBig[0], outputCoords[idChosenBig*dim+1] - coordsCentreBig[1]};

    float goalCoordChosenSmall[2] = {outputCoords[idChosenBig*dim]-edgeCost*unitCentreBigVect[0], outputCoords[idChosenBig*dim+1]-edgeCost*unitCentreBigVect[1]};
    cout << "sizeBig = " << sizeBigHull << " and " << compBig.size() << endl;
    if (sizeBigHull == 1)
    {
      goalCoordChosenSmall[0] = outputCoords[idChosenBig*dim]+edgeCost;
      goalCoordChosenSmall[1] = outputCoords[idChosenBig*dim+1];
    }
    float smallCompMoveVect[2] = {goalCoordChosenSmall[0] - outputCoords[idChosenSmall*dim], goalCoordChosenSmall[1] - outputCoords[idChosenSmall*dim+1]};

    float distBaryPointSmall = compute_dist(&outputCoords[idChosenSmall*dim], coordsCentreSmall);
    float finalPosBarySmall[2] = {goalCoordChosenSmall[0]+unitCentreBigVect[0]*distBaryPointSmall, goalCoordChosenSmall[1]+unitCentreBigVect[1]*distBaryPointSmall};


    float rotationAngle = computeAngle(goalCoordChosenSmall, finalPosBarySmall, coordsCentreSmall, dim==3);
    std::cout << " We chose for big and small: " << idChosenBig << " and " << idChosenSmall << std::endl;

    if (std::isfinite(rotationAngle))
    {
      std::cout << "angle is " << rotationAngle << std::endl;
      printCoords("coordCenter = ", goalCoordChosenSmall);
      printCoords("coordFrom = ", coordsCentreSmall);
      printCoords("coordTo = ", finalPosBarySmall);
    }
    for (size_t curIdSmall : compSmall)
    {//TODO pas bary ? osef il bougera pas !
      std::cout << "pointSmall : " << outputCoords[curIdSmall*dim] << "," << outputCoords[curIdSmall*dim+1] << endl;
      if (curIdSmall == idChosenSmall)
        cout << "\t\t=> CHOSEN small" << endl;
      outputCoords[curIdSmall*dim] += smallCompMoveVect[0];
      outputCoords[curIdSmall*dim+1] += smallCompMoveVect[1];
      if (std::isfinite(rotationAngle))
        rotate(&outputCoords[curIdSmall*dim], goalCoordChosenSmall/*finalPosBarySmall*/, rotationAngle);
    }
    std::cout << "Unit vector bary point = " << unitCentreBigVect[0] << "," << unitCentreBigVect[1] << endl;
    std::cout << "Goal coordinates = " << goalCoordChosenSmall[0] << "," << goalCoordChosenSmall[1] << endl;
    std::cout << "Moving vector coordinates = " << smallCompMoveVect[0] << "," << smallCompMoveVect[1] << endl;

    for (size_t curIdBig : compBig)
    {
      std::cout << "pointBig : " << outputCoords[curIdBig*dim] << "," << outputCoords[curIdBig*dim+1] << endl;
      if (curIdBig == idChosenBig)
        cout << "\t\t=> CHOSEN big" << endl;
    }

    cout << "barycenterBig = " << coordsCentreBig[0] << "," << coordsCentreBig[1] << endl;

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


    cout << " Resulting component = ";
    for (auto x : unionSet)
      cout << x << " ";
    cout << "=======================\n";
    //continue;

    std::cerr << "\tDELETING " << otherRepr << "( and was " <<unionRepr << ")"<< std::endl;
    ufToSets.erase(otherRepr);
    //ufToSets[otherRepr] = unionSet;
  }

  return 0;
}

void ttk::TopologicalMapper::rotatePolygon(std::vector<float> &coords, size_t dim, float* centerCoords, const float angle) const
{
  float xCenter = centerCoords[0], yCenter = centerCoords[1];
  size_t nbPoint = coords.size()/dim;
  std::cout << "center = " << xCenter << "," << yCenter << std::endl;
  for (size_t iPt = 0; iPt < nbPoint; iPt++)
  {
    float &x = coords[iPt*dim], &y = coords[iPt*dim+1];
    float xNew, yNew;
    std::cout << "point = " << x << "," << y << std::endl;
    xNew = (x-xCenter)*cos(angle)-(y-yCenter)*sin(angle)+xCenter;
    yNew = (y-yCenter)*cos(angle)+(x-xCenter)*sin(angle)+yCenter;
    x = xNew;
    y = yNew;
  }
}


// Tests if all points have the same coordinates or are aligned. If so, the convex hull
// is computed by this code. Otherwise, we call qHull if the points are fully 2D.
void ttk::TopologicalMapper::getConvexHull(const std::vector<double>& coords, size_t dim, std::vector<size_t> &idsInHull) const
{
  size_t nbPoint = coords.size()/dim;
  if (nbPoint <= 2)
  {
    idsInHull.push_back(0);
    if (nbPoint == 2)
    {
      double dist2 = (coords[0]-coords[2])*(coords[0]-coords[2]) + (coords[1]-coords[3])*(coords[1]-coords[3]);

      if (dist2 > 1e-7)
        idsInHull.push_back(1);
    }
    return;
  }

  // Testing if all points are colinear

  double dirVect[2] = {0,0};
  bool areColinear = true;

  size_t idFirstDistinct = 1;
  while (idFirstDistinct < nbPoint && abs(dirVect[0]) < 1e-7 && abs(dirVect[1]) < 1e-7)
  {
    dirVect[0] = coords[2*idFirstDistinct]-coords[0];
    dirVect[1] = coords[2*idFirstDistinct+1]-coords[1]; //TODO tester avec deux points mêmes coordonnées !
    if (abs(dirVect[0]) < 1e-7)
      dirVect[0] =  0;
    if (abs(dirVect[1]) < 1e-7)
      dirVect[1] =  0;

    idFirstDistinct++;
  }

  double idMins[2] = {coords[0] < coords[2] ? 0:1, coords[1] < coords[3] ? 0:1};
  double idMaxs[2] = {coords[0] > coords[2] ? 0:1, coords[1] > coords[3] ? 0:1};

  for (size_t iPt = idFirstDistinct; iPt < nbPoint; iPt++)
  {
    double curVect[2] = {coords[2*iPt]-coords[0], coords[2*iPt+1]-coords[1]};
    if (abs(curVect[0]) < 1e-7 && abs(curVect[1]) < 1e-7)
      continue;
    if (abs(curVect[0]*dirVect[1] - curVect[1]*dirVect[0]) > 1e-7)
    {
      areColinear = false;
      break;
    }

    if (coords[2*iPt] < coords[2*idMins[0]])
      idMins[0] = iPt;
    if (coords[2*iPt+1] < coords[2*idMins[1]])
      idMins[1] = iPt;
    if (coords[2*iPt] > coords[2*idMaxs[0]])
      idMaxs[0] = iPt;
    if (coords[2*iPt+1] > coords[2*idMaxs[1]])
      idMaxs[1] = iPt;

  }

  if (areColinear)
  {
    std::cout << "COLINEAR = " << idMins[0] << "," << idMins[1] << "  ;  " << idMaxs[0] << "," << idMaxs[1] << std::endl;
    if (abs(dirVect[0]) > 1e-7)
    {
      idsInHull.push_back(idMins[0]);
      idsInHull.push_back(idMaxs[0]);
    }
    else
    {
      idsInHull.push_back(idMins[1]);
      idsInHull.push_back(idMaxs[1]);
    }
    return;
  }

  char qHullFooStr[2] = "Q";
  try
  {
    orgQhull::Qhull qhullCur;
    qhullCur.runQhull(qHullFooStr, dim, nbPoint, coords.data(), qHullFooStr);
    for (auto u : qhullCur.vertexList())
      idsInHull.push_back(u.id()-1);

    /*
    for (orgQhull::Qhull &qhull : {qhullBig, qhullSmall})
    {
      std::cout << "qhullBig.hullDimension(): " << qhullCur.hullDimension() << "\n";
      std::cout << "qhullBig.volume(): " << qhullCur.volume() << "\n";
      std::cout << "qhullBig.area(): " << qhullCur.area() << "\n";
    }
    */


  }
  catch (orgQhull::QhullError &e)
  {
    printErr("Error with qHull module: " + std::string(e.what()));
    //return -1;
  }
}


