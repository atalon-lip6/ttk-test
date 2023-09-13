/// \ingroup base
/// \class ttk::TopologicalMapper
/// \author Alexandre Talon <alexandre.talon@lip6.fr>
/// \date March 2023
///
///
///TODO
/// This module defines the %TopologicalMapper class that computes a
/// score indicating how good the low dimension distance matrix represents the
/// high dimension one. The score is computed according to the SIM formula.
///
///
/// \b Related \b publication: \n
/// 'Principal Geodesic Analysis of Merge Trees (and Persistence Diagrams)'
/// Mathieu Pont, Jules Vidal and Julien Tierny.
/// IEEE Transactions on Visualization and Computer Graphics, 2023.
///
/// \sa TopologicalMapper

#pragma once

#define VERB 0
#define CHECK 1

#include <vector>
// ttk common includes
#include <Debug.h>
#include <DimensionReduction.h>
#include <UnionFind.h>
#include "libqhullcpp/Qhull.h"
#include <map>
#include <set>

using namespace std;



inline double compute_dist2(const double ptA[], const double ptB[])
{
  double dx = ptB[0]-ptA[0], dy = ptB[1]-ptA[1];

  return dx*dx + dy*dy;
}

inline double compute_dist(const double ptA[], const double ptB[])
{
  double dx = ptB[0]-ptA[0], dy = ptB[1]-ptA[1];

  return sqrt(dx*dx + dy*dy);
}

inline bool are_colinear(double const ptA[], double const ptB[], double const ptC[])
{
  return (abs(ptA[0]*(ptB[1]-ptC[1])+ptB[0]*(ptC[1]-ptA[1])+ptC[0]*(ptA[1]-ptB[1]))) <= 1e-14;
}

void computeUnitVector(double* const coordOrig, double* const coordDest, double* const coordVect);

double computeAngle(double const ptA[], double const ptB[], double const ptC[]);

void rotatePolygon(std::vector<double> &coords, double* centerCoords, const double angle);

void rotate(double ptToRotate[], double const centre[], double angle);

void getPrevPostEdges(const std::vector<size_t> &idsPtsPolygon, size_t idCenter, double* const allCoords, double coordPrev[2], double coordPost[2]);

void rotateMergingCompsBest(const std::vector<size_t> &hull1, const std::vector<size_t> &hull2, const std::set<size_t> &comp1, const std::set<size_t> &comp2, size_t iPt1, size_t iPt2, const std::vector<std::vector<float>> &distMatrix, double* allCoords, size_t angleSamplingFreq, size_t nThread);


void getConvexHull(const std::vector<double>& coords, size_t dim, std::vector<size_t> &idsInHull);

namespace ttk {

  /**
   * The TopologicalMapper class TODO
   */
  class TopologicalMapper : virtual public Debug, public DimensionReduction {

  public:
    TopologicalMapper();
    ~TopologicalMapper();

template<typename T>
int execute(T* outputCoords, const std::vector<std::vector<T>> &distMatrix) const;

  protected:
    bool isQhullEnabled(void) const;

    size_t AngleSamplingFreq{20};
  };



//1. Sort the edges
//2. For each edge uv with cost c_uv

// Get connected components comp(u), comp(v)

// Compute convex hull of comp(u) and of (comp(v))


// Find shortest edge (between two points) e_u of comp(u), same for e_v

// Rotate the two connected components so that they are parallel to the x axis

// Put c_uv distance between the two components




//TODO inputPoints remove
template<typename T>
int ttk::TopologicalMapper::execute(T* outputCoords, const std::vector<std::vector<T>> &distMatrix) const
{
  cout << DBL_DIG << " trop cool\n";
  std::vector<double> edgesMSTBefore, edgesMSTAfter;
  if (!isQhullEnabled())
  {
    printErr("Error, qhull is not enabled. Please see the cmake configuration and enable it, and check that the package is installed on your system.");
  return 1;
  }

  // Only small test
  double aa[2] = {24,27}, bb[2] = {21,25}, cc[2] = {20,30};
  double angleTest = computeAngle(aa, bb, cc);
  //std::cout << "avant " << deg(angleTest) << " (" << angleTest << ") was test angle :D\n";

  double aabis[2] = {-0.707,-2.121}, bbbis[2] = {2.414,1}, ccbis[2] = {3.414,-2};
  double angleTestbis = computeAngle(aabis, bbbis, ccbis);
  //std::cout << "après " << deg(angleTestbis) << " (" << angleTestbis << ") was test angle :D\n";


  size_t n = distMatrix.size();
  cout << "Il y a " << n << " points! :-)\n";
  size_t dim = 2;

  std::vector<double> outputCoordsTmp;
  if (std::is_same<T, float>())
  {
    outputCoordsTmp.resize(2*n);
    for (int i = 0; i < n*2; i++)
      outputCoordsTmp[i] = outputCoords[i];
  }
  double* outputCoordsPtr = std::is_same<T, double>() ? (double*)outputCoords :  outputCoordsTmp.data();
  std::vector<std::pair<double, std::pair<size_t, size_t>>> edgeHeapVect(n*(n+1)/2);//TODO identifiant i,j
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
  }

  for (size_t i = 0; i < 2*n; i++)
    outputCoordsPtr[i] = 0;


  size_t nbEdgesMerged = 0;
  for (const auto &elt : edgeHeapVect)
  {
    if (nbEdgesMerged == n-1)
      break;
    //if (stop)
    //  break;
    //return 0;
    double edgeCost = elt.first;
    size_t u = elt.second.first;
    size_t v = elt.second.second;

    UnionFind *reprU = ufPtrVector[u]->find();
    UnionFind *reprV = ufPtrVector[v]->find();
    if (reprU == reprV) // Already in the same component
    {
      continue;
    }
    nbEdgesMerged++;
#if VERB > 0
    cout << endl << endl;
    std::cout << "considering edge " << u << "-" << v << " : " << edgeCost << endl;
#endif

#if VERB > 1
    for (size_t iPt = 0; iPt < n; iPt++)
    {
      if (abs(outputCoordsPtr[2*iPt])+abs(outputCoordsPtr[2*iPt+1]) < 1e-9)
        continue;
      std::cout << "Coords of " << iPt << ": ";
      for (size_t k = 0; k < 2; k++)
      {
        std::cout << outputCoordsPtr[iPt*+k];
        if ((int)k < (int)dim-1)
          std::cout << ",";
      }
      std::cout << "\n";
    }

    std::cout << "....................................\n\n";
#endif



    edgesMSTBefore.push_back(edgeCost);
    std::set<size_t> &compU = ufToSets[reprU], &compV = ufToSets[reprV];
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


#if VERB > 2
    std::cout << "Component Small = ";
    for (auto x : compSmall)
      std::cout << x << " ";

    std::cout << "\nComponent Big = ";
    for (auto x : compBig)
      std::cout << x << " ";
    cout << "\n FIN DES COMPOSANTES\n";
#endif

    //TODO mettre en fonction ?
    //
    // Retour de la fonction : vertexHull, coordonnéesPointsComp, chosenVert,
    // Input : distance matrix, composante (set), vertexHull, coordonnéesPointsComp, chosenVert?
    for (int idSet = 0; idSet < 2; idSet++)
    {
      //cout << (idSet == 0 ? "Comp small detailed:" : "Comp big detailed:") << "\n";
      std::vector<double> &curPointsSet = pointsSets[idSet];// stores coordinates
      std::vector<size_t> &curHullVerts = idsInHullSets[idSet];
      size_t &idChosenVert = idChosenVerts[idSet];
      size_t idCur = idSet == 0 ? idSmall:idBig;
      std::set<size_t> &curComp = idCur == u ? compU:compV;
      size_t nCur = curComp.size();
      curPointsSet.resize(nCur*2);
      std::vector<size_t> curCompVect(nCur);

      // We retrive the current coordinates of our component.
      size_t cptCur = 0;
      for (int vertId : curComp)
      {
        for (int k = 0; k < 2; k++)
        {
          curPointsSet[2*cptCur+k] = outputCoordsPtr[2*vertId+k];
        }
        curCompVect[cptCur] = vertId;
        cptCur++;
      }
      getConvexHull(curPointsSet, 2, curHullVerts);
      // The ids in curHullVerts are the index of the vertices in the component list, not
      // the real ids of the vertices. The loop just below solves this.
      for (size_t &vert : curHullVerts)
        vert = curCompVect[vert];


#if VERB > 2
      cout << endl;
      std::cout << "\n\nConvex hull : ";
      for (size_t &x : curHullVerts)
      {
        std::cout << x << " ";
      }
      std::cout << std::endl;
#endif


      idChosenVert = curHullVerts[0];
      // We want to select, among all vertices in the convex hull, the one which is
      // closest to the vertex of the edge we work on.
      for (size_t vert : curHullVerts)
      {
        double dist2 = distMatrix[vert][idEdgeVert[idSet]];
        if (dist2 < distMatrix[idChosenVert][idEdgeVert[idSet]])
          idChosenVert = vert;
      }
    }

    std::vector<size_t> &idsInHullSmall = idsInHullSets[0], idsInHullBig = idsInHullSets[1];
    std::vector<double> &pointsBig = pointsSets[1], &pointsSmall = pointsSets[0];
    size_t sizeBigHull = idsInHullSets[1].size(), sizeSmallHull = idsInHullSets[0].size();
    std::vector<double> coordsBigHull(sizeBigHull*2), coordsSmallHull(sizeSmallHull*2);
    for (size_t iHull = 0; iHull < sizeBigHull; iHull++)
    {
      size_t vert = idsInHullBig[iHull];
      coordsBigHull[iHull*2] = outputCoordsPtr[vert*2];
      coordsBigHull[iHull*2+1] = outputCoordsPtr[vert*2+1];
#if VERB > 3
      printCoords(("Vertbig " + to_string(vert)).c_str(), &coordsBigHull[2*iHull]);
#endif
    }
    for (size_t iHull = 0; iHull < sizeSmallHull; iHull++)
    {
      size_t vert = idsInHullSmall[iHull];
      coordsSmallHull[iHull*2] = outputCoordsPtr[vert*2];
      coordsSmallHull[iHull*2+1] = outputCoordsPtr[vert*2+1];
#if VERB > 3
      printCoords(("Vertsmall " + to_string(vert)).c_str(), &coordsSmallHull[2*iHull]);
#endif
    }

    size_t idChosenBig = idChosenVerts[1], idChosenSmall = idChosenVerts[0];
    double coordPrevBig[2], coordPostBig[2];
    double coordPrevSmall[2], coordPostSmall[2];
    getPrevPostEdges(idsInHullSmall, idChosenSmall, outputCoordsPtr, coordPrevSmall, coordPostSmall);
    getPrevPostEdges(idsInHullBig, idChosenBig, outputCoordsPtr, coordPrevBig, coordPostBig);
#if VERB > 3
    printCoords("prevSmall : ", coordPrevSmall);
    printCoords("postSmall : ", coordPostSmall);
#endif
    double coordPtSmall[2] = {outputCoordsPtr[2*idChosenSmall], outputCoordsPtr[2*idChosenSmall+1]};
    double coordPtBig[2] = {outputCoordsPtr[2*idChosenBig], outputCoordsPtr[2*idChosenBig+1]};

    double angleSmall = computeAngle(coordPtSmall, coordPrevSmall, coordPostSmall);
    double angleBig = computeAngle(coordPtBig, coordPrevBig, coordPostBig);
#if VERB > 1
    cout << "UUU The angles are " << deg(angleSmall) << " " << deg(angleBig) << " UUU" << endl;
#endif
    if (angleSmall > M_PI)
    {
      swap(coordPrevSmall, coordPostSmall);
      angleSmall = 2*M_PI-angleSmall;
    }
    if (angleBig > M_PI)
    {
      swap(coordPrevBig, coordPostBig);
      angleBig = 2*M_PI-angleBig;
    }
    double coordsCentreBig[2] = {coordPrevBig[0], coordPrevBig[1]};
    double coordsCentreSmall[2] = {coordPrevSmall[0], coordPrevSmall[1]};
    rotate(coordsCentreSmall, coordPtSmall, angleSmall/2);
    rotate(coordsCentreBig, coordPtBig, angleBig/2);

#if VERB > 1
    cout << "The angles REALLY are " << deg(angleSmall) << " " << deg(angleBig) << endl;
#endif

    double unitCentreBigVect[2], unitCentreSmallVect[2];
    computeUnitVector(&outputCoordsPtr[idChosenBig*2], coordsCentreBig, unitCentreBigVect);
    computeUnitVector(&outputCoordsPtr[idChosenSmall*2], coordsCentreSmall, unitCentreSmallVect);
#if VERB > 3
    printCoords("Le point choisit big : ", &outputCoordsPtr[idChosenBig*2]);
    printCoords("Centre for big comp is : ", coordsCentreBig);
    printCoords("Centre for small comp is : ", coordsCentreSmall);
#endif

    //TODO create vector from two points?
    double goalCoordChosenSmall[2] = {outputCoordsPtr[idChosenBig*2]-edgeCost*unitCentreBigVect[0], outputCoordsPtr[idChosenBig*2+1]-edgeCost*unitCentreBigVect[1]};
    if (sizeBigHull == 1)
    {
      goalCoordChosenSmall[0] = outputCoordsPtr[idChosenBig*2]+edgeCost;
      goalCoordChosenSmall[1] = outputCoordsPtr[idChosenBig*2+1];
      //continue;
    }
    //double smallCompMoveVect[2] = {coordsCentreSmall[0] + smallCompMoveVect[0], coordsCentreSmall[1] +smallCompMoveVect[1]};
    double smallCompMoveVect[2] = {goalCoordChosenSmall[0] - outputCoordsPtr[idChosenSmall*2], goalCoordChosenSmall[1] - outputCoordsPtr[idChosenSmall*2+1]};

    double distBaryPointSmall = compute_dist(&outputCoordsPtr[idChosenSmall*2], coordsCentreSmall);
#if VERB > 3
    printCoords("Le point choisi small : ", &outputCoordsPtr[idChosenSmall*2]);
    printCoords("Barycentre for small comp is : ", coordsCentreSmall);
#endif
    double preFinalPosBarySmall[2] = {coordsCentreSmall[0]+smallCompMoveVect[0], coordsCentreSmall[1]+smallCompMoveVect[1]};
    double prefinalPosBarySmall2[2] = {goalCoordChosenSmall[0]+unitCentreSmallVect[0]*distBaryPointSmall, goalCoordChosenSmall[1]+unitCentreSmallVect[1]*distBaryPointSmall};
    //TODO check identiques
    double finalPosBarySmall[2] = {goalCoordChosenSmall[0]-unitCentreBigVect[0]*distBaryPointSmall, goalCoordChosenSmall[1]-unitCentreBigVect[1]*distBaryPointSmall};

    //TODO : actuellement pour aligner les barycentres, normalement on voudrait aligner les bissectrices je pense
    //
    for (size_t curIdSmall : compSmall)
    {
#if VERB > 3
      std::cout << "pointSmall : " << outputCoordsPtr[curIdSmall*2] << "," << outputCoordsPtr[curIdSmall*2+1] << endl;
#endif
      outputCoordsPtr[curIdSmall*2] += smallCompMoveVect[0];
      outputCoordsPtr[curIdSmall*2+1] += smallCompMoveVect[1];
#if VERB > 3
      std::cout << "\tPre-new-coordinates for " << curIdSmall << " are: " << outputCoordsPtr[curIdSmall*2]<<","<<outputCoordsPtr[curIdSmall*2+1] << "\n";
#endif

      double rotationAngle = computeAngle(goalCoordChosenSmall, preFinalPosBarySmall, finalPosBarySmall);
#if VERB > 1
      std::cout << " We chose for big and small: " << idChosenBig << " and " << idChosenSmall << std::endl;
#endif
      if (nSmall > 1 && std::isfinite(rotationAngle))
      {
#if VERB > 2
        std::cout << "angle is " << rotationAngle << std::endl;
        cout << "cos,sin: " << cos(rotationAngle) << "," << sin(rotationAngle) << endl;
#endif
#if VERB > 4
        printCoords("coordCenter = ", goalCoordChosenSmall);
        printCoords("coordFrom = ", preFinalPosBarySmall);
        printCoords("coordFrom2 = ", prefinalPosBarySmall2);
        printCoords("coordTo = ", finalPosBarySmall);
#endif

        rotate(&outputCoordsPtr[curIdSmall*2], goalCoordChosenSmall, rotationAngle); //Ou - rotatinAngle ?
#if VERB > 4
        std::cout << "\t\tPost-new-coordinates for " << curIdSmall << " are: " << outputCoordsPtr[curIdSmall*2]<<","<<outputCoordsPtr[curIdSmall*2+1] << "\n";
#endif
      }
    }
    if (u == 311 && v == 362)
    {
      //cout << "returning here\n";
      //return 1;
      //break;
    }
    if (nBig > 1)
    {
      rotateMergingCompsBest(idsInHullSmall, idsInHullBig, compSmall, compBig, idChosenSmall, idChosenBig, distMatrix, outputCoordsPtr, this->AngleSamplingFreq, this->threadNumber_);
#if VERB > 4
      std::cout << "\t\tPost-new-coordinates for " << idChosenBig << " are: " << outputCoordsPtr[idChosenBig*2]<<","<<outputCoordsPtr[idChosenBig*2+1] << "\n";
#endif
    }

    double finalDist = compute_dist(&outputCoordsPtr[2*idChosenSmall], &outputCoordsPtr[2*idChosenBig]);

    if (abs(finalDist-edgeCost) > 1e-6)
      cout << "PROBLEM de distances : " << edgeCost << " against " << finalDist << endl;
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


#if VERB > 3
    cout << " Resulting component = ";
    for (auto x : unionSet)
      cout << x << " ";
    cout << "=======================\n";
#endif

    //std::cerr << "\tDELETING " << otherRepr << "( and was " <<unionRepr << ")"<< std::endl;
    ufToSets.erase(otherRepr);
    //ufToSets[otherRepr] = unionSet;
  }
  std::cout << " visited " << nbEdgesMerged << " out of " << n-1 << endl;



  std::vector<std::pair<double, std::pair<size_t, size_t>>> edgeHeapVectAfter;
  for (int u1 = 0; u1 < n; u1++)
  {
    for (int u2 = u1+1; u2 < n; u2++)
    {
      edgeHeapVectAfter.push_back({compute_dist(&outputCoordsPtr[2*u1], &outputCoordsPtr[2*u2]), {u1, u2}});
    }
  }
  sort(edgeHeapVectAfter.begin(), edgeHeapVectAfter.end());


  std::map<UnionFind*, std::set<size_t>> ufToSetsAfter;
  std::vector<UnionFind> ufVectorAfter(n);
  std::vector<UnionFind*> ufPtrVectorAfter(n);
  for (int i = 0; i < n; i++)
  {
    ufPtrVectorAfter[i] = &ufVectorAfter[i];
    ufToSetsAfter[ufPtrVectorAfter[i]].insert(i);
  }

  for (const auto &elt : edgeHeapVectAfter)
  {
    double edgeCost = elt.first;
    size_t u = elt.second.first;
    size_t v = elt.second.second;


    UnionFind *reprU = ufPtrVectorAfter[u]->find();
    UnionFind *reprV = ufPtrVectorAfter[v]->find();
    if (reprU == reprV) // Already in the same component
    {
      continue;
    }
    UnionFind* unionRepr = UnionFind::makeUnion(reprU, reprV);
    UnionFind* otherRepr = (unionRepr == reprU) ? reprV : reprU;
    if (unionRepr != reprU && unionRepr != reprV)
      std::cerr << "NOOOOOOOOOOOOOOO\n";
    std::set<size_t> &unionSet = ufToSetsAfter[unionRepr];
    std::set<size_t> &otherSet = ufToSetsAfter[otherRepr];

    unionSet.insert(otherSet.begin(), otherSet.end());
    //TODO faire des trucs;

    ufPtrVectorAfter[u] = unionRepr;
    ufPtrVectorAfter[v] = unionRepr;

    //std::cerr << "\tDELETING " << otherRepr << "( and was " <<unionRepr << ")"<< std::endl;
    ufToSetsAfter.erase(otherRepr);
    edgesMSTAfter.push_back(edgeCost);
  }

  for (int i = 0; i < edgesMSTBefore.size(); i++)
    if (abs(edgesMSTBefore[i]-edgesMSTAfter[i]) >= 1e-5)
      cout << " ERREUR SUR LARRETE " << i << " ====> " << edgesMSTBefore[i] << " VVSS " << edgesMSTAfter[i] << endl;


  for (size_t iPt = 0; iPt < n; iPt++)
  {
    if (abs(outputCoordsPtr[2*iPt])+abs(outputCoordsPtr[2*iPt+1]) < 1e-9)
      continue;
    std::cout << "Coords of " << iPt << ": ";
    for (size_t k = 0; k < dim; k++)
    {
      std::cout << outputCoordsPtr[iPt*dim+k];
      if ((int)k < (int)dim-1)
        std::cout << ",";
    }
    std::cout << "\n";
  }

  if (std::is_same<T, float>())
  {
    for (int i = 0; i < 2*n; i++)
      outputCoords[i] = outputCoordsPtr[i];
  }
  return 0;
}

} // namespace ttk
template<typename T>
void printCoords(const char prefix[], const T *coords)
{
#if VERB >= 1
  std::cout << prefix <<  coords[0] << "," << coords[1] << "\n";
#endif
}




