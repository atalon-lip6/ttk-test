/// \ingroup base
/// \class ttk::TopoMap
/// \author Alexandre Talon <alexandre.talon@lip6.fr>
/// \date March 2023
///
///
///TODO
/// This module defines the %TopoMap class that computes a
/// score indicating how good the low dimension distance matrix represents the
/// high dimension one. The score is computed according to the SIM formula.
///
///
/// \b Related \b publication: \n
/// "Topomap: A 0-dimensional homology preserving projection of high-dimensional data"\n
/// Harish Doraiswamy, Julien Tierny, Paulo J. S. Silva, Luis Gustavo Nonato, and Claudio Silva\n
/// Proc. of IEEE VIS 2020.\n
/// IEEE Transactions on Visualization and Computer Graphics 27(2): 561-571, 2020.

/// \sa TopoMap

#pragma once

#define VERB 0
#define CHECK 1

// ttk common includes
#include <Debug.h>
#include <Geometry.h>
#include <UnionFind.h>

// STL includes
#include <map>
#include <set>
#include <vector>

// Boost includes
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>



namespace bg = boost::geometry;
BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)

typedef boost::tuple<double, double> Point;
typedef boost::geometry::model::polygon<Point> Polygon;
typedef boost::geometry::model::multi_point<Point> Mpoints;



//TODO remove
static inline double deg(double angle)
{
  return (angle*180)/M_PI;
}
template<typename T>
void printCoords(const char prefix[], const T *coords)
{
  if (prefix == nullptr && coords == nullptr)
    std::cerr << "No compil warning\n";
#if VERB >= 1
  std::cout << prefix <<  coords[0] << "," << coords[1] << "\n";
#endif
}


// cos and sin computations, so we are cautious on testing equalities between double.
static double Epsilon{ttk::Geometry::pow(10.0, -FLT_DIG+2)};
static double EpsilonDBL{ttk::Geometry::pow(10.0, -DBL_DIG+2)};


template<typename T>
static inline T compute_dist2(const T* ptA, const T* ptB)
{
  T dx = ptB[0]-ptA[0], dy = ptB[1]-ptA[1];

  return dx*dx + dy*dy;
}

template<typename T>
static inline T compute_dist(const T* ptA, const T* ptB)
{
  T dx = ptB[0]-ptA[0], dy = ptB[1]-ptA[1];

  return sqrt(dx*dx + dy*dy);
}

template<typename T>
inline bool are_colinear(const T* pptA, const T* pptB, const T* pptC)
{
  double ptA[2] = {pptA[0], pptA[1]}, ptB[2] = {pptB[0], pptB[1]}, ptC[2] = {pptC[0], pptC[1]};
  return fabs(ptA[0]*(ptB[1]-ptC[1])+ptB[0]*(ptC[1]-ptA[1])+ptC[0]*(ptA[1]-ptB[1])) <= EpsilonDBL;
}

// Normalizes a given vector.
template<typename T>
static void computeUnitVector(const T* coordOrig, const T* coordDest, T* coordVect);

// Computes the oriented angle BAC
template<typename T>
static double computeAngle(const T* ptA, const T* ptB, const T* ptC);


// Rotates the first argument by angle around center.
template<typename T>
static void rotate(T* ptToRotate, const T* center, double angle);

// Rotate the set of coords by the given angle, considering the given center as center of the rotation.
template<typename T>
static void rotatePolygon(std::vector<T> &coords, T* centerCoords, const double angle);


template<typename T>
void computeConvexHull(T* allCoords, const std::vector<size_t> &compPtsIds, std::vector<size_t> &idsInHull);



namespace ttk {

  /**
   * The TopoMap class TODO
   */
  class TopoMap : virtual public Debug {

  public:
    TopoMap();
    TopoMap(size_t angularSampleNb,bool checkMST) : AngularSampleNb(angularSampleNb), CheckMST(checkMST) {}
    ~TopoMap();

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
    int execute(T* outputCoords, const std::vector<T> &inputMatrix, bool isDistMat, size_t n) const;
  protected:
    size_t AngularSampleNb{20};
    bool CheckMST{false};
  
  private:
    // Fills its last two arguments with the position of the end of the previous (resp. next) edge so that idCenter is the point at the angle at the intersection of the previous and next edges.
    template<typename T>
    void getPrevNextEdges(const std::vector<size_t> &idsPtsPolygon, size_t idCenter, const T* allCoords, T* coordPrev, T* coordPost) const;

    // Tries to find the best angle of rotation for the two components. Updates the coordiates of their vertices accordingly.
template<typename T>
T rotateMergingCompsBest(const std::vector<size_t> &hull1, const std::vector<size_t> &hull2, const std::vector<size_t> &comp1, const std::vector<size_t> &comp2, size_t iPt1, size_t iPt2, const std::vector<T> &distMatrix, T* allCoords, size_t n, size_t angularSampleNb, size_t nThread) const;

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
int ttk::TopoMap::execute(T* outputCoords, const std::vector<T> &inputMatrix, bool isDistMat, size_t n) const
{
  std::vector<double> edgesMSTBefore, edgesMSTAfter;

  // We first set the minimum value we consider as null. We multiply by 100 the smallest significative value because we use sqrt/cos/sin so we take precautionary measures.
  if (std::is_same<T, double>::value)
    Epsilon = ttk::Geometry::pow(10.0, -DBL_DIG+2);  
  else if (std::is_same<T, float>::value)
    Epsilon = ttk::Geometry::pow(10.0, -FLT_DIG+2);
  else
  {
    printErr("Error, template type must be either double or float.\n");
    return 1;
  }

  std::vector<T> computedDistMatrix(n*n);
  if (!isDistMat)
  {
    size_t dim = inputMatrix.size()/n;
    std::cout << "dim = " << dim << std::endl;
    if (n*dim != inputMatrix.size())
    {
      this->printErr("Error, the coordinates input matrix has invalid size.");
      return 1;
    }
    for (size_t i1 = 0; i1 < n; i1++)
    {
      for (size_t i2 = i1; i2 < n; i2++)
      {
        T dist2 = 0;
        for (size_t k = 0; k < dim; k++)
          dist2 += (inputMatrix[dim*i1+k]-inputMatrix[dim*i2+k])*(inputMatrix[dim*i1+k]-inputMatrix[dim*i2+k]);
        T dist = sqrt(dist2);
        computedDistMatrix[i1*n+i2] = dist;
        computedDistMatrix[i2*n+i1] = dist;
      }
    }
  }
  else
  {
    if (inputMatrix.size() != n*n)
      this->printErr("Invalid size for distance matrix.");
  }

  const std::vector<T> &distMatrix = isDistMat ? inputMatrix : computedDistMatrix;


#if VERB > 0
  std::cout << "Il y a " << n << " points! :-)\n";
#endif

  // 1. Sorting the edges
  std::vector<std::pair<T, std::pair<size_t, size_t>>> edgeHeapVect;
  edgeHeapVect.reserve(n*(n-1)/2);
  for (size_t u1 = 0; u1 < n; u1++)
  {
    for (size_t u2 = u1+1; u2 < n; u2++)
    {
      edgeHeapVect.push_back({distMatrix[u1*n+u2], {u1, u2}});
    }
  }
  sort(edgeHeapVect.begin(), edgeHeapVect.end());


  // ufVectors is the union find vector
  std::map<UnionFind*, std::vector<size_t>> ufToComp;
  std::vector<UnionFind> ufVector(n);
  for (size_t i = 0; i < n; i++)
  {
    ufToComp[&ufVector[i]].push_back(i);
  }

  for (size_t i = 0; i < 2*n; i++)
    outputCoords[i] = 0;


  T finalDistortion = 0;
  // 2. Processing the edges
  size_t nbEdgesMerged = 0;
  for (const auto &elt : edgeHeapVect)
  {
    if (nbEdgesMerged == n-1) // Our spanning tree is finished
      break;
    T edgeCost = elt.first;
    size_t u = elt.second.first;
    size_t v = elt.second.second;

    // 2.a Getting the components
    UnionFind *reprU = ufVector[u].find();
    UnionFind *reprV = ufVector[v].find();
    if (reprU == reprV) // Already in the same component
    {
      continue;
    }
    nbEdgesMerged++;
#if VERB > 0
    std::cout << std::endl << std::endl;
    std::cout << "considering edge " << u << "-" << v << " : " << edgeCost << std::endl;
#endif

#if VERB > 1
    for (size_t iPt = 0; iPt < n; iPt++)
    {
      if (fabs(outputCoords[2*iPt])+fabs(outputCoords[2*iPt+1]) < Epsilon)
        continue;
      std::cout << "Coords of " << iPt << ": ";
      std::cout << outputCoords[iPt*2] << "," << outputCoords[iPt*2+1] << "\n";
    }
    std::cout << "....................................\n\n";
#endif

    edgesMSTBefore.push_back(edgeCost);
    //TODO function
    std::vector<size_t> &compU = ufToComp[reprU], &compV = ufToComp[reprV];
    size_t idSmall = compU.size() < compV.size() ? u:v;
    size_t idBig = idSmall == u ? v:u;
    std::vector<size_t> &compSmall = idSmall == u ? compU:compV;
    std::vector<size_t> &compBig = idSmall == u ? compV:compU;
    size_t nBig = compBig.size();
    size_t nSmall = compSmall.size();
    
    std::vector<size_t> idsInHullSmall, idsInHullBig;



#if VERB > 2
    std::cout << "Component Small = ";
    for (auto x : compSmall)
      std::cout << x << " ";

    std::cout << "\nComponent Big = ";
    for (auto x : compBig)
      std::cout << x << " ";
    std::cout << "\n FIN DES COMPOSANTES\n";
#endif

    //TODO mettre en fonction ?
    //
    // 2.b Computing the convex hull.

    // We retrieve the current coordinates of the big component..
   computeConvexHull(outputCoords, compBig, idsInHullBig);
   computeConvexHull(outputCoords, compSmall, idsInHullSmall);


#if VERB > 2

      std::cout << std::endl;
      std::cout << "\n\nConvex hull small(" << idsInHullSmall.size() <<"): ";
      for (size_t &x : idsInHullSmall)
      {
        std::cout << x << " ";
      }
      std::cout << std::endl;


      std::cout << std::endl;
      std::cout << "\n\nConvex hull big(" <<idsInHullBig.size() << "): ";
      for (size_t &x : idsInHullBig)
      {
        std::cout << x << " ";
      }
      std::cout << std::endl;
#endif

      size_t idChosenBig = idsInHullBig[0], idChosenSmall = idsInHullSmall[0];
 
      // 2.c We want to select, among all vertices in the convex hull, the one which is
      // closest to the vertex of the edge we work on.
      for (size_t vert : idsInHullBig)
      {
        T dist = distMatrix[vert*n+idBig];
        if (dist < distMatrix[idChosenBig*n+idBig])
          idChosenBig = vert;
      }

      for (size_t vert : idsInHullSmall)
      {
        T dist = distMatrix[vert*n+idSmall];
        if (dist < distMatrix[idChosenSmall*n+idSmall])
          idChosenSmall = vert;
      }

    size_t sizeBigHull = idsInHullBig.size(), sizeSmallHull = idsInHullSmall.size();
    std::vector<T> coordsBigHull(sizeBigHull*2), coordsSmallHull(sizeSmallHull*2);
    for (size_t iHull = 0; iHull < sizeBigHull; iHull++)
    {
      size_t vert = idsInHullBig[iHull];
      coordsBigHull[iHull*2] = outputCoords[vert*2];
      coordsBigHull[iHull*2+1] = outputCoords[vert*2+1];
#if VERB > 3
      printCoords(("Vertbig " + std::to_string(vert)).c_str(), &coordsBigHull[2*iHull]);
#endif
    }
    for (size_t iHull = 0; iHull < sizeSmallHull; iHull++)
    {
      size_t vert = idsInHullSmall[iHull];
      coordsSmallHull[iHull*2] = outputCoords[vert*2];
      coordsSmallHull[iHull*2+1] = outputCoords[vert*2+1]; //TODO put in computeconvexhull
#if VERB > 3
      printCoords(("Vertsmall " + std::to_string(vert)).c_str(), &coordsSmallHull[2*iHull]);
#endif
    }

    
    // Identifying the angles we are working on from the convex hulls.
    T coordPrevBig[2], coordPostBig[2];
    T coordPrevSmall[2], coordPostSmall[2];
    getPrevNextEdges(idsInHullSmall, idChosenSmall, outputCoords, coordPrevSmall, coordPostSmall);
    getPrevNextEdges(idsInHullBig, idChosenBig, outputCoords, coordPrevBig, coordPostBig);
#if VERB > 3
    printCoords("prevSmall : ", coordPrevSmall);
    printCoords("postSmall : ", coordPostSmall);
#endif
    T coordPtSmall[2] = {outputCoords[2*idChosenSmall], outputCoords[2*idChosenSmall+1]};
    T coordPtBig[2] = {outputCoords[2*idChosenBig], outputCoords[2*idChosenBig+1]};

    // Computing the angles.
    double angleSmall = computeAngle(coordPtSmall, coordPrevSmall, coordPostSmall);
    double angleBig = computeAngle(coordPtBig, coordPrevBig, coordPostBig);
#if VERB > 1
    std::cout << "UUU The angles are " << deg(angleSmall) << " " << deg(angleBig) << " UUU" << std::endl;
#endif
    //TODO put that in getprevpost?
    if (angleSmall > M_PI)
    {
      std::swap(coordPrevSmall, coordPostSmall);
      angleSmall = 2*M_PI-angleSmall;
    }
    if (angleBig > M_PI)
    {
      std::swap(coordPrevBig, coordPostBig);
      angleBig = 2*M_PI-angleBig;
    }
    T coordscenterBig[2] = {coordPrevBig[0], coordPrevBig[1]};
    T coordscenterSmall[2] = {coordPrevSmall[0], coordPrevSmall[1]};
    // Computing the coordinates of the bisectors.
    rotate(coordscenterSmall, coordPtSmall, angleSmall/2);
    rotate(coordscenterBig, coordPtBig, angleBig/2);

#if VERB > 1
    std::cout << "The angles REALLY are " << deg(angleSmall) << " " << deg(angleBig) << std::endl;
#endif

    T unitcenterBigVect[2], unitcenterSmallVect[2];
    computeUnitVector(&outputCoords[idChosenBig*2], coordscenterBig, unitcenterBigVect);
    computeUnitVector(&outputCoords[idChosenSmall*2], coordscenterSmall, unitcenterSmallVect);
#if VERB > 3
    printCoords("Le point choisit big : ", &outputCoords[idChosenBig*2]);
    printCoords("center for big comp is : ", coordscenterBig);
    printCoords("center for small comp is : ", coordscenterSmall);
#endif

    // Computing the new coordinates for the chosen point on the small components, and then translating the whole component.
    T goalCoordChosenSmall[2] = {outputCoords[idChosenBig*2]-edgeCost*unitcenterBigVect[0], outputCoords[idChosenBig*2+1]-edgeCost*unitcenterBigVect[1]};
    if (sizeBigHull == 1)
    {
      goalCoordChosenSmall[0] = outputCoords[idChosenBig*2]+edgeCost;
      goalCoordChosenSmall[1] = outputCoords[idChosenBig*2+1];
      //continue;
    }
    //double smallCompMoveVect[2] = {coordscenterSmall[0] + smallCompMoveVect[0], coordscenterSmall[1] +smallCompMoveVect[1]};
    T smallCompMoveVect[2] = {goalCoordChosenSmall[0] - outputCoords[idChosenSmall*2], goalCoordChosenSmall[1] - outputCoords[idChosenSmall*2+1]};

    T distBaryPointSmall = compute_dist(&outputCoords[idChosenSmall*2], coordscenterSmall);
#if VERB > 3
    printCoords("Le point choisi small : ", &outputCoords[idChosenSmall*2]);
    printCoords("Barycenter for small comp is : ", coordscenterSmall);
#endif
    T preFinalPosBarySmall[2] = {coordscenterSmall[0]+smallCompMoveVect[0], coordscenterSmall[1]+smallCompMoveVect[1]};
    //TODO check identiques
    T finalPosBarySmall[2] = {goalCoordChosenSmall[0]-unitcenterBigVect[0]*distBaryPointSmall, goalCoordChosenSmall[1]-unitcenterBigVect[1]*distBaryPointSmall};

    // Performing the translation + rotation such that the bisectors in compSmall and compBig are aligned, and the two components are facing each other, not crossing.
    for (size_t curIdSmall : compSmall)
    {
#if VERB > 3
      std::cout << "pointSmall : " << outputCoords[curIdSmall*2] << "," << outputCoords[curIdSmall*2+1] << std::endl;
#endif
      outputCoords[curIdSmall*2] += smallCompMoveVect[0];
      outputCoords[curIdSmall*2+1] += smallCompMoveVect[1];
#if VERB > 3
      std::cout << "\tPre-new-coordinates for " << curIdSmall << " are: " << outputCoords[curIdSmall*2]<<","<<outputCoords[curIdSmall*2+1] << "\n";
#endif

      // 2.e Performing the rotation.
      double rotationAngle = computeAngle(goalCoordChosenSmall, preFinalPosBarySmall, finalPosBarySmall);
#if VERB > 1
      std::cout << " We chose for big and small: " << idChosenBig << " and " << idChosenSmall << std::endl;
#endif
      if (nSmall > 1 && std::isfinite(rotationAngle))
      {
#if VERB > 2
        std::cout << "angle is " << rotationAngle << std::endl;
        std::cout << "cos,sin: " << cos(rotationAngle) << "," << sin(rotationAngle) << std::endl;
#endif
#if VERB > 4
        printCoords("coordCenter = ", goalCoordChosenSmall);
        printCoords("coordFrom = ", preFinalPosBarySmall);
        printCoords("coordFrom2 = ", prefinalPosBarySmall2);
        printCoords("coordTo = ", finalPosBarySmall);
#endif

        rotate(&outputCoords[curIdSmall*2], goalCoordChosenSmall, rotationAngle); //Ou - rotatinAngle ?
#if VERB > 4
        std::cout << "\t\tPost-new-coordinates for " << curIdSmall << " are: " << outputCoords[curIdSmall*2]<<","<<outputCoords[curIdSmall*2+1] << "\n";
#endif
      }
    }

    // 2.f Trying several rotations of the two components and keeping the one which makes the new distance matrix closest to the one provided in the input.
    if (nBig > 1)
    {
      finalDistortion = rotateMergingCompsBest(idsInHullSmall, idsInHullBig, compSmall, compBig, idChosenSmall, idChosenBig, distMatrix, outputCoords, n, this->AngularSampleNb, this->threadNumber_);
#if VERB > 4
      std::cout << "\t\tPost-new-coordinates for " << idChosenBig << " are: " << outputCoords[idChosenBig*2]<<","<<outputCoords[idChosenBig*2+1] << "\n";
#endif
    }

    T finalDist = compute_dist(&outputCoords[2*idChosenSmall], &outputCoords[2*idChosenBig]);

    if (fabs(finalDist-edgeCost) > Epsilon)
      std::cout << "PROBLEM de distances : " << edgeCost << " against " << finalDist << std::endl;
    UnionFind* unionRepr = UnionFind::makeUnion(reprU, reprV);
    UnionFind* otherRepr = (unionRepr == reprU) ? reprV : reprU;
    if (unionRepr != reprU && unionRepr != reprV)
    {
      printErr("Error with the union find module results. Aborting.");
      return 1;
    }
    std::vector<size_t> &unionSet = ufToComp[unionRepr];
    std::vector<size_t> &otherSet = ufToComp[otherRepr];

    unionSet.insert(unionSet.end(), otherSet.begin(), otherSet.end());

#if VERB > 3
    std::cout << " Resulting component = ";
    for (auto x : unionSet)
      std::cout << x << " ";
    std::cout << "=======================\n";
#endif

    ufToComp.erase(otherRepr);
  }
  this->printMsg("Non normalized distance matrix distortion: "+std::to_string(2*finalDistortion));
  std::cout << " visited " << nbEdgesMerged << " out of " << n-1 << std::endl;



  // We check that the lengths of the edges selected to build a minimum spanning tree
  // are preserved by our algorithm, as should be.
  if (CheckMST)
  {
    std::vector<std::pair<T, std::pair<size_t, size_t>>> edgeHeapVectAfter;
    for (size_t u1 = 0; u1 < n; u1++)
    {
      for (size_t u2 = u1+1; u2 < n; u2++)
      {
        edgeHeapVectAfter.push_back({compute_dist(&outputCoords[2*u1], &outputCoords[2*u2]), {u1, u2}});
      }
    }
    sort(edgeHeapVectAfter.begin(), edgeHeapVectAfter.end());


    std::map<UnionFind*, std::vector<size_t>> ufToCompAfter;
    std::vector<UnionFind> ufVectorAfter(n);
    std::vector<UnionFind*> ufPtrVectorAfter(n);
    for (size_t i = 0; i < n; i++)
    {
      ufPtrVectorAfter[i] = &ufVectorAfter[i];
      ufToCompAfter[ufPtrVectorAfter[i]].push_back(i);
    }


    std::cout << " checking " << std::endl;
    for (const auto &elt : edgeHeapVectAfter)
    {
      T edgeCost = elt.first;
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
      std::vector<size_t> &unionSet = ufToCompAfter[unionRepr];
      std::vector<size_t> &otherSet = ufToCompAfter[otherRepr];

      unionSet.insert(unionSet.end(), otherSet.begin(), otherSet.end());
      //TODO faire des trucs;

      ufPtrVectorAfter[u] = unionRepr;
      ufPtrVectorAfter[v] = unionRepr;

      ufToCompAfter.erase(otherRepr);
      edgesMSTAfter.push_back(edgeCost);
    }

    for (size_t i = 0; i < edgesMSTBefore.size(); i++)
      if (fabs(edgesMSTBefore[i]-edgesMSTAfter[i]) >= Epsilon)
        std::cout << " ERREUR SUR LARRETE " << i << " ====> " << edgesMSTBefore[i] << " VVSS " << edgesMSTAfter[i] << std::endl;


  }
  for (size_t iPt = 0; iPt < n; iPt++)
  {
    if (fabs(outputCoords[2*iPt])+fabs(outputCoords[2*iPt+1]) < Epsilon)
      continue;
    std::cout << "Coords of " << iPt << ": ";
    std::cout << outputCoords[iPt*2] << "," << outputCoords[iPt*2+1] << "\n";
  }

  return 0;
}


template<typename T>
void TopoMap::getPrevNextEdges(const std::vector<size_t> &idsPtsPolygon, size_t idCenter, const T* allCoords, T* coordPrev, T* coordPost) const
{
  size_t n = idsPtsPolygon.size();
  size_t iPtPrev = 0, iPtPost = 0;
  bool found = false;

#if VERB > 4
  for (size_t i = 0; i < n; i++)
    printCoords("", &allCoords[idsPtsPolygon[i]*2]);
#endif
  for (size_t i = 0; i < n; i++)
  {
    if (idsPtsPolygon[i] == idCenter)
    {
      found = true;
      iPtPost = idsPtsPolygon[(i+1)%n];
      iPtPrev = idsPtsPolygon[(i+n-1)%n];
      break;
    }
  }
  if (!found)
    printErr("Error, we could not find the edges incident to the point we chose for the rotation of the component.");

  double angle = computeAngle(&allCoords[2*idCenter], &allCoords[2*iPtPrev], &allCoords[2*iPtPost]);
  if (angle > M_PI)
  {
    std::cout << "swaaaaaap\n";
    std::swap(iPtPrev, iPtPost);
  }
  coordPrev[0] = allCoords[2*iPtPrev];
  coordPrev[1] = allCoords[2*iPtPrev+1];
  coordPost[0] = allCoords[2*iPtPost];
  coordPost[1] = allCoords[2*iPtPost+1];
}

// Tries to find the best angle of rotation for the two components. Updates the coordiates of their vertices accordingly.
template<typename T>
T TopoMap::rotateMergingCompsBest(const std::vector<size_t> &hull1, const std::vector<size_t> &hull2, const std::vector<size_t> &comp1, const std::vector<size_t> &comp2, size_t iPt1, size_t iPt2, const std::vector<T> &distMatrix, T* allCoords, size_t n, size_t angularSampleNb, size_t nThread) const
{
  TTK_FORCE_USE(nThread);
  // The distance between the two components.
  T shortestDistPossible = compute_dist(&allCoords[2*iPt1], &allCoords[2*iPt2]);
  T coordPt1[2] = {allCoords[2*iPt1], allCoords[2*iPt1+1]};
  T coordPt2[2] = {allCoords[2*iPt2], allCoords[2*iPt2+1]};
  size_t comp1Size = comp1.size(), comp2Size = comp2.size();

  T coordPrev1[2], coordPost1[2];
  T coordPrev2[2], coordPost2[2];
  getPrevNextEdges(hull1, iPt1, allCoords, coordPrev1, coordPost1);
  getPrevNextEdges(hull2, iPt2, allCoords, coordPrev2, coordPost2);

#if VERB > 3
  std::cout << "Comp 1 = ";
  for (int x : comp1)
    std::cout << x << " ";
  std::cout << std::endl;
  std::cout << "Comp 2 = ";
  for (int x : comp2)
    std::cout << x << " ";
  std::cout << std::endl;
  //std::cout << "prev2, cur2, post2 : " << iPtPrev2 << ", " << iPt2 << ", " << iPtPost2 << std::endl;
#endif

  double angle1 = computeAngle(coordPt1, coordPrev1, coordPost1);
  double angle2 = computeAngle(coordPt2, coordPrev2, coordPost2);
#if VERB > 2
  std::cout << "The angles are " << deg(angle1) << " " << deg(angle2) << std::endl;
#endif
  if (angle1 > M_PI || angle2 > M_PI)
    printErr("Angle, zut quoi!");
#if VERB > 2
  std::cout << "The angles REALLY are " << deg(angle1) << " " << deg(angle2) << std::endl;
#endif

  T coordBissect1[2] = {coordPrev1[0], coordPrev1[1]};
  T coordBissect2[2] = {coordPrev2[0], coordPrev2[1]};
  rotate(coordBissect1, coordPt1, angle1/2);
  rotate(coordBissect2, coordPt2, angle2/2);

  double semiAngle1 = computeAngle(coordPt1, coordPrev1, coordBissect1);
  double semiAngle2 = computeAngle(coordPt2, coordPrev2, coordBissect2);

  if (fabs(semiAngle1 - angle1/2) > Epsilon)
    std::cout << "çavapas angle1 et semi1 : " << semiAngle1 <<  " VS " << angle1/2 << "\n";
  if (fabs(semiAngle2 - angle2/2) > Epsilon)
    std::cout << "çavapas angle2 et semi2 : " << semiAngle2 <<  " VS " << angle2/2 << "\n";

#if VERB > 3
  printCoords("Pt1 : ", coordPt1);
  printCoords("prev1 : ", coordPrev1);
  printCoords("post1 : ", coordPost1);
  printCoords("bissect1 : ", coordBissect1);

  printCoords("Pt2 : ", coordPt2);
  printCoords("prev2 : ", coordPrev2);
  printCoords("post2 : ", coordPost2);
  printCoords("bissect2 : ", coordBissect2);
#endif
  double angleMax1 = M_PI/2-semiAngle1, angleMin1 = -angleMax1; //TODO calculer autrement et comparer
  double angleMax2 = M_PI/2-semiAngle2, angleMin2 = -angleMax2; //TODO calculer autrement et comparer

  if (angleMax1 < angleMin1)
    std::cout << "bizarre, max < min: " << deg(angleMin1) << " and " << deg(angleMax1) << "\n";

  if (angleMax2 < angleMin2)
    std::cout << "bizarre, max < min: " << deg(angleMin2) << " and " << deg(angleMax2) << "\n";
  double step1 = (angleMax1-angleMin1)/angularSampleNb, step2 = (angleMax2-angleMin2)/angularSampleNb;
  double bestAnglePair[2] = {0,0};
  T bestScore = 3.e38; // Near the maximum value for float

 
  std::vector<std::vector<T>> origDistMatrix(comp1Size);
  for (size_t i = 0; i < comp1Size; i++)
  {
    origDistMatrix[i].resize(comp2Size);
    for (size_t j = 0; j < comp2Size; j++)
    {
      origDistMatrix[i][j] = distMatrix[comp1[i]*n+comp2[j]];
    }
  }
  std::vector<T> initialCoords1(2*comp1Size), initialCoords2(2*comp2Size);
  for (size_t i = 0; i < comp1.size(); i++)
  {
    size_t iPt = comp1[i];
    initialCoords1[2*i] = allCoords[2*iPt];
    initialCoords1[2*i+1] = allCoords[2*iPt+1];
  }
  for (size_t i = 0; i < comp2.size(); i++)
  {
    size_t iPt = comp2[i];
    initialCoords2[2*i] = allCoords[2*iPt];
    initialCoords2[2*i+1] = allCoords[2*iPt+1];
  }




#if VERB > 0
  std::cout << "Angles for 1 (small) are to rotate min max " << deg(angleMin1) << " et " << deg(angleMax1) << std::endl;
  std::cout << "Angles for 2 (big) are to rotate min max " << deg(angleMin2) << " et " << deg(angleMax2) << std::endl;
#endif
  size_t nbIter1 = std::isfinite(step1) ? angularSampleNb+1 : 1;
  size_t nbIter2 = std::isfinite(step2) ? angularSampleNb+1 : 1;

  if (step1*nbIter1 < 0.001) // No need to split such a small angle
    nbIter1 = 1;
  if (step2*nbIter2 < 0.001) // No need to split such a small angle
    nbIter2 = 1;
#ifdef TTK_ENABLE_OMP
#pragma omp parallel for num_threads(nThread) shared(allCoords)
#endif
  for (size_t i1 = 0; i1 < nbIter1; i1++)
  {
    std::vector<T> coords1Test(2*comp1Size), coords2Test(2*comp2Size);
    coords1Test = initialCoords1;
    double testAngle1 = angleMin1+step1*i1;
    rotatePolygon(coords1Test, coordPt1, testAngle1);
#if VERB > 2
    std::cout << "\t\t\t\tTesting angle1 " << deg(testAngle1) << std::endl;
#endif
#if CHECK
    if (comp1Size >= 2)
    {
      T *A1, *B1, *C1;
      auto it1 = comp1.begin();
      int cpt1 = 0;
      while (*it1 != iPt1)
      {
        cpt1++;
        ++it1;
      }
      A1 = &coords1Test[2*cpt1];
      if (it1 == comp1.begin())
      {
        ++it1;
        cpt1++;
      }
      else
      {
        --it1;
        cpt1--;
      }
      B1 = &initialCoords1[2*(cpt1)];
      C1 = &coords1Test[2*cpt1];
      double totoAngle = computeAngle(A1,B1,C1);
      if (totoAngle > M_PI)
        totoAngle -= 2*M_PI;
      if (fabs(deg(testAngle1) - deg(totoAngle)) > 0.5)
        std::cout << "ERROR1 : " << deg(testAngle1) << " VS " << deg(totoAngle) << "\n";
    }
#endif

    for (size_t i2 = 0; i2 < nbIter2; i2++)
    {
      coords2Test = initialCoords2;
      double testAngle2 = angleMin2+i2*step2;
#if VERB > 2
      std::cout << "\t\t\t\tTesting angle2 " << deg(testAngle2) << std::endl;
#endif
      rotatePolygon(coords2Test, coordPt2, testAngle2);

      T curScore = 0;
      //TODO dans fonction...
      for (size_t i = 0; i < comp1Size; i++)
      {
        T coordARotate[2] = {coords1Test[2*i],coords1Test[2*i+1]};
        for (size_t j = 0; j < comp2Size; j++)
        {
          T coordBRotate[2] = {coords2Test[2*j],coords2Test[2*j+1]};
          T newDist = compute_dist(coordARotate, coordBRotate);
          curScore += (newDist-origDistMatrix[i][j])*(newDist-origDistMatrix[i][j]);
          if (newDist+Epsilon < shortestDistPossible)
          {
            std::cout << "problem " << newDist << '(' << comp1[i] << ',' << comp2[j] << ')' << " is lower than " << shortestDistPossible << " =====> (" << i1 << "," << i2 << ") (angles iterations)" << std::endl;
            printCoords(" comp1 pt : ", coordARotate);
            printCoords(" comp2 pt : ", coordBRotate);

          }
        }
      }
#if VERB > 7
      std::cout << "New matrix : \n";
      for (int ia = 0; ia < comp1Size+comp2Size; ia++)
      {
        for (int ib = 0; ib < comp1Size+comp2Size; ib++)
        {
          std::cout << newDistMatrix[ia][ib] << "\t\t\t\t";
        }
        std::cout << std::endl;
      }
#endif
#if CHECK
      if (comp2Size >= 2 && comp1Size >= 1)
      {
        T *A, *B, *C;
        auto it2 = comp2.begin();
        int cpt2 = 0;
        while (*it2 != iPt2)
        {
          cpt2++;
          ++it2;
        }
        A = &coords2Test[2*cpt2];
        if (it2 == comp2.begin())
        {
          ++it2;
          cpt2++;
        }
        else
        {
          --it2;
          cpt2--;
        }
        B = &initialCoords2[2*(cpt2)];
        C = &coords2Test[2*cpt2];
        double totoAngle = computeAngle(A,B,C);
        if (totoAngle > M_PI)
          totoAngle -= 2*M_PI;
        if (fabs(deg(testAngle2) - deg(totoAngle)) > 0.5)
          std::cout << "ERROR2 : " << deg(testAngle2) << " VS " << deg(totoAngle) << "\n";
      }
#endif

#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif
{
  //TODO maybe juste update bestScore et bestAnglePair, calculer coords à la fin ?
      if (curScore < bestScore || (curScore == bestScore && (testAngle1 < bestAnglePair[0] || (testAngle1 == bestAnglePair[0] && testAngle2 < bestAnglePair[1])))) // This pair of angle minimises the distortion and the choice is deterministic with threads
      {
#if VERB > 0
        std::cout << "Found new best score = " << curScore << " with angles " << deg(testAngle1) << ", " << deg(testAngle2) << std::endl;
#endif
        //TODO réduire et refaire à la fin
        bestScore = curScore;
        bestAnglePair[0] = std::isfinite(testAngle1) ? testAngle1 : 0;
        bestAnglePair[1] = std::isfinite(testAngle2) ? testAngle2 : 0;
      }
}
    }
  }

  //Todo rotatepolygon
  for (size_t i1 : comp1)
    rotate(&allCoords[2*i1], coordPt1, bestAnglePair[0]);

  for (size_t i2 : comp2)
    rotate(&allCoords[2*i2], coordPt2, bestAnglePair[1]);

#if VERB > 1
  std::cout << "best score is " << bestScore << std::endl;

  std::cout << "The best angles are " << deg(bestAnglePair[0]) << " for 1 and " << deg(bestAnglePair[1]) << " for 2." << std::endl;
#endif

  return bestScore;
}

} // namespace ttk

template<typename T>
static void computeUnitVector(const T* coordOrig, const T* coordDest, T* coordVect)
{
  T tmp[2] = {coordDest[0] - coordOrig[0], coordDest[1] - coordOrig[1]};
  T dist = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1]);
  if (dist < Epsilon)
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

template<typename T>
static double computeAngle(const T* ptA, const T* ptB, const T* ptC)
{
  double angle;
  double vect1[2] = {ptB[0]-ptA[0], ptB[1]-ptA[1]};
  double vect2[2] = {ptC[0]-ptA[0], ptC[1]-ptA[1]};

#if VERB > 7
  printCoords("A : ", ptA);
  printCoords("B : ", ptB);
  printCoords("C : ", ptC);
#endif
  angle = atan2(vect2[1], vect2[0]) - atan2(vect1[1], vect1[0]);
  if (angle < 0)
    angle += 2 * M_PI;
  return angle;
}

template<typename T>
static void rotate(T* ptToRotate, const T* center, double angle)
{
  const double &xCtr = center[0], &yCtr = center[1];
  T &xPt = ptToRotate[0], &yPt = ptToRotate[1];
  const double dx = xPt-xCtr, dy = yPt-yCtr;
  xPt = dx*cos(angle)-dy*sin(angle)+xCtr;
  yPt = dx*sin(angle)+dy*cos(angle)+yCtr;
}

template<typename T>
static void rotatePolygon(std::vector<T> &coords, T* centerCoords, const double angle)
{
  double xCenter = centerCoords[0], yCenter = centerCoords[1];
  size_t nbPoint = coords.size()/2;
#if VERB > 7
  std::cout << "center = " << xCenter << "," << yCenter << std::endl;
#endif
  for (size_t iPt = 0; iPt < nbPoint; iPt++)
  {
    T &x = coords[iPt*2], &y = coords[iPt*2+1];
    double xNew, yNew;
    xNew = (x-xCenter)*cos(angle)-(y-yCenter)*sin(angle)+xCenter;
    yNew = (y-yCenter)*cos(angle)+(x-xCenter)*sin(angle)+yCenter;
#if VERB > 7
    std::cout << "point = " << x << "," << y << " ===> (" << xNew << ", " << yNew << ")" <<  std::endl;
#endif
    x = xNew;
    y = yNew;
  }
}

template<typename T>
void computeConvexHull(T* allCoords, const std::vector<size_t> &compPtsIds, std::vector<size_t> &idsInHull)
{
  //TODO copier coords puis sort
  size_t nbPoint = compPtsIds.size();
  std::vector<T> compCoords(2*nbPoint);
  for (size_t i = 0; i < nbPoint; i++)
  {
    size_t vertId = compPtsIds[i];
    compCoords[2*i] = allCoords[2*vertId];
    compCoords[2*i+1] = allCoords[2*vertId+1];
  }


  if (nbPoint <= 2)
  {
    idsInHull.push_back(compPtsIds[0]);
    if (nbPoint == 2)
    {
      T dist = compute_dist(&compCoords[0], &compCoords[2]);

      if (dist > Epsilon)
        idsInHull.push_back(compPtsIds[1]);
    }
    return;
  }

  // Testing if all points are colinear

  double dirVect[2] = {0,0};
  bool areColinear = true;

  size_t idFirstDistinct = 0;
  while (idFirstDistinct < nbPoint && fabs(dirVect[0]) < Epsilon && fabs(dirVect[1]) < Epsilon)
  {
    idFirstDistinct++;
    dirVect[0] = compCoords[2*idFirstDistinct]-compCoords[0];
    dirVect[1] = compCoords[2*idFirstDistinct+1]-compCoords[1]; //TODO tester avec deux points mêmes coordonnées !
    if (fabs(dirVect[0]) < Epsilon)
      dirVect[0] =  0;
    if (fabs(dirVect[1]) < Epsilon)
      dirVect[1] =  0;
  }

  int idMins[2] = {compCoords[0] < compCoords[2] ? 0:1, compCoords[1] < compCoords[3] ? 0:1};
  int idMaxs[2] = {compCoords[0] > compCoords[2] ? 0:1, compCoords[1] > compCoords[3] ? 0:1};

  const T* pt0 = &compCoords[0];
  const T* ptDistinct = &compCoords[2*idFirstDistinct];

  for (size_t iPt = idFirstDistinct+1; iPt < nbPoint; iPt++)
  {
    T curVect[2] = {compCoords[2*iPt]-compCoords[0], compCoords[2*iPt+1]-compCoords[1]};
    if (fabs(curVect[0]) < Epsilon && fabs(curVect[1]) < Epsilon)
      continue;
    const T* ptCur = &compCoords[2*iPt];
    if (!are_colinear(pt0, ptDistinct, ptCur))
    {
      areColinear = false;
      break;
    }

    if (compCoords[2*iPt] < compCoords[2*idMins[0]])
      idMins[0] = iPt;
    if (compCoords[2*iPt+1] < compCoords[2*idMins[1]])
      idMins[1] = iPt;
    if (compCoords[2*iPt] > compCoords[2*idMaxs[0]])
      idMaxs[0] = iPt;
    if (compCoords[2*iPt+1] > compCoords[2*idMaxs[1]])
      idMaxs[1] = iPt;

  }

  if (areColinear)
  {
#if VERB > 5
    std::cout << "COLINEAR = " << idMins[0] << "," << idMins[1] << "  ;  " << idMaxs[0] << "," << idMaxs[1] << std::endl;
#endif
    if (fabs(dirVect[0]) > Epsilon)
    {
      idsInHull.push_back(compPtsIds[idMins[0]]);
      idsInHull.push_back(compPtsIds[idMaxs[0]]);
    }
    else
    {
      idsInHull.push_back(compPtsIds[idMins[1]]);
      idsInHull.push_back(compPtsIds[idMaxs[1]]);
    }
    return;
  }



  Mpoints multiPoints;
  Polygon hull;
  for (size_t i = 0; i < nbPoint; i++)
  {
    boost::geometry::append(multiPoints, Point(compCoords[2*i], compCoords[2*i+1]));
  }
  boost::geometry::convex_hull(multiPoints, hull);
  for (auto boostPt : hull.outer())
  {
    T coordsCur[2] = {static_cast<T>(boostPt.get<0>()), static_cast<T>(boostPt.get<1>())};
    for (size_t j = 0; j < compCoords.size()/2; j++)
    {
      if (fabs(compCoords[2*j]-coordsCur[0])+fabs(compCoords[2*j+1]-coordsCur[1]) < Epsilon)
      {
        idsInHull.push_back(compPtsIds[j]);
        break;
      }
    }
  }
  // Boost closes the polygon, hence the first and the last vertices are the identical.
  idsInHull.pop_back();
}

