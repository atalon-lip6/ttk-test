#include <TopologicalMapper.h>
#include <UnionFind.h>
#include <vector>
#include <set>
#include <queue>
#include <map>
#include <utility>
#include <algorithm>
#include <Geometry.h>
#include <DistanceMatrixDistorsion.h>
#include "libqhullcpp/Qhull.h"

#include <cmath>

using namespace std;

size_t bestAngleSampleFreq = 5;

inline float deg(float angle)
{
  return (angle*180)/M_PI;
}
void printCoords(const char prefix[], float const *coords)
{
  std::cout << prefix <<  coords[0] << "," << coords[1] << "\n";
}


inline float compute_dist2(const float ptA[], const float ptB[])
{
  float dx = ptB[0]-ptA[0], dy = ptB[1]-ptA[1];

  return dx*dx + dy*dy;
}



float computeAngle(float const ptA[], float const ptB[], float const ptC[], bool is3D = false)
{
  float angle;
  float vect1[2] = {ptB[0]-ptA[0], ptB[1]-ptA[1]};
  float vect2[2] = {ptC[0]-ptA[0], ptC[1]-ptA[1]};
  double dirVect[2] = {0,0};

  if (abs(vect1[0]*vect2[1] - vect2[0]*vect1[1]) < 1e-7)
  {
    cout << " Catched null ANGLE :-)\n";
    return 0;
  }
  angle = atan2(vect2[1], vect2[0]) - atan2(vect1[1], vect1[0]);
  if (angle < 0)
    angle += 2 * M_PI;
  return angle;

  /*
  float angle;
  if (is3D)
  {
    angle = ttk::Geometry::angle(ptA, ptB, ptA, ptC);
  }
  else
  {
    float ptA3D[3] = {ptA[0], ptA[1], 0}, ptB3D[3] = {ptB[0], ptB[1], 0}, ptC3D[3] = {ptC[0], ptC[1], 0};
    printCoords("A", ptA3D);
    printCoords("B", ptB3D);
    printCoords("C", ptC3D);
    if (compute_dist2(ptB3D, ptC3D) < 1e-6) // colinear => better test TODO
      return 0;
    angle = ttk::Geometry::angle(ptA3D, ptB3D, ptA3D, ptC3D);
  }

  return angle;
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

void rotatePolygon(std::vector<float> &coords, size_t dim, float* centerCoords, const float angle)
{
  float xCenter = centerCoords[0], yCenter = centerCoords[1];
  size_t nbPoint = coords.size()/dim;
  //std::cout << "center = " << xCenter << "," << yCenter << std::endl;
  for (size_t iPt = 0; iPt < nbPoint; iPt++)
  {
    float &x = coords[iPt*dim], &y = coords[iPt*dim+1];
    float xNew, yNew;
    //std::cout << "point = " << x << "," << y << std::endl;
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

void rotateMergingCompsBest(const std::vector<size_t> &hull1, const std::vector<size_t> &hull2, const std::set<size_t> &comp1, const std::set<size_t> &comp2, size_t iPt1, size_t iPt2, const std::vector<std::vector<float>> &distMatrix, float* allCoords)
{
  float coordPt1[2] = {allCoords[2*iPt1], allCoords[2*iPt1+1]};
  float coordPt2[2] = {allCoords[2*iPt2], allCoords[2*iPt2+1]};
  size_t n1 = hull1.size(), n2 = hull2.size();

  size_t iPtPrev1, iPtPrev2, iPtPost1, iPtPost2;
  for (size_t i1 = 0; i1 < hull1.size(); i1++)
  {
    if (hull1[i1] == iPt1)
    {
      iPtPost1 = hull1[(i1+1)%n1];
      iPtPrev1 = hull1[(i1+n1-1)%n1];
      break;
    }
  }
  for (size_t i2 = 0; i2 < hull2.size(); i2++)
  {
    if (hull2[i2] == iPt2)
    {
      iPtPost2 = hull2[(i2+1)%n2];
      iPtPrev2 = hull2[(i2 + n2 -1)%n2];
      break;
    }
  }
*/
  float coordPrev1[2], coordPost1[2];
  float coordPrev2[2], coordPost2[2];
  getPrevPostEdges(hull1, iPt1, allCoords, coordPrev1, coordPost1);
  getPrevPostEdges(hull2, iPt2, allCoords, coordPrev2, coordPost2);

  //cout << "prev2, cur2, post2 : " << iPtPrev2 << ", " << iPt2 << ", " << iPtPost2 << endl;

  float angle2 = computeAngle(coordPt2, coordPrev2, coordPost2, false);
  float angle1 = computeAngle(coordPt1, coordPrev1, coordPost1, false);
  cout << "The angles are " << angle1 << " " << angle2 << endl;
  if (angle1 > M_PI)
  {
    swap(coordPrev1, coordPost1);
    angle1 = 2*M_PI-angle1;
  }
  if (angle2 > M_PI)
  {
    swap(coordPrev2, coordPost2);
    angle2 = 2*M_PI-angle2;
  }
  cout << "The angles REALLY are " << angle1 << " " << angle2 << endl;

  float unitVectPrev1[2], unitVectPrev2[2], unitVectPost1[2], unitVectPost2[2];
  computeUnitVector(coordPt1, coordPrev1, unitVectPrev1);
  computeUnitVector(coordPt2, coordPrev2, unitVectPrev2);
  computeUnitVector(coordPt1, coordPrev1, unitVectPrev1);
  computeUnitVector(coordPt2, coordPrev2, unitVectPrev2);

  float coordBissect1[2] = {coordPrev1[0]+ coordPost1[0] - coordPt1[0], coordPrev1[1]+coordPost1[1]-coordPt1[1]};
  float coordBissect2[2] = {coordPrev2[0]+ coordPost2[0]-coordPt2[0], coordPrev2[1]+coordPost2[1]-coordPt2[1]};

  float semiAngle1A = (n1 == 1 || abs(angle1 - M_PI) < 1e-7) ? M_PI/2 : computeAngle(coordPt1, coordPrev1, coordBissect1);
  float semiAngle1B = (n1 == 1 || abs(angle1 - M_PI) < 1e-7) ? M_PI/2 : computeAngle(coordPt1, coordBissect1, coordPost1);
  float semiAngle2A = (n2 == 1 || abs(angle2 - M_PI) < 1e-7) ? M_PI/2 : computeAngle(coordPt2, coordPrev2, coordBissect2);
  float semiAngle2B = (n2 == 1 || abs(angle2 - M_PI) < 1e-7) ? M_PI/2 : computeAngle(coordPt2, coordBissect2, coordPost2);


  cout << " semi 1 = " << semiAngle1A << " VS " << semiAngle1B << endl;
  cout << " semi 2 = " << semiAngle2A << " VS " << semiAngle2B << endl;
  assert(abs(semiAngle1A - semiAngle1B) < 1e-6);
  assert(abs(semiAngle2A - semiAngle2B) < 1e-6);
  printCoords("Pt1 : ", coordPt1);
  printCoords("prev1 : ", coordPrev1);
  printCoords("post1 : ", coordPost1);
  printCoords("bissect1 : ", coordBissect1);

  printCoords("Pt2 : ", coordPt2);
  printCoords("prev2 : ", coordPrev2);
  printCoords("post2 : ", coordPost2);
  printCoords("bissect2 : ", coordBissect2);
  float angleMax1 = M_PI/2-semiAngle1A, angleMin1 = -angleMax1; //TODO calculer autrement et comparer
  float angleMax2 = M_PI/2-semiAngle2A, angleMin2 = -angleMax2; //TODO calculer autrement et comparer


  float step1 = (angleMax1-angleMin1)/bestAngleSampleFreq, step2 = (angleMax2-angleMin2)/bestAngleSampleFreq;
  float bestAnglePair[2] = {0,0};
  float bestScore = -1; //TODO
  float testAngle1 = angleMin1;

  size_t nComp12 = comp1.size()+comp2.size();
  std::vector<size_t> idsComps1And2;
  idsComps1And2.insert(idsComps1And2.begin(), comp1.cbegin(), comp1.cend());
  idsComps1And2.insert(idsComps1And2.end(), comp2.cbegin(), comp2.cend());
  cout << "The concatenation of the two comps is : ";
  for (int x : idsComps1And2)
    cout << x << " ";
  cout << endl;
  std::vector<std::vector<double>> origDistMatrix(n1+n2), newDistMatrix(n1+n2);
  for (size_t i = 0; i < n1+n2; i++)
  {
    newDistMatrix[i].resize(n1+n2, 0);
    origDistMatrix[i].resize(n1+n2);
    for (size_t j = 0; j < n1+n2; j++)
    {
      origDistMatrix[i][j] = distMatrix[idsComps1And2[i]][idsComps1And2[j]];
    }
  }
  std::vector<float> coords1Test(2*comp1.size()), coords2Test(2*comp2.size()); //TODO populate
  std::set<size_t>::const_iterator it = comp1.cbegin();
  for (size_t i = 0; it != comp1.cend(); ++it, i++)
  {
    size_t iPt = *it;
    coords1Test[2*i] = allCoords[2*iPt];
    coords1Test[2*i+1] = allCoords[2*iPt+1];
  }

  it = comp2.cbegin();
  for (size_t i = 0; it != comp2.cend(); ++it, i++)
  {
    size_t iPt = *it;
    coords2Test[2*i] = allCoords[2*iPt];
    coords2Test[2*i+1] = allCoords[2*iPt+1];
  }

  cout << "Angles for 1 (small) are to rotate min max " << angleMin1 << " et " << angleMax1 << endl;
  cout << "Angles for 2 (big) are to rotate min max " << angleMin2 << " et " << angleMax2 << endl;
  ttk::DistanceMatrixDistorsion ttkDistorsion;
  std::vector<double> foo;
  size_t nbIter1 = std::isfinite(step1) ? bestAngleSampleFreq : 1;
  size_t nbIter2 = std::isfinite(step2) ? bestAngleSampleFreq : 1;
  for (size_t i1 = 0; i1 < nbIter1; i1++, testAngle1+= step1)
  {
    cout << "\t\t\t\tTesting angle1 " << testAngle1 << endl;
    if (std::isfinite(step1))
      rotatePolygon(coords1Test, 2, coordPt1, step1);
    float testAngle2 = angleMin2;

    for (size_t i2 = 0; i2 < nbIter2; i2++, testAngle2 += step2)
    {
      cout << "\t\t\t\tTesting angle2 " << testAngle2 << endl;
      if (std::isfinite(step2))
        rotatePolygon(coords2Test, 2, coordPt2, step2);
      //TODO dans fonction...
      for (size_t i = 0; i < n1+n2; i++)
      {
        float coordARotate[2] = {i < n1 ? coords1Test[2*i] : coords2Test[2*(i-n1)], i < n1 ? coords1Test[2*i+1] : coords2Test[2*(n1-i)+1]};
        for (size_t j = 0; j < n1+n2; j++)
        {
          float coordBRotate[2] = {j < n1 ? coords1Test[2*j] : coords2Test[2*(j-n1)], j < n1 ? coords1Test[2*j+1] : coords2Test[2*(n1-j)+1]};
          //float coordBRotate[2] = j < n1 ? {coords1Test[2*j], coords1Test[2*j+1]} : {coords2Test[2*(j-n1)], coords2Test[2*(j-n1)+1]};
          float coord2Rotate[2] = {coords2Test[2*i], coords2Test[2*i+1]};
          newDistMatrix[i][j] = compute_dist(coordARotate, coordBRotate);
        }
        double curScore;
        ttkDistorsion.execute(origDistMatrix, newDistMatrix, curScore, foo);
        if (curScore > bestScore) // This pair of angle minimises the distorsion
        {
          bestScore = curScore;
          bestAnglePair[0] = std::isfinite(testAngle1) ? testAngle1 : 0;
          bestAnglePair[1] = std::isfinite(testAngle2) ? testAngle2 : 0;
          for (int i = 0; i < n1+n2; i++)
          {
            float coordRotated[2] = {i < n1 ? coords1Test[2*i] : coords2Test[2*(i-n1)], i < n1 ? coords1Test[2*i+1] : coords2Test[2*(n1-i)+1]};
            allCoords[2*idsComps1And2[i]] = coordRotated[0];
            allCoords[2*idsComps1And2[i]+1] = coordRotated[1];
          }
        }
      }

      if (step2 < 1e-7) // Angle is null, no need to divide
        break;
    }
    if (step1 < 1e-7)
      break;
  }

  cout << "The best angles are " << bestAnglePair[0] << " for 1 and " << bestAnglePair[1] << " for 2." << endl;
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
int ttk::TopologicalMapper::execute(std::vector<float> &inputPoints, float* outputCoords, const std::vector<std::vector<float>> &distMatrix) const
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

      // We retrive the current coordinates of our component.
      size_t cptCur = 0;
      for (int vertId : curComp)
      {
        for (int k = 0; k < dim; k++)
        {
          curPointsSet[dim*cptCur+k] = outputCoords[dim*vertId+k];
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
      float* coordEdgeVert = &inputPoints[idChosenVert];
      // We want to select, among all vertices in the convex hull, the one which is
      // closest to the vertex of the edge we work on.
      for (size_t vert : curHullVerts)
      {
        float* coordVert = &inputPoints[dim*vert];
        float dist2 = distMatrix[vert][idEdgeVert[idSet]];//  compute_dist2(coordVert, coordEdgeVert);
        std::cout << vert << " is " << dist2 << " close\n";
        if (dist2 < distMatrix[idChosenVert] [idEdgeVert[idSet]])//compute_dist2(&outputCoords[idChosenVert], coordEdgeVert))
          idChosenVert = vert;
      }
      cout << "We chose vertex " << idChosenVert << endl;
    }

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

    float unitCentreBigVect[2], unitCentreSmallVect[2];
    computeUnitVector(&outputCoords[idChosenBig*dim], coordsCentreBig, unitCentreBigVect);
    computeUnitVector(&outputCoords[idChosenSmall*dim], coordsCentreSmall, unitCentreSmallVect);
    printCoords("Le point choisit big : ", &outputCoords[idChosenBig*dim]);
    printCoords("Barycentre for big comp is : ", coordsCentreBig);
    printCoords("Unit vect to centre big : ", unitCentreBigVect);

    //TODO create vector from two points?
    float goalCoordChosenSmall[2] = {outputCoords[idChosenBig*dim]-edgeCost*unitCentreBigVect[0], outputCoords[idChosenBig*dim+1]-edgeCost*unitCentreBigVect[1]};
    if (sizeBigHull == 1)
    {
      goalCoordChosenSmall[0] = outputCoords[idChosenBig*dim]+edgeCost;
      goalCoordChosenSmall[1] = outputCoords[idChosenBig*dim+1];
      //continue;
    }
    float smallCompMoveVect[2] = {goalCoordChosenSmall[0] - outputCoords[idChosenSmall*dim], goalCoordChosenSmall[1] - outputCoords[idChosenSmall*dim+1]};

    float distBaryPointSmall = compute_dist(&outputCoords[idChosenSmall*dim], coordsCentreSmall);
    printCoords("Le point choisi small : ", &outputCoords[idChosenSmall*dim]);
    printCoords("Barycentre for small comp is : ", coordsCentreSmall);
    printCoords("Unit vect to centre small : ", unitCentreSmallVect);
    cout << " and the distBaryPointSmall = " << distBaryPointSmall << endl;
    float finalPosBarySmall[2] = {goalCoordChosenSmall[0]+unitCentreSmallVect[0]*distBaryPointSmall, goalCoordChosenSmall[1]+unitCentreSmallVect[1]*distBaryPointSmall};

    //TODO : actuellement pour aligner les barycentres, normalement on voudrait aligner les bissectrices je pense
    float rotationAngle = computeAngle(goalCoordChosenSmall, coordsCentreBig, finalPosBarySmall, dim==3);
    std::cout << " We chose for big and small: " << idChosenBig << " and " << idChosenSmall << std::endl;

    if (nSmall > 1 && std::isfinite(rotationAngle))
    {
      std::cout << "angle is " << rotationAngle << std::endl;
      cout << "cos,sin: " << cos(rotationAngle) << "," << sin(rotationAngle) << endl;
      printCoords("coordCenter = ", goalCoordChosenSmall);
      printCoords("coordFrom = ", coordsCentreSmall);
      printCoords("coordTo = ", finalPosBarySmall);
    }
    for (size_t curIdSmall : compSmall)
    {
      std::cout << "pointSmall : " << outputCoords[curIdSmall*dim] << "," << outputCoords[curIdSmall*dim+1] << endl;
      if (curIdSmall == idChosenSmall)
        cout << "\t\t=> CHOSEN small" << endl;
      outputCoords[curIdSmall*dim] += smallCompMoveVect[0];
      outputCoords[curIdSmall*dim+1] += smallCompMoveVect[1];
      std::cout << "\tPre-new-coordinates for " << curIdSmall << " are: " << outputCoords[curIdSmall*dim]<<","<<outputCoords[curIdSmall*dim+1] << "\n";
      if (nSmall > 1 && std::isfinite(rotationAngle))
      {
        rotate(&outputCoords[curIdSmall*dim], goalCoordChosenSmall, rotationAngle); //Ou - rotatinAngle ?
      std::cout << "\t\tPost-new-coordinates for " << curIdSmall << " are: " << outputCoords[curIdSmall*dim]<<","<<outputCoords[curIdSmall*dim+1] << "\n";
      }
    }

    if (nBig > 1)
      rotateMergingCompsBest(idsInHullSmall, idsInHullBig, compSmall, compBig, idChosenSmall, idChosenBig, distMatrix, outputCoords);
    std::cout << "Unit vector bary point = " << unitCentreBigVect[0] << "," << unitCentreBigVect[1] << endl;
    std::cout << "Goal coordinates = " << goalCoordChosenSmall[0] << "," << goalCoordChosenSmall[1] << endl;
    std::cout << "Moving vector coordinates = " << smallCompMoveVect[0] << "," << smallCompMoveVect[1] << endl;

    /*
    for (size_t curIdBig : compBig)
    {
      std::cout << "pointBig : " << outputCoords[curIdBig*dim] << "," << outputCoords[curIdBig*dim+1] << endl;
      if (curIdBig == idChosenBig)
        cout << "\t\t=> CHOSEN big" << endl;
    }

    cout << "barycenterBig = " << coordsCentreBig[0] << "," << coordsCentreBig[1] << endl;
    */
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

