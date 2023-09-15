#include <TopologicalMapper.h>
#include <vector>
#include <set>
#include <queue>
#include <map>
#include <utility>
#include <algorithm>
#include <Geometry.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>


#include <cmath>

namespace bg = boost::geometry;
BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)
typedef boost::tuple<double, double> Point;
typedef boost::geometry::model::polygon<Point> Polygon;
typedef boost::geometry::model::multi_point<Point> Mpoints;

ttk::TopologicalMapper::TopologicalMapper() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("TopologicalMapper");
}
ttk::TopologicalMapper::~TopologicalMapper() = default;

// Normalize a given vector.
void computeUnitVector(double* const coordOrig, double* const coordDest, double* const coordVect)
{
  double tmp[2] = {coordDest[0] - coordOrig[0], coordDest[1] - coordOrig[1]};
  double dist = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1]);
  if (dist < EPS)
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

double computeAngle(double const ptA[], double const ptB[], double const ptC[])
{
  double angle;
  double vect1[2] = {ptB[0]-ptA[0], ptB[1]-ptA[1]};
  double vect2[2] = {ptC[0]-ptA[0], ptC[1]-ptA[1]};
  double dirVect[2] = {0,0};

#if VERB > 4
  printCoords("A : ", ptA);
  printCoords("B : ", ptB);
  printCoords("C : ", ptC);
#endif
  angle = atan2(vect2[1], vect2[0]) - atan2(vect1[1], vect1[0]);
  if (angle < 0)
    angle += 2 * M_PI;
  return angle;
}

inline double computeSquaredDistBetweenMatrices(const std::vector<std::vector<double>> &mat1, const std::vector<std::vector<double>> &mat2)
{
  double ret = 0;
  size_t n = mat1.size(), m = mat1[0].size();
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < m; j++)
      ret += (mat1[i][j]-mat2[i][j])*(mat1[i][j]-mat2[i][j]);

  return ret;
}


inline double deg(double angle)
{
  return (angle*180)/M_PI;
}

// Rotate the set of coords by the given angle, considering the given center as center of the rotation.
void rotatePolygon(std::vector<double> &coords, double* centerCoords, const double angle)
{
  double xCenter = centerCoords[0], yCenter = centerCoords[1];
  size_t nbPoint = coords.size()/2;
#if VERB > 3
  std::cout << "center = " << xCenter << "," << yCenter << std::endl;
#endif
  for (size_t iPt = 0; iPt < nbPoint; iPt++)
  {
    double &x = coords[iPt*2], &y = coords[iPt*2+1];
    double xNew, yNew;
    xNew = (x-xCenter)*cos(angle)-(y-yCenter)*sin(angle)+xCenter;
    yNew = (y-yCenter)*cos(angle)+(x-xCenter)*sin(angle)+yCenter;
#if VERB > 3
    std::cout << "point = " << x << "," << y << " ===> (" << xNew << ", " << yNew << ")" <<  std::endl;
#endif
    x = xNew;
    y = yNew;
  }
}

void rotate(double ptToRotate[], double const centre[], double angle)
{
  const double &xCtr = centre[0], &yCtr = centre[1];
  double &xPt = ptToRotate[0], &yPt = ptToRotate[1];
  const double dx = xPt-xCtr, dy = yPt-yCtr;
  xPt = dx*cos(angle)-dy*sin(angle)+xCtr;
  yPt = dx*sin(angle)+dy*cos(angle)+yCtr;
}

// Fills its last two arguments with the position of the end of the previous (resp. next) edge so that idCenter is the point at the angle at the intersection of the previous and next edges.
void getPrevPostEdges(const std::vector<size_t> &idsPtsPolygon, size_t idCenter, double* const allCoords, double coordPrev[2], double coordPost[2])
{
  size_t n = idsPtsPolygon.size();
  size_t iPtPrev, iPtPost;

#if VERB > 4
  for (int i = 0; i < n; i++)
    printCoords("", &allCoords[idsPtsPolygon[i]*2]);
#endif
  for (size_t i = 0; i < n; i++)
  {
    if (idsPtsPolygon[i] == idCenter)
    {
      iPtPost = idsPtsPolygon[(i+1)%n];
      iPtPrev = idsPtsPolygon[(i+n-1)%n];
      break;
    }
  }

  coordPrev[0] = allCoords[2*iPtPrev];
  coordPrev[1] = allCoords[2*iPtPrev+1];
  coordPost[0] = allCoords[2*iPtPost];
  coordPost[1] = allCoords[2*iPtPost+1];
}

// Tries to find the best angle of rotation for the two components. Updates the coordiates of their vertices accordingly.
void rotateMergingCompsBest(const std::vector<size_t> &hull1, const std::vector<size_t> &hull2, const std::set<size_t> &comp1, const std::set<size_t> &comp2, size_t iPt1, size_t iPt2, const std::vector<std::vector<float>> &distMatrix, double* allCoords, size_t angleSamplingFreq, size_t nThread)
{
  TTK_FORCE_USE(nThread);
  // The distance between the two components.
  double shortestDistPossible = compute_dist(&allCoords[2*iPt1], &allCoords[2*iPt2]);
  double coordPt1[2] = {allCoords[2*iPt1], allCoords[2*iPt1+1]};
  double coordPt2[2] = {allCoords[2*iPt2], allCoords[2*iPt2+1]};
  size_t hull1Size = hull1.size(), hull2Size = hull2.size();
  size_t comp1Size = comp1.size(), comp2Size = comp2.size();

  double coordPrev1[2], coordPost1[2];
  double coordPrev2[2], coordPost2[2];
  getPrevPostEdges(hull1, iPt1, allCoords, coordPrev1, coordPost1);
  getPrevPostEdges(hull2, iPt2, allCoords, coordPrev2, coordPost2);

#if VERB > 3
  std::cout << "Comp 1 = ";
  for (int x : comp1)
    std::cout << x << " ";
  std::cout << std::endl;
  std::cout << "Comp 2 = ";
  for (int x : comp2)
    std::cout << x << " ";
  std::cout << std::endl;
  std::cout << "prev2, cur2, post2 : " << iPtPrev2 << ", " << iPt2 << ", " << iPtPost2 << std::endl;
#endif

  double angle1 = computeAngle(coordPt1, coordPrev1, coordPost1);
  double angle2 = computeAngle(coordPt2, coordPrev2, coordPost2);
#if VERB > 2
  std::cout << "The angles are " << deg(angle1) << " " << deg(angle2) << std::endl;
#endif
  if (angle1 > M_PI)
  {
    std::swap(coordPrev1, coordPost1);
    angle1 = 2*M_PI-angle1;
  }
  if (angle2 > M_PI)
  {
    std::swap(coordPrev2, coordPost2);
    angle2 = 2*M_PI-angle2;
  }
#if VERB > 2
  std::cout << "The angles REALLY are " << deg(angle1) << " " << deg(angle2) << std::endl;
#endif

  double coordBissect1[2] = {coordPrev1[0], coordPrev1[1]};
  double coordBissect2[2] = {coordPrev2[0], coordPrev2[1]};
  rotate(coordBissect1, coordPt1, angle1/2);
  rotate(coordBissect2, coordPt2, angle2/2);

  double semiAngle1 = computeAngle(coordPt1, coordPrev1, coordBissect1);
  double semiAngle2 = computeAngle(coordPt2, coordPrev2, coordBissect2);

  if (abs(semiAngle1 - angle1/2) > EPS)
    std::cout << "çavapas angle1 et semi1 : " << semiAngle1 <<  " VS " << angle1/2 << "\n";
  if (abs(semiAngle2 - angle2/2) > EPS)
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
  double step1 = (angleMax1-angleMin1)/angleSamplingFreq, step2 = (angleMax2-angleMin2)/angleSamplingFreq;
  double bestAnglePair[2] = {0,0};
  double bestScore = 1e34;

  std::vector<size_t> idsComp1, idsComp2;;
  idsComp1.insert(idsComp1.begin(), comp1.cbegin(), comp1.cend());
  idsComp2.insert(idsComp2.begin(), comp2.cbegin(), comp2.cend());
  std::vector<std::vector<double>> origDistMatrix(comp1Size);
  for (size_t i = 0; i < comp1Size; i++)
  {
    origDistMatrix[i].resize(comp2Size);
    for (size_t j = 0; j < comp2Size; j++)
    {
      origDistMatrix[i][j] = distMatrix[idsComp1[i]][idsComp2[j]];
    }
  }
  std::vector<double> initialCoords1(2*comp1Size), initialCoords2(2*comp2Size);
  std::set<size_t>::const_iterator it = comp1.cbegin();
  for (size_t i = 0; it != comp1.cend(); ++it, i++)
  {
    size_t iPt = *it;
    initialCoords1[2*i] = allCoords[2*iPt];
    initialCoords1[2*i+1] = allCoords[2*iPt+1];
  }

  it = comp2.cbegin();
  for (size_t i = 0; it != comp2.cend(); ++it, i++)
  {
    size_t iPt = *it;
    initialCoords2[2*i] = allCoords[2*iPt];
    initialCoords2[2*i+1] = allCoords[2*iPt+1];
  }


#if VERB > 0
  std::cout << "Angles for 1 (small) are to rotate min max " << deg(angleMin1) << " et " << deg(angleMax1) << std::endl;
  std::cout << "Angles for 2 (big) are to rotate min max " << deg(angleMin2) << " et " << deg(angleMax2) << std::endl;
#endif
  size_t nbIter1 = std::isfinite(step1) ? angleSamplingFreq+1 : 1;
  size_t nbIter2 = std::isfinite(step2) ? angleSamplingFreq+1 : 1;

  if (step1*nbIter1 < 0.001) // No need to split such a small angle
    nbIter1 = 1;
  if (step2*nbIter2 < 0.001) // No need to split such a small angle
    nbIter2 = 1;
#pragma omp parallel for num_threads(nThread) shared(allCoords)
  for (size_t i1 = 0; i1 < nbIter1; i1++)
  {
    std::vector<double> coords1Test(2*comp1Size), coords2Test(2*comp2Size);
    coords1Test = initialCoords1;
    double testAngle1 = angleMin1+step1*i1;
    rotatePolygon(coords1Test, coordPt1, testAngle1);
#if VERB > 2
    std::cout << "\t\t\t\tTesting angle1 " << deg(testAngle1) << std::endl;
#endif
#if CHECK
    if (comp1Size >= 2)
    {
      double *A1, *B1, *C1;
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
      if (abs(deg(testAngle1) - deg(totoAngle)) > 0.5)
        std::cout << "ERROR1 : " << deg(testAngle1) << " VS " << deg(totoAngle) << "\n";
    }
#endif

    for (size_t i2 = 0; i2 < nbIter2; i2++)
    {
      coords2Test = initialCoords2;
#if VERB > 2
      std::cout << "\t\t\t\tTesting angle2 " << deg(testAngle2) << std::endl;
#endif
      double testAngle2 = angleMin2+i2*step2;
      rotatePolygon(coords2Test, coordPt2, testAngle2);

      double curScore = 0;
      //TODO dans fonction...
      for (size_t i = 0; i < comp1Size; i++)
      {
        double coordARotate[2] = {coords1Test[2*i],coords1Test[2*i+1]};
        for (size_t j = 0; j < comp2Size; j++)
        {
          double coordBRotate[2] = {coords2Test[2*j],coords2Test[2*j+1]};
          double newDist = compute_dist(coordARotate, coordBRotate);
          curScore += abs(newDist-origDistMatrix[i][j])*(newDist-origDistMatrix[i][j]);
          if (newDist+EPS < shortestDistPossible)
          {
            std::cout << "problem " << newDist << '(' << idsComp1[i] << ',' << idsComp2[j] << ')' << " is lower than " << shortestDistPossible << " =====> (" << i1 << "," << i2 << ") (angles iterations)" << std::endl;
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
      //double curScore = computeSquaredDistBetweenMatrices(origDistMatrix, newDistMatrix);
#if CHECK
      if (comp2Size >= 2 && comp1Size >= 1)
      {
        double *A, *B, *C;
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
        if (abs(deg(testAngle2) - deg(totoAngle)) > 0.5)
          std::cout << "ERROR2 : " << deg(testAngle2) << " VS " << deg(totoAngle) << "\n";
      }
#endif

#pragma omp critical
{
  //TODO maybe juste update bestScore et bestAnglePair, calculer coords à la fin ?
      if (curScore < bestScore) // This pair of angle minimises the distortion
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

  for (size_t i1 : idsComp1)
    rotate(&allCoords[2*i1], coordPt1, bestAnglePair[0]);

  for (size_t i2 : idsComp2)
    rotate(&allCoords[2*i2], coordPt2, bestAnglePair[1]);

#if VERB > 1
  std::cout << "best score is " << bestScore << std::endl;

  std::cout << "The best angles are " << deg(bestAnglePair[0]) << " for 1 and " << deg(bestAnglePair[1]) << " for 2." << std::endl;
#endif
}

void getConvexHull(const std::vector<double>& coords, size_t dim, std::vector<size_t> &idsInHull)
{
  //TODO copier coords puis sort
  size_t nbPoint = coords.size()/dim;
  if (nbPoint <= 2)
  {
    idsInHull.push_back(0);
    if (nbPoint == 2)
    {
      double dist = compute_dist(&coords[0], &coords[2]);

      if (dist > EPS) //TODO voir si distance trop grande ?
        idsInHull.push_back(1);
    }
    return;
  }

  Mpoints multiPoints;
  Polygon hull;
  for (size_t i = 0; i < nbPoint; i++)
  {
    boost::geometry::append(multiPoints, Point(coords[2*i], coords[2*i+1]));
  }
  boost::geometry::convex_hull(multiPoints, hull);
  for (auto boostPt : hull.outer())
  {
    double coordsCur[2] = {boostPt.get<0>(), boostPt.get<1>()};
    for (int j = 0; j < coords.size()/2; j++)
    {
      if (abs(coords[2*j]-coordsCur[0])+abs(coords[2*j+1]-coordsCur[1]) < EPS)
      {
        idsInHull.push_back(j);
      }
    }
  }
  // Boost closes the polygon, hence the first and the last vertices are the identical.
  idsInHull.pop_back();
}
