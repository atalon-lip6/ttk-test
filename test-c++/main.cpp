/// \ingroup examples
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date October 2017.
///
/// \brief Minimalist C++-only TTK example pipeline, including:
///  -# The computation of a persistence curve
///  -# The computation of a persistence diagram
///  -# The selection of the most persistent pairs of the diagram
///  -# The pre-simplification of the data according to this selection
///  -# The computation of the Morse-Smale complex on this simplified data
///  -# The storage of the output of this pipeline to disk.
///
/// This example reproduces the Figure 1 of the TTK companion paper:
/// "The Topology ToolKit", J. Tierny, G. Favelier, J. Levine, C. Gueunet, M.
/// Michaux., IEEE Transactions on Visualization and Computer Graphics, Proc.
/// of IEEE VIS 2017.
///
/// See the individual VTK wrappers (core/vtk/) to see how to use each ttk::base
/// (C++-only) TTK component.

// include the local headers
#include <TopologyToolKit.h>

#include <cassert>
#include <cmath>
#include <iostream>

int nbCoord = 4;
int trigSize = 2;

int load(const std::string &inputPath,
         std::vector<float> &pointSet,
         std::vector<long long int> &triangleSetCo,
         std::vector<long long int> &triangleSetOff) {

  // load some terrain from some OFF file.

  if(inputPath.empty())
    return -1;

  ttk::Debug dbg;
  dbg.setDebugLevel(ttk::globalDebugLevel_);
  dbg.setDebugMsgPrefix("main::load");
  dbg.printMsg("Reading input mesh...");

  int vertexNumber = 0, triangleNumber = 0;
  std::string keyword;

  std::ifstream f(inputPath.data(), std::ios::in);
  std::string foo;

  if(!f) {
    dbg.printErr("Cannot read file `" + inputPath + "'!");
    return -1;
  }

  f >> keyword;

  if(keyword != "OFF") {
    dbg.printErr("Input OFF file `" + inputPath + "' seems invalid :(");
    return -2;
  }

  f >> vertexNumber;
  f >> triangleNumber;
  f >> keyword;
  std::cout << vertexNumber << " " << triangleNumber << " " << keyword << std::endl;

  pointSet.resize(nbCoord * vertexNumber);
  triangleSetCo.resize(trigSize * triangleNumber);
  triangleSetOff.resize(triangleNumber + 1);

  for(int i = 0; i < nbCoord * vertexNumber; i++) {
    f >> pointSet[i];
  }

  //std::cout << "LAST = " << pointSet[nbCoord*vertexNumber-1] << "\n";
  int offId = 0;
  int coId = 0;
  for(int i = 0; i < triangleNumber; i++) {
    int cellSize;
    f >> cellSize;
    if(cellSize != trigSize) {
      std::cerr << "cell size " << cellSize << " != " << trigSize << std::endl;
      return -3;
    }
    triangleSetOff[offId++] = coId;
    for(int j = 0; j < trigSize; j++) {
      int cellId;
      f >> cellId;
      triangleSetCo[coId++] = cellId;
    }
    getline(f, foo);
  }
  triangleSetOff[offId] = coId; // the last one

  f.close();

  dbg.printMsg("... done! (read " + std::to_string(vertexNumber) + " vertices, "
               + std::to_string(triangleNumber) + " triangles)");

  return 0;
}

int save(const std::vector<float> &pointSet,
         const std::vector<long long int> &triangleSetCo,
         const std::vector<long long int> &triangleSetOff,
         const std::string &outputPath) {

  // save the simplified terrain in some OFF file
  std::string fileName(outputPath);

  std::ofstream f(fileName.data(), std::ios::out);

  if(!f) {
    ttk::Debug dbg;
    dbg.setDebugLevel(ttk::globalDebugLevel_);
    dbg.setDebugMsgPrefix("main::save");
    dbg.printErr("Could not write output file `" + fileName + "'!");
    return -1;
  }

  const int nbTriangles = triangleSetOff.size() - 1;

  f << "OFF" << std::endl;
  f << pointSet.size() / 3 << " " << nbTriangles << " 0" << std::endl;


  for(int i = 0; i < (int)pointSet.size() / 3; i++) {
    for(int j = 0; j < 3; j++) {
      f << pointSet[3 * i + j];
      f << " ";
    }
    f << std::endl;
  }

  for(int i = 0; i < nbTriangles; i++) {
    int cellSize = triangleSetOff[i + 1] - triangleSetOff[i];
    assert(cellSize == 3);
    f << cellSize << " ";
    for(int j = triangleSetOff[i]; j < triangleSetOff[i + 1]; j++) {
      f << triangleSetCo[j];
      f << " ";
    }
    f << std::endl;
  }

  f.close();

  return 0;
}

int main(int argc, char **argv) {

  std::string inputFilePath;
  ttk::CommandLineParser parser;

  ttk::globalDebugLevel_ = 3;

  // register the arguments to the command line parser
  parser.setArgument("i", &inputFilePath, "Path to input OFF file");
  // parse
  parser.parse(argc, argv);

  std::vector<float> pointSet;
  std::vector<long long int> triangleSetCo, triangleSetOff;
  ttk::Triangulation triangulation;

  // load the input
  load(inputFilePath, pointSet, triangleSetCo, triangleSetOff);
  std::vector<float> pointXYZ;

  for (int i = 0; i < pointSet.size()/nbCoord; i++)
  {
      pointXYZ.push_back(pointSet[nbCoord*i]);
      pointXYZ.push_back(pointSet[nbCoord*i+1]);
      pointXYZ.push_back(pointSet[nbCoord*i+2]);
  }

  std::cout << "avant trig" << std::endl;
  triangulation.setInputPoints(pointXYZ.size() / 3, pointXYZ.data());
  std::cout << "Il y a " << pointXYZ.size()/3 << " points,  chacun avec " << pointXYZ.size()/(pointSet.size()/nbCoord) << " comps " << std::endl;
  long long int triangleNumber = triangleSetOff.size() - 1;
#ifdef TTK_CELL_ARRAY_NEW
  std::cout << "LOL" << std::endl;
  triangulation.setInputCells(
    triangleNumber, triangleSetCo.data(), triangleSetOff.data());
#else
  LongSimplexId *triangleSet;
  CellArray::TranslateToFlatLayout(triangleSetCo, triangleSetOff, triangleSet);
  triangulation.setInputCells(triangleNumber, triangleSet);
#endif
  std::cout << "après trig" << std::endl;

  // NOW, do the TTK processing

  // computing some elevation
  std::vector<float> height(pointSet.size() / nbCoord);
  int vertexId = 0;
  // use the z-coordinate here
  for(int i = nbCoord-1; i < (int)pointSet.size(); i += nbCoord) {
    height[vertexId] = pointSet[i];
    vertexId++;
  }

  // order array: every vertex sorted according to the elevation field
  std::vector<ttk::SimplexId> order(height.size());
  // precondition/fill in the order array according to the elevation field
  ttk::preconditionOrderArray(height.size(), height.data(), order.data());

  /*
  // 2. computing the persistence diagram
  ttk::PersistenceDiagram diagram;
  std::vector<ttk::PersistencePair> diagramOutput;
  diagram.preconditionTriangulation(&triangulation);
  diagram.execute(
    diagramOutput, height.data(), 0, order.data(), &triangulation);

  // 3. computing the persistence curve from the persistence diagram
  ttk::PersistenceCurve curve;
  std::array<ttk::PersistenceCurve::PlotType, 4> outputCurve;
  curve.execute(outputCurve, diagramOutput);

  // 4. selecting the critical point pairs
  std::vector<float> simplifiedHeight = height;
  std::vector<ttk::SimplexId> authorizedCriticalPoints, simplifiedOrder = order;
  for(int i = 0; i < (int)diagramOutput.size(); i++) {
    if(diagramOutput[i].persistence() > 0.05) {
      // 5. selecting the most persistent pairs
      authorizedCriticalPoints.push_back(diagramOutput[i].birth.id);
      authorizedCriticalPoints.push_back(diagramOutput[i].death.id);
    }
  }

  // 6. simplifying the input data to remove non-persistent pairs
  ttk::TopologicalSimplification simplification;
  simplification.preconditionTriangulation(&triangulation);
  simplification.execute<float>(
    height.data(), simplifiedHeight.data(), authorizedCriticalPoints.data(),
    order.data(), simplifiedOrder.data(), authorizedCriticalPoints.size(),
    false, triangulation);

  // assign the simplified values to the input mesh
  for(int i = 0; i < (int)simplifiedHeight.size(); i++) {
    pointSet[3 * i + 2] = simplifiedHeight[i];
  }

  // 7. computing the Morse-Smale complex
  ttk::MorseSmaleComplex morseSmaleComplex;
  // critical points
  ttk::MorseSmaleComplex::OutputCriticalPoints outCriticalPoints{};
  // 1-separatrices
  ttk::MorseSmaleComplex::Output1Separatrices out1Separatrices{};
  // 2-separatrices
  ttk::MorseSmaleComplex::Output2Separatrices out2Separatrices{};
  // segmentation
  std::vector<ttk::SimplexId> ascendingSegmentation(
    triangulation.getNumberOfVertices(), -1),
    descendingSegmentation(triangulation.getNumberOfVertices(), -1),
    mscSegmentation(triangulation.getNumberOfVertices(), -1);
  ttk::MorseSmaleComplex::OutputManifold outSegmentation{
    ascendingSegmentation.data(), descendingSegmentation.data(),
    mscSegmentation.data()};

  morseSmaleComplex.preconditionTriangulation(&triangulation);
  morseSmaleComplex.execute(
    outCriticalPoints, out1Separatrices, out2Separatrices, outSegmentation,
    simplifiedHeight.data(), 0, simplifiedOrder.data(), triangulation);

  // save the output
  save(pointSet, triangleSetCo, triangleSetOff, "output.off");
*/

  // 8. Compute a mapper
  size_t nPoint = pointSet.size()/nbCoord;
  ttk::Mapper mapper;
  mapper.enableReEmbed();
  std::vector<int> outputBucket(nPoint), outputConnComp(nPoint), connCompBucket;
  std::vector<std::array<float, 3>> compBaryCoords;
  std::vector<std::set<SimplexId>> compArcs;
  std::vector<float> outputPointsCoords(nPoint*3);

  ttk::Mapper::Matrix distMat(nPoint, nPoint);

  for (size_t i1 = 0; i1 < nPoint; i1++)
  {
    for (size_t i2 = 0; i2 < nPoint; i2++)
    {
      float dx = (pointSet[3*i1]-pointSet[3*i2]);
      float dy = (pointSet[3*i1+1]-pointSet[3*i2+1]);
      float dz = (pointSet[3*i1+2]-pointSet[3*i2+2]);
      distMat.get(i1, i2) = sqrt(dx*dx+dy*dy+dz*dz);
    }
  }
  std::cerr << "Preconditioning triangulation" << std::endl;
  mapper.preconditionTriangulation(&triangulation);
  std::cerr << "Finishing preconditioning" << std::endl;
  mapper.execute(outputBucket.data(), outputConnComp.data(), compBaryCoords, connCompBucket, compArcs, outputPointsCoords.data(), distMat, height.data(), triangulation);


  return 0;
}
