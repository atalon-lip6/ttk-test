/// \defgroup examples examples
/// \brief The Topology ToolKit - Example programs.
/// @{
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date October 2017.
///
/// \brief Minimalist VTK-based TTK example pipeline, including:
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

#include <CommandLineParser.h>
#include <vector>
#include <string>

#include <ttkMapper.h>
#include <Mapper.h>
#include <ttkPersistenceCurve.h>
#include <ttkPersistenceDiagram.h>
#include <ttkTopologicalSimplification.h>

#include <vtkNew.h>
#include <vtkTableWriter.h>
#include <vtkThreshold.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

using namespace std;

int main(int argc, char **argv) {

  ttk::CommandLineParser parser;
  ttk::globalDebugLevel_ = 2;

  std::string inputFilePath;

  parser.setArgument("i", &inputFilePath, "Path to input VTU file");
  parser.parse(argc, argv);

  // 1. loading the input data
  vtkNew<vtkXMLUnstructuredGridReader> reader{};
  reader->SetFileName(inputFilePath.data());

  // 2. computing the persistence diagram
  vtkNew<ttkMapper> mapper{};
  mapper->SetInputConnection(reader->GetOutputPort());
  mapper->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "Elevation");

  mapper->SetSelectMatrixWithRegexp(true);
  mapper->SetDistanceMatrixRegexp("Distance.*");

  ttk::Mapper::Matrix inputDistMat{};
  mapper->SetNumberOfBuckets(500);
  mapper->SetReEmbedMapper(true);
  vtkNew<vtkXMLUnstructuredGridWriter> mapperWriter{};
  mapperWriter->SetInputConnection(mapper->GetOutputPort());
  mapperWriter->SetFileName("mapper.vtu");
  mapperWriter->Write();
  std::cout << "YEAH\n";

  return 0;
}

/// @}
