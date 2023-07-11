/// \ingroup vtk
/// \class ttkMapper
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date May 2022.
///
/// \brief TTK Mapper VTK-filter.
///
/// This class generates a mapper from a data-set with a given number
/// of buckets. Data-set vertices and edges are placed into buckets
/// from their input scalar field value then connected components are
/// detected per bucket.
///
/// \param Input vtkDataSet
/// \param Output Nodes, Arcs (both vtkUnstructuredGrid), Segmentation
/// (vtkDataSet)
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).

#pragma once

// VTK includes

// VTK Module
#include <ttkMapperModule.h>

// ttk code includes
#include <Mapper.h>
#include <ttkAlgorithm.h>
#include <ttkMacros.h>
#include <vtkNew.h>
#include <vtkIntArray.h>
#include <vtkUnstructuredGrid.h>

class TTKMAPPER_EXPORT ttkMapper : public ttkAlgorithm, protected ttk::Mapper {
public:
  static ttkMapper *New();

  vtkTypeMacro(ttkMapper, ttkAlgorithm);

  void SetNumberOfBuckets(int nbBuckets)
  {
    this->NumberOfBuckets = nbBuckets;
    this->firstTimeReembed_ = true;
    this->Modified();
  }
  vtkGetMacro(NumberOfBuckets, int);

  void SetReEmbedMapper(bool enable)
  {
    this->ReEmbedMapper = enable;
    this->firstTimeReembed_ = true;
    this->Modified();
  }
  vtkGetMacro(ReEmbedMapper, bool);

  void SetSelectMatrixWithRegexp (bool enable)
  {
    this->SelectMatrixWithRegexp = enable;
    this->firstTimeReembed_ = true;
    this->Modified();
  }
  vtkGetMacro(SelectMatrixWithRegexp, bool);

  void SetDistanceMatrixRegexp(const std::string &matName)
  {
    this->DistanceMatrixRegexp = matName;
    this->firstTimeReembed_ = true;
    this->Modified();
  }
  vtkGetMacro(DistanceMatrixRegexp, std::string);

  void SetDistMat(const std::string &s) {
    this->DistanceMatrix.emplace_back(s);
    this->Modified();
  }
  void ClearDistMat() {
    this->DistanceMatrix.clear();
    this->Modified();
  }

  void SetLowerDimension(unsigned int lowDim)//(LOWER_DIMENSION lowDim)
  {
    this->LowerDimension = (LOWER_DIMENSION)lowDim;
    this->firstTimeReembed_ = true;
    this->Modified();
  }
  vtkGetEnumMacro(LowerDimension, LOWER_DIMENSION);

  void SetReductionAlgo(int redAlgo)//(REDUCTION_ALGO redAlgo)
  {
    this->ReductionAlgo = (REDUCTION_ALGO)redAlgo;
    this->firstTimeReembed_ = true;
    this->Modified();
  }
  //ttkSetEnumMacro(ReductionAlgo, REDUCTION_ALGO);
  vtkGetEnumMacro(ReductionAlgo, REDUCTION_ALGO);

  void SetReembedMethod(int reembedMethod)//REEMBED_METHOD reembedMethod)
  {
    this->ReembedMethod = (REEMBED_METHOD)reembedMethod;
    this->firstTimeReembed_ = true;
    this->Modified();
  }
  vtkGetEnumMacro(ReembedMethod, REEMBED_METHOD);

  vtkGetMacro(DilatationCoeff, double);
  void SetDilatationCoeff(double coeff)
  {
    this->DilatationCoeff = coeff;
    ModifiedDilatation();
  }
  
  void ModifiedDilatation() {
    if (!firstTimeReembed_)
    {
      std::cerr <<"Modifying, yeah :-)";
      needPartialUpdate_ = DILATATION;
      ttkAlgorithm::Modified(); //Useless?
    }
  }


  vtkGetMacro(AlphaCoeff, double);
  void SetAlphaCoeff(double coeff)
  {
    this->AlphaCoeff = coeff;
    ModifiedAlpha();
  }
  
  void ModifiedAlpha() {
    if (!firstTimeReembed_)
    {
      std::cerr <<"Modifying, yeah :-)";
      needPartialUpdate_ = ALPHA;
      ttkAlgorithm::Modified(); //Useless?
    }
    else if (this->AlphaCoeff < 1e-6)
    {
      needPartialUpdate_ = NOPE;
      firstTimeReembed_ = true;
    }
  }

protected:
  ttkMapper();

  ~ttkMapper() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  // domain vertices coordinates in high dimension
  enum NEED_UPDATE{NOPE, DILATATION, ALPHA};
  bool SelectMatrixWithRegexp{false};
  std::string DistanceMatrixRegexp{".*"};
  std::vector<std::string> DistanceMatrix{};
  int needPartialUpdate_{NOPE}; // To avoid to recompute everything when not necessary.
  bool firstTimeReembed_{true}; // To avoid to recompute everything when not necessary. //Check sze f arrays...
  std::vector<float> pointsCoordsBackup_;
  vtkNew<vtkIntArray> connCompPrev_{}, bucketPrev_{};
  vtkNew<vtkUnstructuredGrid> arcsPrev_{}, nodesPrev_{};
};
