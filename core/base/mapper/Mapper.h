// \ingroup base
/// \class ttk::Mapper
/// \author (mainly) Pierre Guillou <pierre.guillou@lip6.fr>
/// \author (a bit) Alexandre Talon <alexandre.talon@lip6.fr>
/// \date March 2023.
///
/// \brief TTK processing package for mapper.
///
/// This class generates a mapper from a data-set with a given number
/// of buckets. Data-set vertices and edges are placed into buckets
/// from their input scalar field value then connected components are
/// detected per bucket.
///
/// \param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \sa ttkMapper.cpp %for a usage example.

#pragma once

#include <DimensionReduction.h>
#include <LDistanceMatrix.h>
#include <Triangulation.h>
#include <UnionFind.h>

#include <array>
#include <limits>
#include <set>

inline double max(double &a, double &b) {
  if(a < b)
    return b;
  return a;
}

namespace ttk {

  class Mapper : virtual public Debug {
  public:
    Mapper();
    ~Mapper() { printErr("Destroying :'(");}
    enum class LOWER_DIMENSION {
      LOWER_DIM_2D = 2,
      LOWER_DIM_3D = 3,
    };
    using REDUCTION_ALGO = DimensionReduction::METHOD;

    enum class REEMBED_METHOD {
      ARCS_GEODESIC,
      BUCKETS_GEODESIC,
    };
    /**
     * @brief Utility class to store matrices/2-dimensional arrays
     * inside a single buffer
     *
     * @todo Move it to core/base/common near FlatJaggedArray?
     */
    class Matrix {
    public:
      inline void
        alloc(const size_t nRows, const size_t nCols, const double fill = 0.0) {
        this->nRows_ = nRows;
        this->nCols_ = nCols;
        this->matrix_.resize(this->nRows_ * this->nCols_, fill);
      }
      inline Matrix() = default;
      inline Matrix(const size_t nRows, const size_t nCols) {
        this->alloc(nRows, nCols);
      }
      inline Matrix(const size_t nRows, const size_t nCols, const double fill) {
        this->alloc(nRows, nCols, fill);
      }

      inline double &get(const size_t r, const size_t c) {
        return this->matrix_[r * this->nCols_ + c];
      }
      inline const double &get(const size_t r, const size_t c) const {
        return this->matrix_[r * this->nCols_ + c];
      }
      inline const std::vector<double> &data() const {
        return this->matrix_;
      }
      inline size_t nRows() const {
        return this->nRows_;
      }
      inline size_t nCols() const {
        return this->nCols_;
      }
      inline void fill(const double val) {
        std::fill(this->matrix_.begin(), this->matrix_.end(), val);
      }
      inline void toCSV(const std::string &fName) const {
        std::ofstream out(fName);
        // headers
        for(size_t i = 0; i < this->nRows(); ++i) {
          out << "Column" << std::setfill('0') << std::setw(4) << i;
          if(i == this->nRows() - 1) {
            out << '\n';
          } else {
            out << ',';
          }
        }
        // data
        for(size_t i = 0; i < this->nRows(); ++i) {
          for(size_t j = 0; j < this->nCols(); ++j) {
            out << this->get(i, j);
            if(j == this->nCols() - 1) {
              out << '\n';
            } else {
              out << ',';
            }
          }
        }
      }

    private:
      std::vector<double> matrix_{};
      size_t nRows_{};
      size_t nCols_{};
    };

    inline void
      preconditionTriangulation(AbstractTriangulation *const triangulation) {
      if(triangulation != nullptr) {
        triangulation->preconditionEdges();
        triangulation->preconditionVertexEdges();
      }
    }

    /**
     * @brief Update the projected points when the dilatation
     * coefficient is modified. Avoids to recompute the whole
     * mapper when there is no need for it, hence spees the
     * calculations in this case.
     *
     * @param[out] outputPointsCoords stores the new coordinates:
     *  outputPointsCoords[3*i+k] stores the k-th coordinate for
     *  the i-th point.
     * @param[in] pointsPrev pointsPrev contains the previous
     * coordinates of each point, used to update them.
     *
     * @return 0 in case of success.
     */
    int updateNonCentroidsCoords(float *const outputPointsCoords, const std::vector<float> &pointsPrev);

void updateNonCentroidPointsAlpha(float *outputPointsCoords,
                                               //SimplexId* const centroidId,
                                               SimplexId* const connCompId,
                                               const Matrix &highDimDistMat,
                                               const size_t nbPoint,
                                               const double alpha);

    /**
     * @brief Compute mapper
     *
     * @param[out] outputBucket associate its Bucket id to each vertex
     * @param[out] outputConnComps Component id for each vertex (in
     * the strict sense i.e. vertices belonging to the corresponding
     * bucket)
     * @param[out] compBaryCoords 3D barycenters coordinates
     * @param[out] compArcs Arcs between connected components
     * @param[in] distMat High-dimension input distance matrix (for
     * re-embedding)
     * @param[in] inputSf Input scalar field (on vertices)
     * @param[in] triangulation Triangulation
     *
     * @return 0 in case of success.
     */
    template <typename dataType, typename triangulationType>
    int execute(int *const outputBucket,
                int *const outputConnComp,
                std::vector<std::array<float, 3>> &compBaryCoords,
                std::vector<int> &connCompBucket,
                std::vector<std::set<SimplexId>> &compArcs,
                float *const outputPointsCoords,
                const Matrix &distMat,
                const dataType *const inputSf,
                const triangulationType &triangulation);

  private:
    /**
     * @brief Map vertices & edges to buckets
     *
     * @param[out] vertsBucket Bucket id for each vertex
     * @param[out] bucketEdges List of edges in bucket per bucket
     * @param[in] inputSf Input scalar field (on vertices)
     * @param[in] triangulation Triangulation
     *
     * @return 0 in case of success.
     */
    template <typename dataType, typename triangulationType>
    int findBuckets(int *const vertsBucket,
                    std::vector<std::vector<SimplexId>> &bucketEdges,
                    const dataType *const inputSf,
                    const triangulationType &triangulation) const;

    /**
     * @brief Detect connected components for a given bucket.
     *
     * @param[in] bucketId Id of current bucket
     * @param[out] connComps Component id for each vertex
     * @param[in] vertsBucket Bucket id for each vertex
     * @param[out] connCompEdges List of edges per connected component
     * @param[in] bucketEdges List of edges in the current bucket (reused
     * storage)
     * @param[in] edgeGidToLid Global to compact edge mapping (reused storage)
     * @param[in] uf Local edge Union-Find vector (reused storage)
     * @param[out] connCompIds Local edge connected component id (reused
     * storage)
     * @param[in] triangulation Triangulation
     *
     * @return 0 in case of success.
     */
    template <typename triangulationType>
    int processBucket(const int bucketId,
                      int *const connComps,
                      const int *const vertsBucket,
                      std::vector<std::vector<SimplexId>> &connCompEdges,
                      std::vector<SimplexId> &bucketEdges,
                      std::vector<SimplexId> &edgeGidToLid,
                      std::vector<ttk::UnionFind> &uf,
                      std::vector<int> &connCompIds,
                      const triangulationType &triangulation) const;

    /**
     * @brief Compute connected component barycenter.
     *
     * @param[out] baryCoords 3D coordinates for current component
     * @param[in] connCompEdges List of edges of the current component
     * @param[in] triangulation Triangulation
     *
     * @return 0 in case of success.
     */
    template <typename triangulationType>
    int computeCompBarycenter(std::array<float, 3> &baryCoords,
                              const std::vector<SimplexId> &connCompEdges,
                              const triangulationType &triangulation) const;

    template <typename triangulationType>
    int computeCompHighDimBarycenter(
      double *const highDimBaryCoords,
      const Matrix &inputHighDimCoords,
      const std::vector<SimplexId> &connCompVertices,
      const size_t nHighDims,
      const triangulationType &triangulation) const;

    /**
     * @brief Call DimensionReduction
     *
     * @param[out] outputCoords Embedding coordinates (one array per component)
     * @param[in] mat Distance matrix to be reduced
     * @param[in] isDistanceMatrix If @p mat is a distance matrix
     * @param[in] method DimensionReduction method
     *
     * @return 0 in case of success
     */
    int reduceMatrix(std::vector<std::vector<double>> &outputCoords,
                     const Matrix &mat,
                     const bool isDistanceMatrix,
                     const ttk::DimensionReduction::METHOD method
                     = ttk::DimensionReduction::METHOD::MDS) const;

    /**
     * @brief Compute connected components centroid
     *
     * The centroid is the connected component vertex that minimizes
     * the square distance to every other vertices in the component.
     *
     * @param[out] centroidId Centroid vertex identifier per connected component
     * @param[in] distMat Input distance matrix between all vertices
     * @param[in] connCompVertices Set of vertices per connected
     * component (in the broad sense i.e. vertices of edges that cross
     * the connected component)
     * @param[in] outputConnComp Connected component id per vertex (in
     * the strict sense i.e. vertices belonging to the corresponding
     * bucket)
     */
    void computeCompCentroid(
      std::vector<SimplexId> &centroidId,
      const Matrix &distMat,
      const std::vector<std::set<SimplexId>> &connCompVertices,
      const int *const outputConnComp) const;

    /**
     * @brief Extract sub-matrix from distance matrix
     *
     * @param[out] subDistMat Extracted sub-matrix
     * @param[in] vertsId List of rows/columns identifiers
     * @param[in] distMat Input distance matrix
     */
    void extractSubDistMat(Matrix &subDistMat,
                           const std::vector<SimplexId> &vertsId,
                           const Matrix &distMat) const;

    //TODO doc
    void computeGlobalWeightedDistMatrix(Matrix &globalDistMat,
        const Matrix &centroidDistMat,
        const std::vector<int> &centroidId,
        const int *const outputConnComp,
        const Matrix &highDimDistMat,
        const size_t nbPoint,
        const double alpha) const;
    
    /**
     * @brief Re-embed by reducing centroids & individual connected components
     *
     * @param[out] compBaryCoords Output centroid 3D coordinates
     * @param[out] outputPointsCoords Output dataset vertices 3D coordinates
     * @param[in] distMat Input distance matrix between dataset vertices
     * @param[in] outputConnComp Component id per vertex (computed by execute())
     * @param[in] compArcs Arcs between components (computed by execute())
     * @param[in] connCompEdges Edges per component (computed by execute())
     * @param[in] triangulation Input mesh
     */
    template <typename triangulationType>
    int reEmbedMapper(std::vector<std::array<float, 3>> &compBaryCoords,
                      float *const outputPointsCoords,
                      const Matrix &distMat,
                      const int *const outputConnComp,
                      const std::vector<std::set<SimplexId>> &compArcs,
                      const std::vector<std::vector<SimplexId>> &connCompEdges,
                      const std::vector<int> &connCompBucket,
                      const triangulationType &triangulation);

  protected:
    int NumberOfBuckets{10};
    LOWER_DIMENSION LowerDimension{LOWER_DIMENSION::LOWER_DIM_2D};
    REDUCTION_ALGO ReductionAlgo{REDUCTION_ALGO::MDS};
    REEMBED_METHOD ReembedMethod{REEMBED_METHOD::ARCS_GEODESIC};
    bool ReEmbedMapper{false};
    double DilatationCoeff{0.4};
    double AlphaCoeff{0.0};

    // Variables used to update the coordinates when the dilatation
    // coefficient is changed.
    double prevDilatationCoeff_{-1};
    std::vector<double> compSpecialCoeffToSave_{};
    Matrix prevCentroidDistMat_;
    std::vector<int> prevCentroidId_;
  };

} // namespace ttk

int ttk::Mapper::updateNonCentroidsCoords(float *const outputPointsCoords, const std::vector<float> &pointsPrev)
{
  const size_t nbPoint = compSpecialCoeffToSave_.size()/3;
  if (pointsPrev.size() != 3*nbPoint)
  {
    printErr("Error in updating the dilatation coefficient only. We want to update the coordinates of "
        + std::to_string(nbPoint) + " points but the data we saved to do so concerns "
        + std::to_string(compSpecialCoeffToSave_.size()) + " points.");
    return 1;
  }

  size_t dim = LowerDimension == LOWER_DIMENSION::LOWER_DIM_2D ? 2:3;

  for (size_t i = 0; i < nbPoint; i++)
  {
    if (i < 5)
    {
      std::cerr << i << " pPrev => " << pointsPrev[3*i] << "," << pointsPrev[3*i+1] << "," << pointsPrev[3*i+2] << std::endl;
      std::cerr << " oh la qui voilà = " << compSpecialCoeffToSave_[3*i+0] << " - " << compSpecialCoeffToSave_[3*i+1] << std::endl;
   }
    for (size_t iDim = 0; iDim < dim; iDim++)
    {
      outputPointsCoords[3*i+iDim] = pointsPrev[3*i+iDim] + (DilatationCoeff - prevDilatationCoeff_) * compSpecialCoeffToSave_[3*i+iDim];
    }
    if (i < 5)
      std::cerr << i << " => " << outputPointsCoords[3*i] << "," << outputPointsCoords[3*i+1] << "," << outputPointsCoords[3*i+2] << std::endl;
  }
  prevDilatationCoeff_ = DilatationCoeff;


  return 0;
}

void ttk::Mapper::updateNonCentroidPointsAlpha(float *outputPointsCoords,
                                               //SimplexId* const centroidId,
                                               SimplexId* const connCompId,
                                               const Matrix &highDimDistMat,
                                               const size_t nbPoint,
                                               const double alpha)
{
  printErr("Updating alpha :-)");
  Matrix globalDistMat;
  computeGlobalWeightedDistMatrix(globalDistMat, prevCentroidDistMat_, prevCentroidId_, connCompId, highDimDistMat, nbPoint, alpha);
  /*
  globalDistMat.alloc(nbPoint, nbPoint, 0);
  for (size_t i1 = 0; i1 < nbPoint; i1++)
  {
    for (size_t i2 = 0; i2 < i1; i2++)
    {
      double squareCentroidDist = 0;
      size_t comp1 = connCompId[i1], comp2 = connCompId[i2];
      size_t centroid1 = prevCentroidId_[comp1], centroid2 = prevCentroidId_[comp2];
      globalDistMat.get(i1, i2) = (1-alpha)*(highDimDistMat.get(i1,centroid1)+prevCentroidDistMat_.get(comp1, comp2)+highDimDistMat.get(i2, centroid2))+alpha*highDimDistMat.get(i1, i2);
      globalDistMat.get(i2, i1) = (1-alpha)*(highDimDistMat.get(i1,centroid1)+prevCentroidDistMat_.get(comp1, comp2)+highDimDistMat.get(i2, centroid2))+alpha*highDimDistMat.get(i1, i2);
    }
  }
  */
    size_t dim = LowerDimension == LOWER_DIMENSION::LOWER_DIM_2D ? 2:3;
    std::vector<std::vector<double>> coordsAll;
    reduceMatrix(coordsAll, globalDistMat, true, this->ReductionAlgo);
    for (size_t iPt = 0; iPt < globalDistMat.nRows(); iPt++)
    {
      for (size_t iDim = 0; iDim < dim; iDim++)
      {
        outputPointsCoords[3*iPt+iDim] = coordsAll[iDim][iPt];
      }
    }

  /*
    for (size_t iCentr = 0; iCentr < compBaryCoords.size(); iCentr++)
    {
      size_t iVert = centroidId[iCentr];
      for (size_t iDim = 0; iDim < 3; iDim++)
      {
        compBaryCoords[iVert][iDim] = outputConnComp[3*iVert+iDim];
      }
    }*/

}



// template functions
template <typename dataType, typename triangulationType>
int ttk::Mapper::execute(int *const outputBucket,
                         int *const outputConnComp,
                         std::vector<std::array<float, 3>> &compBaryCoords,
                         std::vector<int> &connCompBucket,
                         std::vector<std::set<SimplexId>> &compArcs,
                         float *const outputPointsCoords,
                         const Matrix &distMat,
                         const dataType *const inputSf,
                         const triangulationType &triangulation) {

  Timer tm{}, tmsec{};

  std::vector<std::vector<SimplexId>> bucketEdges(this->NumberOfBuckets);
  // We fill outputBucket: associates to each vertex its bucket id,
  // and bucketEdes array: associates to each bucket the list of edges crossing
  // this bucket.
  this->findBuckets(outputBucket, bucketEdges, inputSf, triangulation);

  // allocate memory
  compBaryCoords.reserve(2 * this->NumberOfBuckets);
  compArcs.reserve(2 * this->NumberOfBuckets);

  std::vector<SimplexId> edgeGidToLid(triangulation.getNumberOfEdges(), -1);
  std::vector<ttk::UnionFind> uf{};
  std::vector<int> connCompIds{};
  std::vector<std::vector<std::vector<SimplexId>>> connCompsEdgesPerBucket(
    this->NumberOfBuckets);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_) \
  firstprivate(edgeGidToLid, uf, connCompIds)
#endif // TTK_ENABLE_OPENMP
  // We compute the connected components inside buckets.
  for(int i = 0; i < this->NumberOfBuckets; ++i) {
    this->processBucket(i, outputConnComp, outputBucket,
                        connCompsEdgesPerBucket[i], bucketEdges[i],
                        edgeGidToLid, uf, connCompIds, triangulation);
  }

  this->printMsg("Detected connected components", 1.0, tmsec.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);
  tmsec.reStart();

  size_t nConnComps{};
  // (bucket id, local component id) -> global component id
  std::map<std::pair<int, size_t>, int> connCompLidToGid{};
  for(int i = 0; i < this->NumberOfBuckets; ++i) {
    for(size_t j = 0; j < connCompsEdgesPerBucket[i].size(); ++j) {
      connCompLidToGid[{i, j}] = nConnComps++;
    }
  }

  // offset component id
  const auto nVerts{triangulation.getNumberOfVertices()};
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nVerts; ++i) {
    outputConnComp[i] = connCompLidToGid[{outputBucket[i], outputConnComp[i]}];
  }

  this->printMsg("Generated segmentation", 1.0, tmsec.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);
  tmsec.reStart();

  std::vector<std::vector<SimplexId>> connCompEdges(nConnComps);
  // aggregate component edges
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < this->NumberOfBuckets; ++i) {
    for(size_t j = 0; j < connCompsEdgesPerBucket[i].size(); ++j) {
      connCompEdges[connCompLidToGid[{i, j}]]
        = std::move(connCompsEdgesPerBucket[i][j]);
    }
  }

  // find barycenter coordinates
  compBaryCoords.resize(nConnComps);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  // We compute the barycentre of each connected component. This will be used to
  // reduce the component to its barycentre: we take the barycentre of the
  // middle of each edge in this component.
  for(size_t i = 0; i < nConnComps; ++i) {
    this->computeCompBarycenter(
      compBaryCoords[i], connCompEdges[i], triangulation);
  }

  this->printMsg("Found barycenter coordinates", 1.0, tmsec.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);
  tmsec.reStart();

  // transposition of connCompEdges
  std::vector<std::vector<int>> edgeConnComps(triangulation.getNumberOfEdges());
  for(size_t i = 0; i < connCompEdges.size(); ++i) {
    for(const auto e : connCompEdges[i]) {
      edgeConnComps[e].emplace_back(i);
    }
  }

  // global component id -> bucket id
  connCompBucket.resize(nConnComps);
  for(const auto &p : connCompLidToGid) {
    connCompBucket[p.second] = p.first.first;
  }

  // We detect the arcs between connected components:
  // for each edge in our base graph, if it belongs to two components of
  // contiguous buckets, then we put an edge between the barycenters that
  // represents them.
  compArcs.resize(nConnComps);

  for(size_t i = 0; i < edgeConnComps.size(); ++i) {
    // an arc edge should at least cover two connected components
    if(edgeConnComps[i].size() < 2) {
      continue;
    }
    for(size_t j = 1; j < edgeConnComps[i].size(); ++j) {
      const auto currCompId{edgeConnComps[i][j]};
      for(size_t k = 0; k < j; ++k) {
        const auto prevCompId{edgeConnComps[i][k]};
        // We compute an adjacency list, putting only the backlink (idx ->idx-1)
        if(connCompBucket[prevCompId] == connCompBucket[currCompId] - 1) {
          compArcs[currCompId].emplace(prevCompId);
        }
      }
    }
  }

  this->printMsg("Detected arcs between connected components", 1.0,
                 tmsec.getElapsedTime(), this->threadNumber_,
                 debug::LineMode::NEW, debug::Priority::DETAIL);

  this->printMsg(
    "Found " + std::to_string(nConnComps) + " connected components", 1,
    tm.getElapsedTime(), this->threadNumber_);

  if(this->ReEmbedMapper) {
    // Resetting the coordinates so that if we passe from 3D to 2D, the third
    // coordinate will not remain there.
    for(auto &t : compBaryCoords) {
      for(float &x : t)
        x = 0;
    }
    this->reEmbedMapper(compBaryCoords, outputPointsCoords, distMat,
                        outputConnComp, compArcs, connCompEdges, connCompBucket,
                        triangulation);
  }

    return 0;
}

template <typename dataType, typename triangulationType>
int ttk::Mapper::findBuckets(int *const vertsBucket,
                             std::vector<std::vector<SimplexId>> &bucketEdges,
                             const dataType *const inputSf,
                             const triangulationType &triangulation) const {

  Timer tm{};

  if(this->NumberOfBuckets < 0) {
    this->printErr("Number of buckets should not be negative");
    return -1;
  }

  const auto nVerts{triangulation.getNumberOfVertices()};
  const auto sfRange = std::minmax_element(inputSf, inputSf + nVerts);

  std::vector<double> buckets(this->NumberOfBuckets + 1);
  buckets.front() = *sfRange.first;
  buckets.back() = *sfRange.second;
  const double bucketLength{(buckets.back() - buckets.front())
                            / this->NumberOfBuckets};
  for(int i = 1; i < this->NumberOfBuckets; ++i) {
    buckets[i] = buckets[i - 1] + bucketLength;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < nVerts; ++i) {
    for(int j = 0; j < this->NumberOfBuckets; ++j) {
      if(inputSf[i] < buckets[j + 1]) {
        vertsBucket[i] = j;
        break;
      }
      vertsBucket[i] = this->NumberOfBuckets - 1;
    }
  }

  const auto nEdges{triangulation.getNumberOfEdges()};

  // find edges covering bucket
  for(SimplexId e = 0; e < nEdges; ++e) {
    SimplexId v0{}, v1{};
    triangulation.getEdgeVertex(e, 0, v0);
    triangulation.getEdgeVertex(e, 1, v1);
    auto b0 = vertsBucket[v0];
    auto b1 = vertsBucket[v1];
    if(b0 > b1) {
      std::swap(b0, b1);
    }
    for(int i = b0; i < b1 + 1; ++i) {
      bucketEdges[i].emplace_back(e);
    }
  }

  this->printMsg("Placed vertices & edges in buckets", 1.0, tm.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::DETAIL);

  return 0;
}

// On associe aussi aux sommets (et arêtes) les id des composantes connexes
// auxquelles ils apparaissent. Et attention, deux types d'appartenance pour les
// sommets : appartient à un bucket si une arête incidente le croise (v1), ou
// bien s'il est contenu lui-même dans le bucket (v2)
template <typename triangulationType>
int ttk::Mapper::processBucket(
  const int bucketId,
  int *const connComps,
  const int *const vertsBucket,
  std::vector<std::vector<SimplexId>> &connCompEdges,
  std::vector<SimplexId> &bucketEdges,
  std::vector<SimplexId> &edgeGidToLid,
  std::vector<ttk::UnionFind> &uf,
  std::vector<int> &connCompIds,
  const triangulationType &triangulation) const {

  // init uf, connCompIds vectors
  uf.resize(bucketEdges.size());
  connCompIds.resize(bucketEdges.size());
  for(size_t i = 0; i < bucketEdges.size(); ++i) {
    edgeGidToLid[bucketEdges[i]] = i;
    uf[i].setParent(&uf[i]);
    uf[i].setRank(i);
    connCompIds[i] = -1;
  }

  // find connected components in edges
  for(size_t i = 0; i < bucketEdges.size(); ++i) {
    const auto e{bucketEdges[i]};
    for(SimplexId j = 0; j < 2; ++j) {
      SimplexId v{};
      triangulation.getEdgeVertex(e, j, v);
      for(SimplexId k = 0; k < triangulation.getVertexEdgeNumber(v); ++k) {
        SimplexId oe{};
        triangulation.getVertexEdge(v, k, oe);
        const auto l{edgeGidToLid[oe]};
        if(oe >= e || l == -1) {
          continue;
        }
        ttk::UnionFind::makeUnion(&uf[i], &uf[l]);
      }
    }
  }

  size_t nLocConnComps{};
  // find connected component ids
  for(size_t i = 0; i < uf.size(); ++i) {
    const auto root{uf[i].find()};
    if(root != nullptr) {
      const auto rank{root->getRank()};
      if(connCompIds[rank] == -1) {
        connCompIds[rank] = nLocConnComps++;
      }
      connCompIds[i] = connCompIds[rank];
    }
  }

  // store component edges ids (useful to get arcs)
  connCompEdges.resize(nLocConnComps);
  for(auto &vec : connCompEdges) {
    // random guess...
    vec.reserve(bucketEdges.size() / nLocConnComps + 1);
  }

  for(size_t i = 0; i < bucketEdges.size(); ++i) {
    const auto e{bucketEdges[i]};
    const auto compId{connCompIds[i]};
    if(compId == -1) {
      continue;
    }
    connCompEdges[compId].emplace_back(e);
    // also put component id on (relevant) edge vertices
    for(SimplexId j = 0; j < 2; ++j) {
      SimplexId v{};
      triangulation.getEdgeVertex(e, j, v);
      if(vertsBucket[v] == bucketId) {
        connComps[v] = compId;
      }
    }
  }

  // cleanup
  for(const auto e : bucketEdges) {
    edgeGidToLid[e] = -1;
  }
  bucketEdges.clear();

  return 0;
}

template <typename triangulationType>
int ttk::Mapper::computeCompBarycenter(
  std::array<float, 3> &baryCoords,
  const std::vector<SimplexId> &connCompEdges,
  const triangulationType &triangulation) const {

  size_t nVertsInComp{};

  for(const auto e : connCompEdges) {
    for(SimplexId j = 0; j < 2; ++j) {
      SimplexId v{};
      triangulation.getEdgeVertex(e, j, v);
      // get barycenter from all vertices from all considered edges
      std::array<float, 3> pt{};
      triangulation.getVertexPoint(v, pt[0], pt[1], pt[2]);
      baryCoords[0] += pt[0];
      baryCoords[1] += pt[1];
      baryCoords[2] += pt[2];
      nVertsInComp++;
    }
  }

  // store barycenter coordinates
  baryCoords[0] /= nVertsInComp;
  baryCoords[1] /= nVertsInComp;
  baryCoords[2] /= nVertsInComp;

  return 0;
}

// Never called so far.
template <typename triangulationType>
int ttk::Mapper::computeCompHighDimBarycenter(
  double *const highDimBaryCoords,
  const Matrix &inputHighDimCoords,
  const std::vector<SimplexId> &connCompVertices,
  const size_t nHighDims,
  const triangulationType &triangulation) const {

  for(const auto v : connCompVertices) {
    std::array<float, 3> pt{};
    triangulation.getVertexPoint(v, pt[0], pt[1], pt[2]);

    for(size_t i = 0; i < nHighDims; ++i) {
      highDimBaryCoords[i] += inputHighDimCoords.get(v, i);
    }
  }

  for(size_t i = 0; i < nHighDims; ++i) {
    highDimBaryCoords[i] /= connCompVertices.size();
  }

  return 0;
}

template <typename triangulationType>
int ttk::Mapper::reEmbedMapper(
  std::vector<std::array<float, 3>> &compBaryCoords,
  float *const outputPointsCoords,
  const Matrix &distMat,
  const int *const outputConnComp,
  const std::vector<std::set<SimplexId>> &compArcs,
  const std::vector<std::vector<SimplexId>> &connCompEdges,
  const std::vector<int> &connCompBucket,
  const triangulationType &triangulation) {

  Timer tm{};
  size_t dim = LowerDimension == LOWER_DIMENSION::LOWER_DIM_2D ? 2:3;

  compSpecialCoeffToSave_.clear();
  compSpecialCoeffToSave_.resize(3*triangulation.getNumberOfVertices(), 0.0);
  std::cerr << "Resized the compSpecial to " + std::to_string(compSpecialCoeffToSave_.size());

  // 1. extract vertices in component edges. A vertex is considered in a bucket
  // if it lies insied or if it is the extremity of an edge which crosses that
  // bucked.
  std::vector<std::set<SimplexId>> connCompVertices(connCompEdges.size());
  for(size_t i = 0; i < connCompEdges.size(); ++i) {
    const auto &cce{connCompEdges[i]};
    auto &ccv{connCompVertices[i]};
    for(size_t j = 0; j < cce.size(); ++j) {
      SimplexId v0{}, v1{};
      triangulation.getEdgeVertex(cce[j], 0, v0);
      triangulation.getEdgeVertex(cce[j], 1, v1);
      ccv.emplace(v0);
      ccv.emplace(v1);
    }
  }

  // 2. sort vertices per connected component. Here a vertex belongs to a
  // connected component at a stricter sense: it must belong to the component
  // bucket.
  std::vector<std::vector<SimplexId>> connCompVertsStrict(
    connCompVertices.size());
  for(SimplexId i = 0; i < triangulation.getNumberOfVertices(); ++i) {
    connCompVertsStrict[outputConnComp[i]].emplace_back(i);
  }

  // 3. compute connected components centroids from input distance matrix: one
  // centroid per component. The centroid is the existing point with minimum
  // average distance to the other points. Difference with a barycentre: the
  // barycentre might not belong to the input points.
  std::vector<SimplexId> centroidId{};
  this->computeCompCentroid(
    centroidId, distMat, connCompVertices, outputConnComp);

  // 4+5: try to obtain an embedding (low-dimension coordinates) from the list
  // of cendroids and the distance matrix.

  // 4. extract distance matrix between centroids. TODO choices
  const double DOUBLE_MAX = std::numeric_limits<double>::max();
  Matrix centroidsDistMat(compArcs.size(), compArcs.size(), DOUBLE_MAX);
  size_t nComp = compArcs.size();

  std::vector<std::pair<std::pair<int, int>, double>> initAdjListWithCosts;
  initAdjListWithCosts.reserve(
    2 * nComp); // Guess because linked to about 2 buckets (prev and next ones)

  /*
  double distMax = 0;
  for (size_t i1 = 0; i1 < nComp; i1++)
  {
    for (size_t i2 = i1+1; i2 < nComp; i2++)
    {
      double cur = distMat.get(i1,i2);
      distMax = max(distMax, cur);
    }
  }*/

  if(ReembedMethod == REEMBED_METHOD::ARCS_GEODESIC) {
    for(size_t i = 0; i < nComp; ++i) {
      for(const auto el : compArcs[i]) {
        double curDist = distMat.get(centroidId[i], centroidId[el]);
        centroidsDistMat.get(i, el) = curDist;
        centroidsDistMat.get(el, i) = curDist;
        initAdjListWithCosts.push_back({{i, el}, curDist});
      }
    }
  }

  else if(ReembedMethod == REEMBED_METHOD::BUCKETS_GEODESIC) {
    for(size_t i = 0; i < nComp; ++i) {
      for(size_t el = i + 1; el < nComp; ++el) {
        if(abs(connCompBucket[i] - connCompBucket[el]) == 1) {
          double curDist = distMat.get(centroidId[i], centroidId[el]);
          centroidsDistMat.get(i, el) = curDist;
          centroidsDistMat.get(el, i) = curDist;
          initAdjListWithCosts.push_back({{i, el}, curDist});
        }
      }
    }
  } else {
    this->printErr(" The reembed method must be chosen between arcs or buckets "
                   "geodesic. The current value is "
                   + std::to_string((int)ReembedMethod) + ".");
    return 0;
  }

  this->printMsg(
    "Computed initial arcs", 1, tm.getElapsedTime(), this->threadNumber_);
  tm.reStart();

  // TODO faster algorithm for sparse graph
  for(size_t u = 0; u < nComp; u++)
    centroidsDistMat.get(u, u) = 0;

  for(size_t inter = 0; inter < nComp; inter++) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t u = 0; u < nComp; u++) {
      for(size_t v = 0; v < nComp; v++) {
        if(centroidsDistMat.get(u, inter) + centroidsDistMat.get(inter, v)
           < centroidsDistMat.get(u, v))
          centroidsDistMat.get(u, v)
            = centroidsDistMat.get(u, inter) + centroidsDistMat.get(inter, v);
      }
    }
  }

  /*
  for(size_t i = 0; i < compArcs.size(); ++i) {
    for(const auto el : compArcs[i]) {
    //for(size_t el = 0; el < compArcs.size(); ++el) {
      if (abs(connCompBucket[i] - connCompBucket[el]) <= 1)
      {
        centroidsDistMat.get(i, el) = distMat.get(centroidId[i],
centroidId[el]); centroidsDistMat.get(el, i) = distMat.get(centroidId[i],
centroidId[el]);
      }
    }
  }
}
*/
  prevCentroidDistMat_ = centroidsDistMat;
  prevCentroidId_ = centroidId;
  this->printMsg(
    "Computed geodesics", 1.0, tm.getElapsedTime(), this->threadNumber_);
  tm.reStart();

  // 5. get an embedding of the distance matrix between the centroids
  std::vector<std::vector<double>> embedCentroids{};
  this->reduceMatrix(
    embedCentroids, centroidsDistMat, true, this->ReductionAlgo);
  this->printErr(std::to_string(embedCentroids.size()) + " : dimension\n");

  // Conversion between storage layout x0,y0,z0,x1,y1,z1... and x0,
  // x1...,y0,y1...,z0,z1...
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < embedCentroids.size(); ++i) {
    for(size_t j = 0; j < embedCentroids[i].size(); ++j) {
      compBaryCoords[j][i] = embedCentroids[i][j];
    }
  }

  // matrice à partir des coordonnées calculées
  // TODO: voir si on peut extraire à partir matrice de l'entrée
  // 6. compute a distance matrix from the embedded centroids
  std::vector<const float *> inputs(compBaryCoords.size());
  for(size_t i = 0; i < compBaryCoords.size(); ++i) {
    inputs[i] = compBaryCoords[i].data();
  }
  std::vector<std::vector<double>> centroidsEmbedDistMat{};
  ttk::LDistanceMatrix ldm{};
  if(this->LowerDimension == LOWER_DIMENSION::LOWER_DIM_2D) {
    ldm.execute(centroidsEmbedDistMat, inputs, 2);
  } else {
    ldm.execute(centroidsEmbedDistMat, inputs, 3);
  }

  this->printMsg(". Re-embedded centroids", 1.0, tm.getElapsedTime(),
                 this->threadNumber_, debug::LineMode::NEW,
                 debug::Priority::PERFORMANCE);
  //return 0; //To avoid computing segmentation for faster tests.
  


  if (AlphaCoeff > 1e-6)
  {
    std::cout << "projecting GLOBAL !" << std::endl;
    Matrix globalWeightedMatrix;
    this->computeGlobalWeightedDistMatrix(globalWeightedMatrix,
                                    centroidsDistMat,
                                    centroidId,
                                    outputConnComp,
                                    distMat,
                                    distMat.nRows(),
                                    0.0);

    std::vector<std::vector<double>> coordsAll;
    reduceMatrix(coordsAll, globalWeightedMatrix, true, this->ReductionAlgo);
    
    for (size_t iPt = 0; iPt < distMat.nRows(); iPt++)
    {
      size_t idComp = outputConnComp[iPt];
        size_t idCentroid = centroidId[idComp];
      
      for (size_t iDim = 0; iDim < dim; iDim++)
      {
        double delta = coordsAll[iDim][iPt]-coordsAll[iDim][idCentroid];
        if (!std::isfinite(coordsAll[iDim][iPt]))
        {
            std::cerr << "NO, INFINITE\n";
            delta = 0;
        }
        //std::cout << delta << " + " << embedCentroids[iDim][idCentroid] << "    --->  " << idCentroid << "," << outputConnComp[iPt] << std::endl;
        float newCoord = delta;// + embedCentroids[iDim][idComp];
        //outputPointsCoords[3*iPt+iDim] = coordsAll[iDim][iPt];
        //std::cout << outputPointsCoords[3*iPt+iDim] << "\n";
        outputPointsCoords[3*iPt+iDim] = coordsAll[iDim][idCentroid] + newCoord;
      }
    }
     //Pour faire coller les barycentres à la projection globale nouvellement calculée
    for (size_t iComp = 0; iComp < compBaryCoords.size(); iComp++)
    {
      size_t iVert = centroidId[iComp];
      //printErr("ICentr = " + std::to_string(iComp) + " and iVert = " + std::to_string(iVert));
      for (size_t iDim = 0; iDim < dim; iDim++)
      {
        //compBaryCoords.at(iComp).at(iDim) = coordsAll[iDim][iVert];//outputConnComp[3*iVert+iDim];
        compBaryCoords.at(iComp).at(iDim) = outputConnComp[3*iVert+iDim];
      }
    }
    printErr("NOT RETURNING YEAH\n");
    //return 0;
  }
  else
{
  std::cout << " COMPUTING PROJECTION NORMAL " << std::endl;
  // 7. get an embedding for all vertices of each connected component
  // Not in parallel because it calls some Python code. Parallelising the calls
  // to Python causes errors.
  for(size_t i = 0; i < connCompEdges.size(); ++i) {

    Timer tmcomp{};

    if(connCompVertsStrict[i].empty()) {
      continue;
    }

    if(connCompVertsStrict[i].size() == 1) {
      // lock point to centroid
      for(size_t k = 0; k < 3; ++k) {
        outputPointsCoords[3 * connCompVertsStrict[i][0] + k]
          = compBaryCoords[i][k];
      }
      for (size_t j = 0; j < dim; j++)
        compSpecialCoeffToSave_[3*connCompVertsStrict[i][0]+j] = 0;
      continue;
    }
    //TODO dilatation en fait coefficient adaptable si lock sur centroide...

    // return 0;

    Matrix distMatConnComp{};
    std::vector<std::vector<double>> embedConnComp{};
    this->extractSubDistMat(distMatConnComp, connCompVertsStrict[i], distMat);
    this->reduceMatrix(
        embedConnComp, distMatConnComp, true, this->ReductionAlgo);
    // get max distance between points in sub-distance matrix
    double compDiag{};
    for(size_t j = 0; j < distMatConnComp.nCols() - 1; ++j) {
      for(size_t k = j + 1; k < distMatConnComp.nCols(); ++k) {
        if(distMatConnComp.get(j, k) > compDiag) {
          compDiag = distMatConnComp.get(j, k);
        }
      }
    }

    // get min distance to neighbor in centroidsEmbedDistMat
    auto maxDistNeigh = std::numeric_limits<double>::max();
    for(size_t j = 0; j < centroidsEmbedDistMat[i].size(); ++j) {
      if(j == i) {
        continue;
      }
      if(centroidsEmbedDistMat[i][j] < maxDistNeigh) {
        maxDistNeigh = centroidsEmbedDistMat[i][j];
      }
    }

    // resize & shift component reduced coordinates
    // after this operation, it is centered on 0,0
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t j = 0; j < embedCentroids.size(); ++j) {
      for (size_t iPtComp = 0; iPtComp < embedConnComp[j].size(); iPtComp++)
      {
        auto &coords = //outputPointsCoords[3*connCompVertsStrict[i][iPtComp]+j]
                       embedConnComp[j][iPtComp];
        compSpecialCoeffToSave_[3*connCompVertsStrict[i][iPtComp]+j] = coords*maxDistNeigh/compDiag;
        std::cout << "\t\tlà worth " << compSpecialCoeffToSave_[3*connCompVertsStrict[i][iPtComp]+j] << std::endl;
        coords *= DilatationCoeff * maxDistNeigh / compDiag;
        coords += embedCentroids[j][i];
      }
    }

    // store reduced components in output vector
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(size_t j = 0; j < connCompVertsStrict[i].size(); ++j) {
      for(size_t k = 0; k < embedConnComp.size(); ++k) {

        if(std::isfinite(embedConnComp[k][j])) { // to avoid nan problems
          outputPointsCoords[3 * connCompVertsStrict[i][j] + k]
          = embedConnComp[k][j];

        // connCompVertsStrict[i][j] : global id of the jth point in the ith
        // connected component
        } else { // in case the value is infinite, we take the centroi coordinate
          outputPointsCoords[3 * connCompVertsStrict[i][j] + k]
            = compBaryCoords[i][k]; // = embedCentroids[k][i]
        }
      }
    }

    this->printMsg(".. Re-embedded component " + std::to_string(i), 1.0,
        tmcomp.getElapsedTime(), this->threadNumber_,
        debug::LineMode::NEW, debug::Priority::DETAIL);
  }
}

  std::cout << "lol" << std::endl;
  this->printMsg(
    "Re-embedded mapper", 1.0, tm.getElapsedTime(), this->threadNumber_);

  prevDilatationCoeff_ = DilatationCoeff;

  return 0;
}
