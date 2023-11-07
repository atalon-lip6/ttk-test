#include <PeriodicImplicitTriangulation.h>
using namespace std;
using namespace ttk;

template <size_t card>
ttk::PeriodicImplicitTriangulation<card>::PeriodicImplicitTriangulation()
  : dimensionality_{-1}, cellNumber_{}, vertexNumber_{}, edgeNumber_{},
    triangleNumber_{}, tetrahedronNumber_{}, isAccelerated_{} {
  this->setDebugMsgPrefix("PeriodicImplicitTriangulation");
  this->hasPeriodicBoundaries_ = true;
}

template<size_t card>
ttk::PeriodicImplicitTriangulation<card>::~PeriodicImplicitTriangulation<card>() = default;

template<size_t card>
int ttk::PeriodicImplicitTriangulation<card>::setInputGrid(const float &xOrigin,
                                                const float &yOrigin,
                                                const float &zOrigin,
                                                const float &xSpacing,
                                                const float &ySpacing,
                                                const float &zSpacing,
                                                const ttk::SimplexId &xDim,
                                                const ttk::SimplexId &yDim,
                                                const ttk::SimplexId &zDim) {

  // Dimensionality //
  if(xDim < 1 or yDim < 1 or zDim < 1)
    dimensionality_ = -1;
  else if(xDim > 1 and yDim > 1 and zDim > 1)
    dimensionality_ = 3;
  else if((xDim > 1 and yDim > 1) or (yDim > 1 and zDim > 1)
          or (xDim > 1 and zDim > 1))
    dimensionality_ = 2;
  else if(xDim > 1 or yDim > 1 or zDim > 1)
    dimensionality_ = 1;
  else
    dimensionality_ = 0;

  // Essentials //
  this->origin_[0] = xOrigin;
  this->origin_[1] = yOrigin;
  this->origin_[2] = zOrigin;
  this->spacing_[0] = xSpacing;
  this->spacing_[1] = ySpacing;
  this->spacing_[2] = zSpacing;
  this->dimensions_[0] = xDim;
  this->dimensions_[1] = yDim;
  this->dimensions_[2] = zDim;
  nbvoxels_[0] = xDim - 1;
  nbvoxels_[1] = yDim - 1;
  nbvoxels_[2] = zDim - 1;

  if(dimensionality_ == 3) {
    // VertexShift
    this->vshift_[0] = xDim;
    this->vshift_[1] = xDim * yDim;
    // EdgeSetDimensions
    esetdims_[0] = xDim * yDim * zDim;
    esetdims_[1] = xDim * yDim * zDim;
    esetdims_[2] = xDim * yDim * zDim;
    esetdims_[3] = xDim * yDim * zDim;
    esetdims_[4] = xDim * yDim * zDim;
    esetdims_[5] = xDim * yDim * zDim;
    esetdims_[6] = xDim * yDim * zDim;
    // EdgeSetShift
    this->esetshift_[0] = esetdims_[0];
    for(int k = 1; k < 7; ++k)
      this->esetshift_[k] = this->esetshift_[k - 1] + esetdims_[k];
    // EdgeShift
    this->eshift_[0] = xDim;
    this->eshift_[1] = xDim * yDim;
    this->eshift_[2] = xDim;
    this->eshift_[3] = xDim * yDim;
    this->eshift_[4] = xDim;
    this->eshift_[5] = xDim * yDim;
    this->eshift_[6] = xDim;
    this->eshift_[7] = xDim * yDim;
    this->eshift_[8] = xDim;
    this->eshift_[9] = xDim * yDim;
    this->eshift_[10] = xDim;
    this->eshift_[11] = xDim * yDim;
    this->eshift_[12] = xDim;
    this->eshift_[13] = xDim * yDim;
    // TriangleSetDimensions
    tsetdims_[0] = xDim * yDim * zDim * 2;
    tsetdims_[1] = xDim * yDim * zDim * 2;
    tsetdims_[2] = xDim * yDim * zDim * 2;
    tsetdims_[3] = xDim * yDim * zDim * 2;
    tsetdims_[4] = xDim * yDim * zDim * 2;
    tsetdims_[5] = xDim * yDim * zDim * 2;
    // TriangleSetShift
    tsetshift_[0] = tsetdims_[0];
    for(int k = 1; k < 6; ++k)
      tsetshift_[k] = tsetshift_[k - 1] + tsetdims_[k];
    // TriangleShift
    tshift_[0] = xDim * 2;
    tshift_[1] = xDim * yDim * 2;
    tshift_[2] = xDim * 2;
    tshift_[3] = xDim * yDim * 2;
    tshift_[4] = xDim * 2;
    tshift_[5] = xDim * yDim * 2;
    tshift_[6] = xDim * 2;
    tshift_[7] = xDim * yDim * 2;
    tshift_[8] = xDim * 2;
    tshift_[9] = xDim * yDim * 2;
    tshift_[10] = xDim * 2;
    tshift_[11] = xDim * yDim * 2;
    // TetrahedronShift
    tetshift_[0] = xDim * 6;
    tetshift_[1] = xDim * yDim * 6;

    this->wrap_[0] = xDim;
    this->wrap_[1] = xDim * yDim;
    this->wrap_[2] = xDim * yDim * zDim;

    // Numbers
    vertexNumber_ = xDim * yDim * zDim;
    this->edgeNumber_ = 0;
    for(int k = 0; k < 7; ++k)
      this->edgeNumber_ += esetdims_[k];
    triangleNumber_ = 0;
    for(int k = 0; k < 6; ++k)
      triangleNumber_ += tsetdims_[k];
    tetrahedronNumber_ = xDim * yDim * zDim * 6;
    cellNumber_ = tetrahedronNumber_;

    checkAcceleration();
  } else if(dimensionality_ == 2) {
    // dimensions selectors
    if(xDim == 1) {
      Di_ = 1;
      this->Dj_ = 2;
    } else if(yDim == 1) {
      Di_ = 0;
      this->Dj_ = 2;
    } else {
      Di_ = 0;
      this->Dj_ = 1;
    }
    // VertexShift
    this->vshift_[0] = this->dimensions_[Di_];
    // EdgeSetDimensions
    esetdims_[0] = this->dimensions_[Di_] * this->dimensions_[this->Dj_];
    esetdims_[1] = this->dimensions_[Di_] * this->dimensions_[this->Dj_];
    esetdims_[2] = this->dimensions_[Di_] * this->dimensions_[this->Dj_];
    // EdgeSetShift
    this->esetshift_[0] = esetdims_[0];
    for(int k = 1; k < 3; ++k)
      this->esetshift_[k] = this->esetshift_[k - 1] + esetdims_[k];
    // EdgeShift
    this->eshift_[0] = this->dimensions_[Di_];
    this->eshift_[2] = this->dimensions_[Di_];
    this->eshift_[4] = this->dimensions_[Di_];
    // TriangleShift
    tshift_[0] = this->dimensions_[Di_] * 2;

    this->wrap_[0] = this->dimensions_[Di_];
    this->wrap_[1] = this->dimensions_[Di_] * this->dimensions_[this->Dj_];

    // Numbers
    vertexNumber_ = this->dimensions_[Di_] * this->dimensions_[this->Dj_];
    this->edgeNumber_ = 0;
    for(int k = 0; k < 3; ++k)
      this->edgeNumber_ += esetdims_[k];
    triangleNumber_ = this->dimensions_[Di_] * this->dimensions_[this->Dj_] * 2;
    cellNumber_ = triangleNumber_;

    checkAcceleration();
  } else if(dimensionality_ == 1) {
    // dimensions selectors
    for(int k = 0; k < 3; ++k) {
      if(this->dimensions_[k] > 1) {
        Di_ = k;
        break;
      }
    }

    // Numbers
    vertexNumber_ = this->dimensions_[Di_];
    this->edgeNumber_ = vertexNumber_;
    cellNumber_ = this->edgeNumber_;
  }

  return 0;
}

template<size_t card>
int ttk::PeriodicImplicitTriangulation<card>::checkAcceleration() {
  isAccelerated_ = false;

  unsigned long long int msb[3];
  if(dimensionality_ == 3) {
    bool allDimensionsArePowerOfTwo = true;
    for(int k = 0; k < 3; ++k)
      if(!isPowerOfTwo(this->dimensions_[k], msb[k]))
        allDimensionsArePowerOfTwo = false;

    if(allDimensionsArePowerOfTwo) {
      mod_[0] = this->dimensions_[0] - 1;
      mod_[1] = this->dimensions_[0] * this->dimensions_[1] - 1;
      div_[0] = msb[0];
      div_[1] = msb[0] + msb[1];
      isAccelerated_ = true;
    }
  } else if(dimensionality_ == 2) {
    bool const isDi = isPowerOfTwo(this->dimensions_[Di_], msb[Di_]);
    bool const isDj = isPowerOfTwo(this->dimensions_[this->Dj_], msb[this->Dj_]);
    bool const allDimensionsArePowerOfTwo = (isDi and isDj);

    if(allDimensionsArePowerOfTwo) {
      mod_[0] = this->dimensions_[Di_] - 1;
      div_[0] = msb[Di_];
      isAccelerated_ = true;
    }
  }

  if(isAccelerated_) {
    this->printMsg("Accelerated getVertex*() requests.", debug::Priority::INFO);
  }

  return 0;
}

template<size_t card>
bool ttk::PeriodicImplicitTriangulation<card>::isPowerOfTwo(unsigned long long int v,
                                                 unsigned long long int &r) {
  if(v && !(v & (v - 1))) {
    r = 0;
    while(v >>= 1)
      r++;
    return true;
  }
  return false;
}

template<size_t card>
bool ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(
  isVertexOnBoundary)(const ttk::SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#else
  TTK_FORCE_USE(vertexId);
#endif

  return false;
}

template<size_t card>
bool ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(
  isEdgeOnBoundary)(const ttk::SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= this->edgeNumber_)
    return -1;
#else
  TTK_FORCE_USE(edgeId);
#endif

  return false;
}

template<size_t card>
bool ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(
  isTriangleOnBoundary)(const ttk::SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
#else
  TTK_FORCE_USE(triangleId);
#endif

  return false;
}

template<size_t card>
ttk::SimplexId ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(
  getVertexNeighborNumber)(const ttk::SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#else
  TTK_FORCE_USE(vertexId);
#endif

  if(dimensionality_ == 3) {
    return 14; // abcdefgh
  } else if(dimensionality_ == 2) {
    return 6; // abcd
  } else if(dimensionality_ == 1) {
    return 2; // ab
  }

  return -1;
}

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  getVertexNeighbor)(const ttk::SimplexId &vertexId,
                     const int &localNeighborId,
                     ttk::SimplexId &neighborId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localNeighborId < 0
     or localNeighborId >= getVertexNeighborNumber(vertexId))
    return -1;
#endif

  neighborId = -1;
  const auto &p = this->underlying().getVertexCoords(vertexId);

  if(this->dimensionality_ == 3) {
    neighborId = this->getVertexNeighbor3d(p.data(), vertexId, localNeighborId);
  } else if(this->dimensionality_ == 2) {
    neighborId = this->getVertexNeighbor2d(p.data(), vertexId, localNeighborId);
  } else if(this->dimensionality_ == 1) {
    // ab
    if(vertexId > 0 and vertexId < this->nbvoxels_[this->Di_]) {
      if(localNeighborId == 0)
        neighborId = vertexId + 1;
      else
        neighborId = vertexId - 1;
    } else if(vertexId == 0) {
      if(localNeighborId == 0)
        neighborId = vertexId + 1; // a
      else
        neighborId = this->nbvoxels_[this->Di_];
    } else {
      if(localNeighborId == 0)
        neighborId = 0; // a
      else
        neighborId = vertexId - 1; // b
    }
  }

  return 0;
}

template<size_t card>
const std::vector<std::vector<ttk::SimplexId>> *
  ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(
    getVertexNeighbors)() {
  if(this->vertexNeighborList_.empty()) {
    Timer t;
    this->vertexNeighborList_.resize(vertexNumber_);
    for(ttk::SimplexId i = 0; i < vertexNumber_; ++i) {
      this->vertexNeighborList_[i].resize(getVertexNeighborNumber(i));
      for(ttk::SimplexId j = 0; j < (ttk::SimplexId)this->vertexNeighborList_[i].size(); ++j)
        this->getVertexNeighbor(i, j, this->vertexNeighborList_[i][j]);
    }

    this->printMsg("Built " + std::to_string(vertexNumber_) + " vertex neighbors.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->vertexNeighborList_;
}

template<size_t card>
ttk::SimplexId ttk::PeriodicImplicitTriangulation<card>::getVertexEdgeNumberInternal(
  const ttk::SimplexId &vertexId) const {
  return getVertexNeighborNumber(vertexId);
}

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::getVertexEdgeInternal(
  const ttk::SimplexId &vertexId, const int &localEdgeId, ttk::SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localEdgeId < 0 or localEdgeId >= getVertexEdgeNumberInternal(vertexId))
    return -1;
#endif

  edgeId = -1;
  const auto &p = this->underlying().getVertexCoords(vertexId);

  if(this->dimensionality_ == 3) {
    edgeId = this->getVertexEdge3d(p.data(), localEdgeId);
  } else if(this->dimensionality_ == 2) {
    edgeId = this->getVertexEdge2d(p.data(), localEdgeId);
  } else if(this->dimensionality_ == 1) {
    if(vertexId > 0 and vertexId < this->nbvoxels_[this->Di_]) {
      // ab
      edgeId = localEdgeId == 0 ? vertexId : vertexId - 1;
    } else if(vertexId == 0) {
      // a
      edgeId = localEdgeId == 0 ? vertexId : 0;
    } else {
      // b
      edgeId = localEdgeId == 0 ? 0 : vertexId - 1;
    }
  }

  return 0;
}

template<size_t card>
const std::vector<std::vector<ttk::SimplexId>> *
  ttk::PeriodicImplicitTriangulation<card>::getVertexEdgesInternal() {
  if(this->vertexEdgeList_.empty()) {
    Timer t;

    this->vertexEdgeList_.resize(vertexNumber_);
    for(ttk::SimplexId i = 0; i < vertexNumber_; ++i) {
      this->vertexEdgeList_[i].resize(getVertexEdgeNumberInternal(i));
      for(ttk::SimplexId j = 0; j < (ttk::SimplexId)this->vertexEdgeList_[i].size(); ++j)
        this->getVertexEdgeInternal(i, j, this->vertexEdgeList_[i][j]);
    }

    this->printMsg("Built " + std::to_string(vertexNumber_) + " vertex edges.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->vertexEdgeList_;
}

template<size_t card>
ttk::SimplexId ttk::PeriodicImplicitTriangulation<card>::getVertexTriangleNumberInternal(
  const ttk::SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#else
  TTK_FORCE_USE(vertexId);
#endif

  if(this->dimensionality_ == 3) {
    return 36; // abcdefgh
  }

  return 0;
}

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::getVertexTriangleInternal(
  const ttk::SimplexId &vertexId,
  const int &localTriangleId,
  ttk::SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localTriangleId < 0
     or localTriangleId >= getVertexTriangleNumberInternal(vertexId))
    return -1;
#endif
  triangleId = -1;

  const auto &p = this->underlying().getVertexCoords(vertexId);

  if(this->dimensionality_ == 3) {
    triangleId = this->getVertexTriangle3d(p.data(), localTriangleId);
  }

  return 0;
}

template<size_t card>
const std::vector<std::vector<ttk::SimplexId>> *
  ttk::PeriodicImplicitTriangulation<card>::getVertexTrianglesInternal() {
  if(this->vertexTriangleList_.empty()) {
    Timer t;

    this->vertexTriangleList_.resize(vertexNumber_);
    for(ttk::SimplexId i = 0; i < vertexNumber_; ++i) {
      this->vertexTriangleList_[i].resize(getVertexTriangleNumberInternal(i));
      for(ttk::SimplexId j = 0; j < (ttk::SimplexId)this->vertexTriangleList_[i].size(); ++j)
        this->getVertexTriangleInternal(i, j, this->vertexTriangleList_[i][j]);
    }

    this->printMsg("Built " + std::to_string(vertexNumber_) + " vertex triangles.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->vertexTriangleList_;
}

template<size_t card>
ttk::SimplexId ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(
  getVertexLinkNumber)(const ttk::SimplexId &vertexId) const {
  return getVertexStarNumber(vertexId);
}

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  getVertexLink)(const ttk::SimplexId &vertexId,
                 const int &localLinkId,
                 ttk::SimplexId &linkId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localLinkId < 0 or localLinkId >= getVertexLinkNumber(vertexId))
    return -1;
#endif

  linkId = -1;
  const auto &p = this->underlying().getVertexCoords(vertexId);

  if(this->dimensionality_ == 3) {
    linkId = this->getVertexLink3d(p.data(), localLinkId);
  } else if(this->dimensionality_ == 2) {
    linkId = this->getVertexLink2d(p.data(), localLinkId); // abcd
  }

  return 0;
}

template<size_t card>
const std::vector<std::vector<ttk::SimplexId>> *
  ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getVertexLinks)() {
  if(this->vertexLinkList_.empty()) {
    Timer t;

    this->vertexLinkList_.resize(vertexNumber_);
    for(ttk::SimplexId i = 0; i < vertexNumber_; ++i) {
      this->vertexLinkList_[i].resize(getVertexLinkNumber(i));
      for(ttk::SimplexId j = 0; j < (ttk::SimplexId)this->vertexLinkList_[i].size(); ++j)
        this->getVertexLink(i, j, this->vertexLinkList_[i][j]);
    }

    this->printMsg("Built " + std::to_string(vertexNumber_) + " vertex links.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->vertexLinkList_;
}

template<size_t card>
ttk::SimplexId ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(
  getVertexStarNumber)(const ttk::SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#else
  TTK_FORCE_USE(vertexId);
#endif

  if(this->dimensionality_ == 3) {
    return 24; // abcdefgh
  } else if(this->dimensionality_ == 2) {
    return 6; // abcd
  }

  return 0;
}

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  getVertexStar)(const ttk::SimplexId &vertexId,
                 const int &localStarId,
                 ttk::SimplexId &starId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localStarId < 0 or localStarId >= getVertexStarNumber(vertexId))
    return -1;
#endif

  starId = -1;
  const auto &p = this->underlying().getVertexCoords(vertexId);

  if(this->dimensionality_ == 3) {
    starId = this->getVertexStar3d(p.data(), localStarId);
  } else if(this->dimensionality_ == 2) {
    starId = this->getVertexStar2d(p.data(), localStarId);
  }

  return 0;
}

template<size_t card>
const std::vector<std::vector<ttk::SimplexId>> *
  ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getVertexStars)() {
  if(this->vertexStarList_.empty()) {
    Timer t;
    this->vertexStarList_.resize(vertexNumber_);
    for(ttk::SimplexId i = 0; i < vertexNumber_; ++i) {
      this->vertexStarList_[i].resize(getVertexStarNumber(i));
      for(ttk::SimplexId j = 0; j < (ttk::SimplexId)this->vertexStarList_[i].size(); ++j)
        this->getVertexStar(i, j, this->vertexStarList_[i][j]);
    }

    this->printMsg("Built " + std::to_string(vertexNumber_) + " vertex stars.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->vertexStarList_;
}

#ifdef TTK_ENABLE_MPI

template<size_t card>
void ttk::PeriodicImplicitTriangulation<card>::setIsBoundaryPeriodic(
  std::array<unsigned char, 6> boundary) {
  this->isBoundaryPeriodic = boundary;
}

template<size_t card>
void ttk::PeriodicImplicitTriangulation<card>::createMetaGrid(const double *const bounds) {
  // only works with 2 processes or more
  if(!ttk::isRunningWithMPI()) {
    return;
  }

  // no need to create it anew?
  if(this->metaGrid_ != nullptr) {
    return;
  }

  // Reorganize bounds to only execute Allreduce twice
  std::array<double, 6> tempBounds = {
    bounds[0], bounds[2], bounds[4], bounds[1], bounds[3], bounds[5],
  };

  for(int i = 0; i < 3; i++) {
    if(this->dimensionality_ > i) {
      if(this->isBoundaryPeriodic[2 * i] == 1) {
        tempBounds[i] += this->spacing_[i];
      }
      if(this->isBoundaryPeriodic[2 * i + 1] == 1) {
        tempBounds[3 + i] -= this->spacing_[i];
      }
    }
  }

  std::array<double, 6> tempGlobalBounds{};
  // Compute and send to all processes the lower bounds of the data set
  MPI_Allreduce(tempBounds.data(), tempGlobalBounds.data(), 3, MPI_DOUBLE,
                MPI_MIN, ttk::MPIcomm_);
  // Compute and send to all processes the higher bounds of the data set
  MPI_Allreduce(&tempBounds[3], &tempGlobalBounds[3], 3, MPI_DOUBLE, MPI_MAX,
                ttk::MPIcomm_);

  // re-order tempGlobalBounds
  std::array<double, 6> globalBounds{
    tempGlobalBounds[0], tempGlobalBounds[3], tempGlobalBounds[1],
    tempGlobalBounds[4], tempGlobalBounds[2], tempGlobalBounds[5],
  };

  const std::array<ttk::SimplexId, 3> dimensions = {
    static_cast<ttk::SimplexId>(
      std::round((globalBounds[1] - globalBounds[0]) / this->spacing_[0]))
      + 1,
    static_cast<ttk::SimplexId>(
      std::round((globalBounds[3] - globalBounds[2]) / this->spacing_[1]))
      + 1,
    static_cast<ttk::SimplexId>(
      std::round((globalBounds[5] - globalBounds[4]) / this->spacing_[2]))
      + 1,
  };

  this->localGridOffset_ = {
    static_cast<ttk::SimplexId>(
      std::round((this->origin_[0] - globalBounds[0]) / this->spacing_[0])),
    static_cast<ttk::SimplexId>(
      std::round((this->origin_[1] - globalBounds[2]) / this->spacing_[1])),
    static_cast<ttk::SimplexId>(
      std::round((this->origin_[2] - globalBounds[4]) / this->spacing_[2])),
  };

  this->metaGrid_ = std::make_shared<PeriodicNoPreconditions<card>>();
  this->metaGrid_->setInputGrid(globalBounds[0], globalBounds[1],
                                globalBounds[2], this->spacing_[0],
                                this->spacing_[1], this->spacing_[2],
                                dimensions[0], dimensions[1], dimensions[2]);
}

template<size_t card>
int ttk::PeriodicImplicitTriangulation<card>::preconditionDistributedCells() {
  if(this->hasPreconditionedDistributedCells_) {
    return 0;
  }
  if(!ttk::hasInitializedMPI()) {
    return -1;
  }
  if(this->metaGrid_ == nullptr) {
    return 0;
  }
  if(this->cellGhost_ == nullptr) {
    if(ttk::isRunningWithMPI()) {
      this->printErr("Missing cell ghost array!");
    }
    return -3;
  }

  Timer tm{};

  // there are 6 tetrahedra per cubic cell (and 2 triangles per square)
  const int nTetraPerCube{this->dimensionality_ == 3 ? 6 : 2};

  // number of local cells (with ghost cells but without the additional periodic
  // cells)
  const auto nLocCells{(this->dimensions_[0] - 1) * (this->dimensions_[1] - 1)
                       * (this->dimensions_[2] - 1) * nTetraPerCube};

  std::vector<unsigned char> fillCells(nLocCells / nTetraPerCube);

  this->neighborCellBBoxes_.resize(ttk::MPIsize_);
  auto &localBBox{this->neighborCellBBoxes_[ttk::MPIrank_]};
  // "good" starting values?
  localBBox = {
    this->localGridOffset_[0] + this->dimensions_[0], this->localGridOffset_[0],
    this->localGridOffset_[1] + this->dimensions_[1], this->localGridOffset_[1],
    this->localGridOffset_[2] + this->dimensions_[2], this->localGridOffset_[2],
  };
  const auto &dims{this->metaGrid_->getGridDimensions()};

  for(ttk::SimplexId lcid = 0; lcid < nLocCells; ++lcid) {
    // only keep non-ghost cells
    if(this->cellGhost_[lcid / nTetraPerCube] != 0) {
      continue;
    }
    // local vertex coordinates
    std::array<ttk::SimplexId, 3> p{};
    if(this->dimensionality_ == 3) {
      this->tetrahedronToPosition(lcid, p.data());
    } else if(this->dimensionality_ == 2) {
      this->triangleToPosition2d(lcid, p.data());
      // compatibility with tetrahedronToPosition; fix a bounding box
      // error in the first axis
      p[0] /= 2;
    }

    // global vertex coordinates
    p[0] += this->localGridOffset_[0];
    p[1] += this->localGridOffset_[1];
    p[2] += this->localGridOffset_[2];

    if(p[0] < localBBox[0]) {
      localBBox[0] = std::max(p[0], static_cast<ttk::SimplexId>(0));
    }
    if(p[0] > localBBox[1]) {
      localBBox[1] = min(p[0], dims[0]);
    }
    if(p[1] < localBBox[2]) {
      localBBox[2] = std::max(p[1], static_cast<ttk::SimplexId>(0));
    }
    if(p[1] > localBBox[3]) {
      localBBox[3] = std::min(p[1], dims[1]);
    }
    if(p[2] < localBBox[4]) {
      localBBox[4] = std::max(p[2], static_cast<ttk::SimplexId>(0));
    }
    if(p[2] > localBBox[5]) {
      localBBox[5] = std::min(p[2], dims[2]);
    }
  }
  localBBox[0] -= isBoundaryPeriodic[0];
  if(this->dimensionality_ > 1) {
    localBBox[2] -= isBoundaryPeriodic[2];
    if(this->dimensionality_ > 2)
      localBBox[4] -= isBoundaryPeriodic[4];
  }

  localBBox[1]++;
  localBBox[3]++;
  localBBox[5]++;

  for(size_t i = 0; i < this->neighborRanks_.size(); ++i) {
    const auto neigh{this->neighborRanks_[i]};
    MPI_Sendrecv(this->neighborCellBBoxes_[ttk::MPIrank_].data(), 6,
                 ttk::getMPIType(ttk::SimplexId{}), neigh, ttk::MPIrank_,
                 this->neighborCellBBoxes_[neigh].data(), 6,
                 ttk::getMPIType(ttk::SimplexId{}), neigh, neigh, ttk::MPIcomm_,
                 MPI_STATUS_IGNORE);
  }

  this->hasPreconditionedDistributedCells_ = true;

  return 0;
}

template<size_t card>
std::array<ttk::SimplexId, 3> ttk::PeriodicImplicitTriangulation<card>::getVertGlobalCoords(
  const ttk::SimplexId lvid) const {
  // local vertex coordinates
  std::array<ttk::SimplexId, 3> p{};
  if(this->dimensionality_ == 3) {
    this->vertexToPosition(lvid, p.data());
  } else if(this->dimensionality_ == 2) {
    this->vertexToPosition2d(lvid, p.data());
  }
  // global vertex coordinates
  p[0] += this->localGridOffset_[0];
  p[1] += this->localGridOffset_[1];
  p[2] += this->localGridOffset_[2];

  const auto &dims{this->metaGrid_->getGridDimensions()};

  p[0] = (p[0] + dims[0]) % dims[0];
  if(this->dimensionality_ > 1) {
    p[1] = (p[1] + dims[1]) % dims[1];
    if(this->dimensionality_ > 2)
      p[2] = (p[2] + dims[2]) % dims[2];
  }

  return p;
}

template<size_t card>
std::array<ttk::SimplexId, 3> ttk::PeriodicImplicitTriangulation<card>::getVertLocalCoords(
  const ttk::SimplexId gvid) const {
  // global vertex coordinates
  std::array<ttk::SimplexId, 3> pGlobal{};
  if(this->dimensionality_ == 3) {
    this->metaGrid_->vertexToPosition(gvid, pGlobal.data());
  } else if(this->dimensionality_ == 2) {
    this->metaGrid_->vertexToPosition2d(gvid, pGlobal.data());
  }
  std::array<ttk::SimplexId, 3> p{pGlobal};
  // local vertex coordinates
  p[0] -= this->localGridOffset_[0];
  p[1] -= this->localGridOffset_[1];
  p[2] -= this->localGridOffset_[2];

  const auto &dims{this->getGridDimensions()};

  if(p[0] >= 0 && p[1] >= 0 && p[2] >= 0 && p[0] <= dims[0] - 1
     && p[1] <= dims[1] - 1 && p[2] <= dims[2] - 1) {
    return p;
  }
  for(int i = 0; i < 3; i++) {
    if((p[i] < 0 || p[i] > dims[i] - 1) && pGlobal[i] == 0) {
      p[i] = dims[i] - 1;
    }
    if((p[i] < 0 || p[i] > dims[i] - 1)
       && pGlobal[i] == this->metaGrid_->dimensions_[i] - 1) {
      p[i] = 0;
    }
  }

  if(p[0] >= 0 && p[1] >= 0 && p[2] >= 0 && p[0] <= dims[0] - 1
     && p[1] <= dims[1] - 1 && p[2] <= dims[2] - 1) {
    if(this->vertexGhost_[p[0] + p[1] * dims[0] + p[2] * dims[0] * dims[1]]
       != 0) {
      return p;
    }
  }
  return std::array<ttk::SimplexId, 3>{-1, -1, -1};
}

template<size_t card>
int ttk::PeriodicImplicitTriangulation<card>::getCellRankInternal(
  const ttk::SimplexId lcid) const {

  const int nTetraPerCube{this->dimensionality_ == 3 ? 6 : 2};
  const auto locCubeId{lcid / nTetraPerCube};

  if(this->cellGhost_[locCubeId] == 0) {
    return ttk::MPIrank_;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(this->neighborRanks_.empty()) {
    this->printErr("Empty neighborsRanks_!");
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  const auto nVertsCell{this->getCellVertexNumber(lcid)};
  std::vector<bool> inRank(nVertsCell);
  std::map<int, int> neighborOccurrences;
  for(const auto neigh : this->neighborRanks_) {
    std::fill(inRank.begin(), inRank.end(), false);
    const auto &bbox{this->neighborCellBBoxes_[neigh]};
    for(ttk::SimplexId i = 0; i < nVertsCell; ++i) {
      ttk::SimplexId v{};
      this->getCellVertex(lcid, i, v);
      if(this->vertexGhost_[v] == 0) {
        inRank[i] = true;
      } else {
        const auto p{this->getVertGlobalCoords(v)};
        if(p[0] >= bbox[0] && p[0] <= bbox[1] && p[1] >= bbox[2]
           && p[1] <= bbox[3] && p[2] >= bbox[4] && p[2] <= bbox[5]) {
          inRank[i] = true;
        }
      }
    }
    if(std::all_of(
         inRank.begin(), inRank.end(), [](const bool v) { return v; })) {
      return neigh;
    }
    neighborOccurrences[neigh]
      = std::accumulate(inRank.begin(), inRank.end(), 0);
  }

  auto pr = std::max_element(
    std::begin(neighborOccurrences), std::end(neighborOccurrences),
    [](const std::pair<int, int> &p1, const std::pair<int, int> &p2) {
      return p1.second < p2.second;
    });
  return pr->first;
}

#endif // TTK_ENABLE_MPI

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  getVertexPoint)(const ttk::SimplexId &vertexId,
                  float &x,
                  float &y,
                  float &z) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#endif

  if(this->dimensionality_ == 3) {
    const auto &p = this->underlying().getVertexCoords(vertexId);

    x = this->origin_[0] + this->spacing_[0] * p[0];
    y = this->origin_[1] + this->spacing_[1] * p[1];
    z = this->origin_[2] + this->spacing_[2] * p[2];
  } else if(this->dimensionality_ == 2) {
    const auto &p = this->underlying().getVertexCoords(vertexId);

    if(this->dimensions_[0] > 1 and this->dimensions_[1] > 1) {
      x = this->origin_[0] + this->spacing_[0] * p[0];
      y = this->origin_[1] + this->spacing_[1] * p[1];
      z = this->origin_[2];
    } else if(this->dimensions_[1] > 1 and this->dimensions_[2] > 1) {
      x = this->origin_[0];
      y = this->origin_[1] + this->spacing_[1] * p[0];
      z = this->origin_[2] + this->spacing_[2] * p[1];
    } else if(this->dimensions_[0] > 1 and this->dimensions_[2] > 1) {
      x = this->origin_[0] + this->spacing_[0] * p[0];
      y = this->origin_[1];
      z = this->origin_[2] + this->spacing_[2] * p[1];
    }
  } else if(this->dimensionality_ == 1) {
    if(this->dimensions_[0] > 1) {
      x = this->origin_[0] + this->spacing_[0] * vertexId;
      y = this->origin_[1];
      z = this->origin_[2];
    } else if(this->dimensions_[1] > 1) {
      x = this->origin_[0];
      y = this->origin_[1] + this->spacing_[1] * vertexId;
      z = this->origin_[2];
    } else if(this->dimensions_[2] > 1) {
      x = this->origin_[0];
      y = this->origin_[1];
      z = this->origin_[2] + this->spacing_[2] * vertexId;
    }
  }

  return 0;
}

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::getEdgeVertexInternal(
  const ttk::SimplexId &edgeId,
  const int &localVertexId,
  ttk::SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= this->edgeNumber_)
    return -1;
  if(localVertexId < 0 or localVertexId >= 2)
    return -2;
#endif

  using CRTPClass = ttk::PeriodicImplicitTriangulationCRTP<card, Derived>;
  vertexId = -1;
  const auto &p = this->underlying().getEdgeCoords(edgeId);
  const ttk::SimplexId wrapXRight = (p[0] == this->nbvoxels_[0] ? -this->wrap_[0] : 0);
  const ttk::SimplexId wrapYBottom = (p[1] == this->nbvoxels_[1] ? -this->wrap_[1] : 0);
  const ttk::SimplexId wrapZFront = (p[2] == this->nbvoxels_[2] ? -this->wrap_[2] : 0);
  const auto a = p[0] + this->underlying().getEdgeVertexAccelerated(edgeId);

  switch(this->underlying().getEdgePosition(edgeId)) {
    case CRTPClass::EdgePosition::L_3D:
      vertexId = a + (localVertexId == 0 ? 0 : (1 + wrapXRight));
      break;
    case CRTPClass::EdgePosition::H_3D:
      vertexId = a + (localVertexId == 0 ? 0 : (this->vshift_[0] + wrapYBottom));
      break;
    case CRTPClass::EdgePosition::P_3D:
      vertexId = a + (localVertexId == 0 ? 0 : (this->vshift_[1] + wrapZFront));
      break;
    case CRTPClass::EdgePosition::D1_3D:
      vertexId = a
                 + (localVertexId == 0 ? (1 + wrapXRight)
                                       : (this->vshift_[0] + wrapYBottom));
      break;
    case CRTPClass::EdgePosition::D2_3D:
      vertexId = a
                 + (localVertexId == 0
                      ? 0
                      : (this->vshift_[0] + wrapYBottom + this->vshift_[1] + wrapZFront));
      break;
    case CRTPClass::EdgePosition::D3_3D:
      vertexId
        = a
          + (localVertexId == 0 ? (1 + wrapXRight) : (this->vshift_[1] + wrapZFront));
      break;
    case CRTPClass::EdgePosition::D4_3D:
      vertexId = a
                 + (localVertexId == 0
                      ? (1 + wrapXRight)
                      : (this->vshift_[0] + wrapYBottom + this->vshift_[1] + wrapZFront));
      break;
    case CRTPClass::EdgePosition::L_2D:
      vertexId = a + (localVertexId == 0 ? 0 : (1 + wrapXRight));
      break;
    case CRTPClass::EdgePosition::H_2D:
      vertexId = a + (localVertexId == 0 ? 0 : (this->vshift_[0] + wrapYBottom));
      break;
    case CRTPClass::EdgePosition::D1_2D:
      vertexId = a
                 + (localVertexId == 0 ? (1 + wrapXRight)
                                       : (this->vshift_[0] + wrapYBottom));
      break;
    case CRTPClass::EdgePosition::FIRST_EDGE_1D:
      vertexId = localVertexId == 0 ? 0 : 1;
      break;
    case CRTPClass::EdgePosition::LAST_EDGE_1D:
      vertexId = localVertexId == 0 ? edgeId : 0;
      break;
    case CRTPClass::EdgePosition::CENTER_1D:
      vertexId = localVertexId == 0 ? edgeId : edgeId + 1;
      break;
    default:
      break;
  }

  return 0;
}

template<size_t card>
const std::vector<std::array<ttk::SimplexId, 2>> *
  ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getEdges)() {

  if(this->edgeList_.empty()) {
    Timer t;

    this->edgeList_.resize(this->edgeNumber_);
    for(ttk::SimplexId i = 0; i < this->edgeNumber_; ++i) {
      ttk::SimplexId id0, id1;
      this->getEdgeVertexInternal(i, 0, id0);
      this->getEdgeVertexInternal(i, 1, id1);
      this->edgeList_[i] = {id0, id1};
    }

    this->printMsg(
      "Built " + std::to_string(this->edgeNumber_) + " edges.", 1, t.getElapsedTime(), 1);
  }

  return &this->edgeList_;
}

template <size_t card, typename Derived>
ttk::SimplexId
  ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::getEdgeTriangleNumberInternal(
    const ttk::SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= this->edgeNumber_)
    return -1;
#endif

  using CRTPClass = ttk::PeriodicImplicitTriangulationCRTP<card, Derived>;
  //using EdgePosition = typename ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::EdgePosition;
  switch(this->underlying().getEdgePosition(edgeId)) {
    case CRTPClass::EdgePosition::L_3D:
    case CRTPClass::EdgePosition::H_3D:
    case CRTPClass::EdgePosition::P_3D:
    case CRTPClass::EdgePosition::D4_3D:
      return 6;
    case CRTPClass::EdgePosition::D1_3D:
    case CRTPClass::EdgePosition::D2_3D:
    case CRTPClass::EdgePosition::D3_3D:
      return 4;
    case CRTPClass::EdgePosition::L_2D:
    case CRTPClass::EdgePosition::H_2D:
    case CRTPClass::EdgePosition::D1_2D:
      return 2;
    default:
      return 0;
  }
}

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::getEdgeTriangleInternal(
  const ttk::SimplexId &edgeId,
  const int &localTriangleId,
  ttk::SimplexId &triangleId) const {

  using CRTPClass = ttk::PeriodicImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(localTriangleId < 0
     or localTriangleId >= getEdgeTriangleNumberInternal(edgeId))
    return -1;
#endif

  triangleId = -1;
  const auto &p = this->underlying().getEdgeCoords(edgeId);

  switch(this->underlying().getEdgePosition(edgeId)) {
    case CRTPClass::EdgePosition::L_3D:
      triangleId = this->getEdgeTriangle3dL(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::H_3D:
      triangleId = this->getEdgeTriangle3dH(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::P_3D:
      triangleId = this->getEdgeTriangle3dP(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::D1_3D:
      triangleId = this->getEdgeTriangle3dD1(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::D2_3D:
      triangleId = this->getEdgeTriangle3dD2(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::D3_3D:
      triangleId = this->getEdgeTriangle3dD3(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::D4_3D:
      triangleId = this->getEdgeTriangle3dD4(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::L_2D:
      triangleId = this->getEdgeTriangle2dL(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::H_2D:
      triangleId = this->getEdgeTriangle2dH(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::D1_2D:
      triangleId = this->getEdgeTriangle2dD1(p.data(), localTriangleId);
      break;
    default:
      break;
  }

  return 0;
}

template<size_t card>
const std::vector<std::vector<ttk::SimplexId>> *
  ttk::PeriodicImplicitTriangulation<card>::getEdgeTrianglesInternal() {
  if(this->edgeTriangleList_.empty()) {
    Timer t;

    this->edgeTriangleList_.resize(this->edgeNumber_);
    for(ttk::SimplexId i = 0; i < this->edgeNumber_; ++i) {
      this->edgeTriangleList_[i].resize(this->getEdgeTriangleNumberInternal(i));
      for(ttk::SimplexId j = 0; j < (ttk::SimplexId)this->edgeTriangleList_[i].size(); ++j)
        this->getEdgeTriangleInternal(i, j, this->edgeTriangleList_[i][j]);
    }

    this->printMsg("Built " + std::to_string(this->edgeNumber_) + " edge triangles.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->edgeTriangleList_;
}

template<size_t card>
ttk::SimplexId ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(
  getEdgeLinkNumber)(const ttk::SimplexId &edgeId) const {

  return this->getEdgeStarNumber(edgeId);
}

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  getEdgeLink)(const ttk::SimplexId &edgeId,
               const int &localLinkId,
               ttk::SimplexId &linkId) const {

  using CRTPClass = ttk::PeriodicImplicitTriangulationCRTP<card, Derived>;
  //using EdgePosition = typename ttk::PeriodicImplicitTriangulationCRTP<card,ttk::PeriodicImplicitTriangulation<card>>::EdgePosition;
#ifndef TTK_ENABLE_KAMIKAZE
  if(localLinkId < 0 or localLinkId >= getEdgeLinkNumber(edgeId))
    return -1;
#endif

  linkId = -1;
  const auto &p = this->underlying().getEdgeCoords(edgeId);

  switch(this->underlying().getEdgePosition(edgeId)) {
    case CRTPClass::EdgePosition::L_3D:
      linkId = this->getEdgeLinkL(p.data(), localLinkId);
      break;
    case CRTPClass::EdgePosition::H_3D:
      linkId = this->getEdgeLinkH(p.data(), localLinkId);
      break;
    case CRTPClass::EdgePosition::P_3D:
      linkId = this->getEdgeLinkP(p.data(), localLinkId);
      break;
    case CRTPClass::EdgePosition::D1_3D:
      linkId = this->getEdgeLinkD1(p.data(), localLinkId);
      break;
    case CRTPClass::EdgePosition::D2_3D:
      linkId = this->getEdgeLinkD2(p.data(), localLinkId);
      break;
    case CRTPClass::EdgePosition::D3_3D:
      linkId = this->getEdgeLinkD3(p.data(), localLinkId);
      break;
    case CRTPClass::EdgePosition::D4_3D:
      linkId = this->getEdgeLinkD4(p.data(), localLinkId);
      break;
    case CRTPClass::EdgePosition::L_2D:
      linkId = this->getEdgeLink2dL(p.data(), localLinkId);
      break;
    case CRTPClass::EdgePosition::H_2D:
      linkId = this->getEdgeLink2dH(p.data(), localLinkId);
      break;
    case CRTPClass::EdgePosition::D1_2D:
      linkId = this->getEdgeLink2dD1(p.data(), localLinkId);
      break;
    default:
      break;
  }

  return 0;
}

template<size_t card>
const std::vector<std::vector<ttk::SimplexId>> *
  ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getEdgeLinks)() {

  if(this->edgeLinkList_.empty()) {
    Timer t;

    this->edgeLinkList_.resize(this->edgeNumber_);
    for(ttk::SimplexId i = 0; i < this->edgeNumber_; ++i) {
      this->edgeLinkList_[i].resize(this->getEdgeLinkNumber(i));
      for(ttk::SimplexId j = 0; j < (ttk::SimplexId)this->edgeLinkList_[i].size(); ++j)
        this->getEdgeLink(i, j, this->edgeLinkList_[i][j]); //TODO PAS SÛR cf argument ttk...internal
    }

    this->printMsg("Built " + std::to_string(this->edgeNumber_) + " edge links.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->edgeLinkList_;
}

template <size_t card, typename Derived>
ttk::SimplexId
  ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
    getEdgeStarNumber)(const ttk::SimplexId &edgeId) const {

  using CRTPClass = ttk::PeriodicImplicitTriangulationCRTP<card, Derived>;
  //using EdgePosition = typename ttk::PeriodicImplicitTriangulationCRTP<card,ttk::PeriodicImplicitTriangulation<card>>::EdgePosition;
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= this->edgeNumber_)
    return -1;
#endif

  switch(this->underlying().getEdgePosition(edgeId)) {
    case CRTPClass::EdgePosition::L_3D:
    case CRTPClass::EdgePosition::H_3D:
    case CRTPClass::EdgePosition::P_3D:
    case CRTPClass::EdgePosition::D4_3D:
      return 6;
    case CRTPClass::EdgePosition::D1_3D:
    case CRTPClass::EdgePosition::D2_3D:
    case CRTPClass::EdgePosition::D3_3D:
      return 4;
    case CRTPClass::EdgePosition::L_2D:
    case CRTPClass::EdgePosition::H_2D:
    case CRTPClass::EdgePosition::D1_2D:
      return 2;
    default:
      return 0;
  }

  return 0;
}

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  getEdgeStar)(const ttk::SimplexId &edgeId,
               const int &localStarId,
               ttk::SimplexId &starId) const {
  using CRTPClass = ttk::PeriodicImplicitTriangulationCRTP<card, Derived>;
  //using EdgePosition = typename ttk::PeriodicImplicitTriangulationCRTP<card,ttk::PeriodicImplicitTriangulation<card>>::EdgePosition;
#ifndef TTK_ENABLE_KAMIKAZE
  if(localStarId < 0 or localStarId >= getEdgeStarNumber(edgeId))
    return -1;
#endif

  starId = -1;
  const auto &p = this->underlying().getEdgeCoords(edgeId);

  switch(this->underlying().getEdgePosition(edgeId)) {
    case CRTPClass::EdgePosition::L_3D:
      starId = this->getEdgeStarL(p.data(), localStarId);
      break;
    case CRTPClass::EdgePosition::H_3D:
      starId = this->getEdgeStarH(p.data(), localStarId);
      break;
    case CRTPClass::EdgePosition::P_3D:
      starId = this->getEdgeStarP(p.data(), localStarId);
      break;
    case CRTPClass::EdgePosition::D1_3D:
      starId = this->getEdgeStarD1(p.data(), localStarId);
      break;
    case CRTPClass::EdgePosition::D2_3D:
      starId = this->getEdgeStarD2(p.data(), localStarId);
      break;
    case CRTPClass::EdgePosition::D3_3D:
      starId = this->getEdgeStarD3(p.data(), localStarId);
      break;
    case CRTPClass::EdgePosition::D4_3D:
      starId
        = p[2] * this->tshift_[1] + p[1] * this->tshift_[0] + p[0] * 6 + localStarId;
      break;
    case CRTPClass::EdgePosition::L_2D:
      starId = this->getEdgeStar2dL(p.data(), localStarId);
      break;
    case CRTPClass::EdgePosition::H_2D:
      starId = this->getEdgeStar2dH(p.data(), localStarId);
      break;
    case CRTPClass::EdgePosition::D1_2D:
      starId = p[0] * 2 + p[1] * this->tshift_[0] + localStarId;
      break;
    default:
      break;
  }

  return 0;
}

template<size_t card>
const std::vector<std::vector<ttk::SimplexId>> *
  ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getEdgeStars)() {

  if(this->edgeStarList_.empty()) {
    Timer t;

    this->edgeStarList_.resize(this->edgeNumber_);
    for(ttk::SimplexId i = 0; i < this->edgeNumber_; ++i) {
      this->edgeStarList_[i].resize(this->getEdgeStarNumber(i));
      for(ttk::SimplexId j = 0; j < (ttk::SimplexId)this->edgeStarList_[i].size(); ++j)
        this->getEdgeStar(i, j, this->edgeStarList_[i][j]);
    }

    this->printMsg("Built " + std::to_string(this->edgeNumber_) + " edge stars.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->edgeStarList_;
}

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::getTriangleVertexInternal(
  const ttk::SimplexId &triangleId,
  const int &localVertexId,
  ttk::SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
  if(localVertexId < 0 or localVertexId >= 3)
    return -2;
#endif

  vertexId = -1;
  const auto &p = this->underlying().getTriangleCoords(triangleId);
  const ttk::SimplexId wrapXRight = (p[0] / 2 == this->nbvoxels_[this->Di_]) ? -this->wrap_[0] : 0;
  const ttk::SimplexId wrapYBottom = (p[1] == this->nbvoxels_[this->Dj_]) ? -this->wrap_[1] : 0;

  using CRTPClass = ttk::PeriodicImplicitTriangulationCRTP<card, Derived>;
  //using TrianglePosition = typename ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::TrianglePosition;
  switch(this->underlying().getTrianglePosition(triangleId)) {
    case CRTPClass::TrianglePosition::F_3D:
      vertexId = this->getTriangleVertexF(p.data(), localVertexId);
      break;
    case CRTPClass::TrianglePosition::H_3D:
      vertexId = this->getTriangleVertexH(p.data(), localVertexId);
      break;
    case CRTPClass::TrianglePosition::C_3D:
      vertexId = this->getTriangleVertexC(p.data(), localVertexId);
      break;
    case CRTPClass::TrianglePosition::D1_3D:
      vertexId = this->getTriangleVertexD1(p.data(), localVertexId);
      break;
    case CRTPClass::TrianglePosition::D2_3D:
      vertexId = this->getTriangleVertexD2(p.data(), localVertexId);
      break;
    case CRTPClass::TrianglePosition::D3_3D:
      vertexId = this->getTriangleVertexD3(p.data(), localVertexId);
      break;
    case CRTPClass::TrianglePosition::TOP_2D:
      if(localVertexId == 0) {
        vertexId = p[0] / 2 + p[1] * this->vshift_[0];
      } else if(localVertexId == 1) {
        vertexId = p[0] / 2 + p[1] * this->vshift_[0] + 1 + wrapXRight;
      } else if(localVertexId == 2) {
        vertexId = p[0] / 2 + p[1] * this->vshift_[0] + this->vshift_[0] + wrapYBottom;
      }
      break;
    case CRTPClass::TrianglePosition::BOTTOM_2D:
      if(localVertexId == 0) {
        vertexId = p[0] / 2 + p[1] * this->vshift_[0] + 1 + wrapXRight;
      } else if(localVertexId == 1) {
        vertexId = p[0] / 2 + p[1] * this->vshift_[0] + this->vshift_[0] + 1 + wrapXRight
                   + wrapYBottom;
      } else if(localVertexId == 2) {
        vertexId = p[0] / 2 + p[1] * this->vshift_[0] + this->vshift_[0] + wrapYBottom;
      }
      break;
    default:
      break;
  }

  return 0;
}

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::getTriangleEdgeInternal(
  const ttk::SimplexId &triangleId,
  const int &localEdgeId,
  ttk::SimplexId &edgeId) const {
  using CRTPClass = ttk::PeriodicImplicitTriangulationCRTP<card, Derived>;
  //using TrianglePosition = typename ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::TrianglePosition;
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
  if(localEdgeId < 0 or localEdgeId >= 3)
    return -2;
#endif

  edgeId = -1;
  const auto &p = this->underlying().getTriangleCoords(triangleId);
  const ttk::SimplexId wrapXRight = (p[0] / 2 == this->nbvoxels_[this->Di_]) ? -this->wrap_[0] : 0;
  const ttk::SimplexId wrapYBottom = (p[1] == this->nbvoxels_[this->Dj_]) ? -this->wrap_[1] : 0;
  const ttk::SimplexId id = triangleId % 2;

  switch(this->underlying().getTrianglePosition(triangleId)) {
    case CRTPClass::TrianglePosition::F_3D:
      edgeId = (id == 1) ? this->getTriangleEdgeF_1(p.data(), localEdgeId)
                         : this->getTriangleEdgeF_0(p.data(), localEdgeId);
      break;
    case CRTPClass::TrianglePosition::H_3D:
      edgeId = (id == 1) ? this->getTriangleEdgeH_1(p.data(), localEdgeId)
                         : this->getTriangleEdgeH_0(p.data(), localEdgeId);
      break;
    case CRTPClass::TrianglePosition::C_3D:
      edgeId = (id == 1) ? this->getTriangleEdgeC_1(p.data(), localEdgeId)
                         : this->getTriangleEdgeC_0(p.data(), localEdgeId);
      break;
    case CRTPClass::TrianglePosition::D1_3D:
      edgeId = (id == 1) ? this->getTriangleEdgeD1_1(p.data(), localEdgeId)
                         : this->getTriangleEdgeD1_0(p.data(), localEdgeId);
      break;
    case CRTPClass::TrianglePosition::D2_3D:
      edgeId = (id == 1) ? this->getTriangleEdgeD2_1(p.data(), localEdgeId)
                         : this->getTriangleEdgeD2_0(p.data(), localEdgeId);
      break;
    case CRTPClass::TrianglePosition::D3_3D:
      edgeId = (id == 1) ? this->getTriangleEdgeD3_1(p.data(), localEdgeId)
                         : this->getTriangleEdgeD3_0(p.data(), localEdgeId);
      break;
    case CRTPClass::TrianglePosition::TOP_2D:
      if(localEdgeId == 0) {
        edgeId = p[0] / 2 + p[1] * this->eshift_[0];
      } else if(localEdgeId == 1) {
        edgeId = this->esetshift_[0] + p[0] / 2 + p[1] * this->eshift_[2];
      } else if(localEdgeId == 2) {
        edgeId = this->esetshift_[1] + p[0] / 2 + p[1] * this->eshift_[4];
      }
      break;
    case CRTPClass::TrianglePosition::BOTTOM_2D:
      if(localEdgeId == 0) {
        edgeId = p[0] / 2 + (p[1] + 1) * this->eshift_[0] + wrapYBottom;
      } else if(localEdgeId == 1) {
        edgeId
          = this->esetshift_[0] + (p[0] + 1) / 2 + p[1] * this->eshift_[2] + wrapXRight;
      } else if(localEdgeId == 2) {
        edgeId = this->esetshift_[1] + p[0] / 2 + p[1] * this->eshift_[4];
      }
      break;
    default:
      break;
  }

  return 0;
}

template<size_t card>
int ttk::PeriodicImplicitTriangulation<card>::getTriangleEdgesInternal(
  std::vector<std::vector<ttk::SimplexId>> &edges) const {
  edges.resize(triangleNumber_);
  for(ttk::SimplexId i = 0; i < triangleNumber_; ++i) {
    edges[i].resize(3);
    for(int j = 0; j < 3; ++j)
      this->getTriangleEdgeInternal(i, j, edges[i][j]);
  }
  return 0;
}

template<size_t card>
const std::vector<std::vector<ttk::SimplexId>> *
  ttk::PeriodicImplicitTriangulation<card>::getTriangleEdgesInternal() {
  if(this->triangleEdgeVector_.empty()) {
    Timer t;

    getTriangleEdgesInternal(this->triangleEdgeVector_);

    this->printMsg("Built " + std::to_string(triangleNumber_) + " triangle edges.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->triangleEdgeVector_;
}

template<size_t card>
const std::vector<std::array<ttk::SimplexId, 3>> *
  ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getTriangles)() {

  if(this->triangleList_.empty()) {
    Timer t;

    this->triangleList_.resize(triangleNumber_);
    for(ttk::SimplexId i = 0; i < triangleNumber_; ++i) {
      for(int j = 0; j < 3; ++j)
        this->getTriangleVertexInternal(i, j, this->triangleList_[i][j]);
    }

    this->printMsg("Built " + std::to_string(triangleNumber_) + " triangles.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->triangleList_;
}

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  getTriangleLink)(const ttk::SimplexId &triangleId,
                   const int &localLinkId,
                   ttk::SimplexId &linkId) const {
  //using TrianglePosition = typename ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::TrianglePosition;
#ifndef TTK_ENABLE_KAMIKAZE
  if(localLinkId < 0 or localLinkId >= getTriangleLinkNumber(triangleId))
    return -1;
#endif

  linkId = -1;
  const auto &p = this->underlying().getTriangleCoords(triangleId);

  switch(this->underlying().getTrianglePosition(triangleId)) {
    case Derived::TrianglePosition::F_3D:
      linkId = this->getTriangleLinkF(p.data(), localLinkId);
      break;
    case Derived::TrianglePosition::H_3D:
      linkId = this->getTriangleLinkH(p.data(), localLinkId);
      break;
    case Derived::TrianglePosition::C_3D:
      linkId = this->getTriangleLinkC(p.data(), localLinkId);
      break;
    case Derived::TrianglePosition::D1_3D:
      linkId = this->getTriangleLinkD1(p.data(), localLinkId);
      break;
    case Derived::TrianglePosition::D2_3D:
      linkId = this->getTriangleLinkD2(p.data(), localLinkId);
      break;
    case Derived::TrianglePosition::D3_3D:
      linkId = this->getTriangleLinkD3(p.data(), localLinkId);
      break;
    default:
      break;
  }

  return 0;
}

template<size_t card>
ttk::SimplexId ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(
  getTriangleLinkNumber)(const ttk::SimplexId &triangleId) const {

  return this->getTriangleStarNumber(triangleId);
}

template<size_t card>
const std::vector<std::vector<ttk::SimplexId>> *
  ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(
    getTriangleLinks)() {

  if(this->triangleLinkList_.empty()) {
    Timer t;

    this->triangleLinkList_.resize(triangleNumber_);
    for(ttk::SimplexId i = 0; i < triangleNumber_; ++i) {
      this->triangleLinkList_[i].resize(getTriangleLinkNumber(i));
      for(ttk::SimplexId j = 0; j < (ttk::SimplexId)this->triangleLinkList_[i].size(); ++j)
        this->getTriangleLink(i, j, this->triangleLinkList_[i][j]);
    }

    this->printMsg("Built " + std::to_string(triangleNumber_) + " triangle links.", 1,
             t.getElapsedTime(), 1);
  }
  return &this->triangleLinkList_;
}

template <size_t card, typename Derived>
ttk::SimplexId
  ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
    getTriangleStarNumber)(const ttk::SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
#else
  TTK_FORCE_USE(triangleId);
#endif

  if(this->dimensionality_ == 3) {
    return 2;
  }
  return 0;
}

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  getTriangleStar)(const ttk::SimplexId &triangleId,
                   const int &localStarId,
                   ttk::SimplexId &starId) const {
  using CRTPClass = ttk::PeriodicImplicitTriangulationCRTP<card, Derived>;
  //using TrianglePosition = typename ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::TrianglePosition;
#ifndef TTK_ENABLE_KAMIKAZE
  if(localStarId < 0 or localStarId >= getTriangleStarNumber(triangleId))
    return -1;
#endif

  starId = -1;
  const auto &p = this->underlying().getTriangleCoords(triangleId);

  switch(this->underlying().getTrianglePosition(triangleId)) {
    case CRTPClass::TrianglePosition::F_3D:
      starId = this->getTriangleStarF(p.data(), localStarId);
      break;
    case CRTPClass::TrianglePosition::H_3D:
      starId = this->getTriangleStarH(p.data(), localStarId);
      break;
    case CRTPClass::TrianglePosition::C_3D:
      starId = this->getTriangleStarC(p.data(), localStarId);
      break;
    case CRTPClass::TrianglePosition::D1_3D:
      starId = this->getTriangleStarD1(p.data(), localStarId);
      break;
    case CRTPClass::TrianglePosition::D2_3D:
      starId = this->getTriangleStarD2(p.data(), localStarId);
      break;
    case CRTPClass::TrianglePosition::D3_3D:
      starId = this->getTriangleStarD3(p.data(), localStarId);
      break;
    default:
      break;
  }

  return 0;
}

template<size_t card>
const std::vector<std::vector<ttk::SimplexId>> *
  ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(
    getTriangleStars)() {

  if(this->triangleStarList_.empty()) {
    Timer t;

    this->triangleStarList_.resize(triangleNumber_);
    for(ttk::SimplexId i = 0; i < triangleNumber_; ++i) {
      this->triangleStarList_[i].resize(this->getTriangleStarNumber(i));
      for(ttk::SimplexId j = 0; j < (ttk::SimplexId)this->triangleStarList_[i].size(); ++j)
        this->getTriangleStar(i, j, this->triangleStarList_[i][j]);
    }

    this->printMsg("Built " + std::to_string(triangleNumber_) + " triangle stars.", 1,
             t.getElapsedTime(), 1);
  }
  return &this->triangleStarList_;
}

template<size_t card>
ttk::SimplexId ttk::PeriodicImplicitTriangulation<card>::getTriangleNeighborNumber(
  const ttk::SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
#endif

  TTK_FORCE_USE(triangleId);
  if(this->dimensionality_ == 2) {
    return 3;
  }

  return 0;
}

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::getTriangleNeighbor(
  const ttk::SimplexId &triangleId,
  const int &localNeighborId,
  ttk::SimplexId &neighborId) const {
  using CRTPClass = ttk::PeriodicImplicitTriangulationCRTP<card, Derived>;
  //using TrianglePosition = typename ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::TrianglePosition;
#ifndef TTK_ENABLE_KAMIKAZE
  if(localNeighborId < 0
     or localNeighborId >= getTriangleNeighborNumber(triangleId))
    return -1;
#endif

  neighborId = -1;
  const auto &p = this->underlying().getTriangleCoords(triangleId);

  switch(this->underlying().getTrianglePosition(triangleId)) {
    case CRTPClass::TrianglePosition::BOTTOM_2D:

      if(p[0] / 2 == this->nbvoxels_[this->Di_] and p[1] == this->nbvoxels_[this->Dj_]) {
        if(localNeighborId == 0) {
          neighborId = triangleId - 1;
        } else if(localNeighborId == 1) {
          neighborId = triangleId + 1 - this->wrap_[0] * 2;
        } else if(localNeighborId == 2) {
          neighborId = triangleId + this->tshift_[0] - 1 - this->wrap_[1] * 2;
        }

      } else if(p[0] / 2 == this->nbvoxels_[this->Di_]) {
        if(localNeighborId == 0) {
          neighborId = triangleId - 1;
        } else if(localNeighborId == 1) {
          neighborId = triangleId + 1 - this->wrap_[0] * 2;
        } else if(localNeighborId == 2) {
          neighborId = triangleId + this->tshift_[0] - 1;
        }

      } else if(p[1] == this->nbvoxels_[this->Dj_]) {
        if(localNeighborId == 0) {
          neighborId = triangleId - 1;
        } else if(localNeighborId == 1) {
          neighborId = triangleId + 1;
        } else if(localNeighborId == 2) {
          neighborId = triangleId + this->tshift_[0] - 1 - this->wrap_[1] * 2;
        }

      } else {
        if(localNeighborId == 0) {
          neighborId = triangleId - 1;
        } else if(localNeighborId == 1) {
          neighborId = triangleId + 1;
        } else if(localNeighborId == 2) {
          neighborId = triangleId + this->tshift_[0] - 1;
        }
      }
      break;

    case CRTPClass::TrianglePosition::TOP_2D:

      if(p[0] / 2 == 0 and p[1] == 0) {
        if(localNeighborId == 0) {
          neighborId = triangleId + 1;
        } else if(localNeighborId == 1) {
          neighborId = triangleId - 1 + this->wrap_[0] * 2;
        } else if(localNeighborId == 2) {
          neighborId = triangleId - this->tshift_[0] + 1 + this->wrap_[1] * 2;
        }

      } else if(p[0] / 2 == 0) {
        if(localNeighborId == 0) {
          neighborId = triangleId + 1;
        } else if(localNeighborId == 1) {
          neighborId = triangleId - 1 + this->wrap_[0] * 2;
        } else if(localNeighborId == 2) {
          neighborId = triangleId - this->tshift_[0] + 1;
        }

      } else if(p[1] == 0) {
        if(localNeighborId == 0) {
          neighborId = triangleId + 1;
        } else if(localNeighborId == 1) {
          neighborId = triangleId - 1;
        } else if(localNeighborId == 2) {
          neighborId = triangleId - this->tshift_[0] + 1 + this->wrap_[1] * 2;
        }

      } else {
        if(localNeighborId == 0) {
          neighborId = triangleId + 1;
        } else if(localNeighborId == 1) {
          neighborId = triangleId - 1;
        } else if(localNeighborId == 2) {
          neighborId = triangleId - this->tshift_[0] + 1;
        }
      }
      break;
    default:
      break;
  }

  return 0;
}

template <size_t card>
int ttk::PeriodicImplicitTriangulation<card>::getTriangleNeighbors(
  std::vector<std::vector<ttk::SimplexId>> &neighbors) {
  neighbors.resize(triangleNumber_);
  for(ttk::SimplexId i = 0; i < triangleNumber_; ++i) {
    neighbors[i].resize(getTriangleNeighborNumber(i));
    for(ttk::SimplexId j = 0; j < (ttk::SimplexId)neighbors[i].size(); ++j)
      getTriangleNeighbor(i, j, neighbors[i][j]);
  }
  return 0;
}

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::getTetrahedronVertex(
  const ttk::SimplexId &tetId, const int &localVertexId, ttk::SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(tetId < 0 or tetId >= tetrahedronNumber_)
    return -1;
  if(localVertexId < 0 or localVertexId >= 4)
    return -2;
#endif

  vertexId = -1;

  if(this->dimensionality_ == 3) {
    const auto &p = this->underlying().getTetrahedronCoords(tetId);
    const ttk::SimplexId id = tetId % 6;

    switch(id) {
      case 0:
        vertexId = this->getTetrahedronVertexABCG(p.data(), localVertexId);
        break;
      case 1:
        vertexId = this->getTetrahedronVertexBCDG(p.data(), localVertexId);
        break;
      case 2:
        vertexId = this->getTetrahedronVertexABEG(p.data(), localVertexId);
        break;
      case 3:
        vertexId = this->getTetrahedronVertexBEFG(p.data(), localVertexId);
        break;
      case 4:
        vertexId = this->getTetrahedronVertexBFGH(p.data(), localVertexId);
        break;
      case 5:
        vertexId = this->getTetrahedronVertexBDGH(p.data(), localVertexId);
        break;
    }
  }
  return 0;
}

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::getTetrahedronEdge(
  const ttk::SimplexId &tetId, const int &localEdgeId, ttk::SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(tetId < 0 or tetId >= tetrahedronNumber_)
    return -1;
  if(localEdgeId < 0 or localEdgeId >= 6)
    return -2;
#endif

  edgeId = -1;

  if(this->dimensionality_ == 3) {
    const auto &p = this->underlying().getTetrahedronCoords(tetId);
    const ttk::SimplexId id = tetId % 6;

    switch(id) {
      case 0:
        edgeId = this->getTetrahedronEdgeABCG(p.data(), localEdgeId);
        break;
      case 1:
        edgeId = this->getTetrahedronEdgeBCDG(p.data(), localEdgeId);
        break;
      case 2:
        edgeId = this->getTetrahedronEdgeABEG(p.data(), localEdgeId);
        break;
      case 3:
        edgeId = this->getTetrahedronEdgeBEFG(p.data(), localEdgeId);
        break;
      case 4:
        edgeId = this->getTetrahedronEdgeBFGH(p.data(), localEdgeId);
        break;
      case 5:
        edgeId = this->getTetrahedronEdgeBDGH(p.data(), localEdgeId);
        break;
    }
  }

  return 0;
}

template<size_t card>
int ttk::PeriodicImplicitTriangulation<card>::getTetrahedronEdges(
  std::vector<std::vector<ttk::SimplexId>> &edges) const {
  edges.resize(tetrahedronNumber_);
  for(ttk::SimplexId i = 0; i < tetrahedronNumber_; ++i) {
    edges[i].resize(6);
    for(int j = 0; j < 6; ++j)
      this->getTetrahedronEdge(i, j, edges[i][j]);
  }

  return 0;
}

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::getTetrahedronTriangle(
  const ttk::SimplexId &tetId,
  const int &localTriangleId,
  ttk::SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(tetId < 0 or tetId >= tetrahedronNumber_)
    return -1;
  if(localTriangleId < 0 or localTriangleId >= 4)
    return -2;
#endif

  triangleId = -1;

  if(this->dimensionality_ == 3) {
    const auto &p = this->underlying().getTetrahedronCoords(tetId);
    const ttk::SimplexId id = tetId % 6;

    switch(id) {
      case 0:
        triangleId = this->getTetrahedronTriangleABCG(p.data(), localTriangleId);
        break;
      case 1:
        triangleId = this->getTetrahedronTriangleBCDG(p.data(), localTriangleId);
        break;
      case 2:
        triangleId = this->getTetrahedronTriangleABEG(p.data(), localTriangleId);
        break;
      case 3:
        triangleId = this->getTetrahedronTriangleBEFG(p.data(), localTriangleId);
        break;
      case 4:
        triangleId = this->getTetrahedronTriangleBFGH(p.data(), localTriangleId);
        break;
      case 5:
        triangleId = this->getTetrahedronTriangleBDGH(p.data(), localTriangleId);
        break;
    }
  }

  return 0;
}

template<size_t card>
int ttk::PeriodicImplicitTriangulation<card>::getTetrahedronTriangles(
  std::vector<std::vector<ttk::SimplexId>> &triangles) const {
  triangles.resize(tetrahedronNumber_);
  for(ttk::SimplexId i = 0; i < tetrahedronNumber_; ++i) {
    triangles[i].resize(4);
    for(int j = 0; j < 4; ++j)
      this->getTetrahedronTriangle(i, j, triangles[i][j]);
  }

  return 0;
}

template<size_t card>
ttk::SimplexId ttk::PeriodicImplicitTriangulation<card>::getTetrahedronNeighborNumber(
  const ttk::SimplexId &tetId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(tetId < 0 or tetId >= tetrahedronNumber_)
    return -1;
#endif

  TTK_FORCE_USE(tetId);
  if(this->dimensionality_ == 3) {
    return 4;
  }

  return 0;
}

template <size_t card, typename Derived>
int ttk::PeriodicImplicitTriangulationCRTP<card, Derived>::getTetrahedronNeighbor(
  const ttk::SimplexId &tetId,
  const int &localNeighborId,
  ttk::SimplexId &neighborId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localNeighborId < 0
     or localNeighborId >= this->getTetrahedronNeighborNumber(tetId))
    return -1;
#endif

  neighborId = -1;

  if(this->dimensionality_ == 3) {
    const auto &p = this->underlying().getTetrahedronCoords(tetId);
    const ttk::SimplexId id = tetId % 6;

    switch(id) {
      case 0:
        neighborId
          = this->getTetrahedronNeighborABCG(tetId, p.data(), localNeighborId);
        break;
      case 1:
        neighborId
          = this->getTetrahedronNeighborBCDG(tetId, p.data(), localNeighborId);
        break;
      case 2:
        neighborId
          = this->getTetrahedronNeighborABEG(tetId, p.data(), localNeighborId);
        break;
      case 3:
        neighborId
          = this->getTetrahedronNeighborBEFG(tetId, p.data(), localNeighborId);
        break;
      case 4:
        neighborId
          = this->getTetrahedronNeighborBFGH(tetId, p.data(), localNeighborId);
        break;
      case 5:
        neighborId
          = this->getTetrahedronNeighborBDGH(tetId, p.data(), localNeighborId);
        break;
    }
  }

  return 0;
}

template<size_t card>
int ttk::PeriodicImplicitTriangulation<card>::getTetrahedronNeighbors(
  std::vector<std::vector<ttk::SimplexId>> &neighbors) {
  neighbors.resize(tetrahedronNumber_);
  for(ttk::SimplexId i = 0; i < tetrahedronNumber_; ++i) {
    neighbors[i].resize(getTetrahedronNeighborNumber(i));
    for(ttk::SimplexId j = 0; j < (ttk::SimplexId)neighbors[i].size(); ++j)
      this->getTetrahedronNeighbor(i, j, neighbors[i][j]);
  }

  return 0;
}

template<size_t card>
ttk::SimplexId ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(
  getCellVertexNumber)(const ttk::SimplexId &ttkNotUsed(cellId)) const {

  return this->dimensionality_ + 1;
}

template<size_t card>
int ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getCellVertex)(
  const ttk::SimplexId &cellId,
  const int &localVertexId,
  ttk::SimplexId &vertexId) const {

  if(this->dimensionality_ == 3)
    this->getTetrahedronVertex(cellId, localVertexId, vertexId);
  else if(this->dimensionality_ == 2)
    this->getTriangleVertexInternal(cellId, localVertexId, vertexId);
  else if(this->dimensionality_ == 1)
    this->getEdgeVertexInternal(cellId, localVertexId, vertexId);

  return 0;
}

template<size_t card>
ttk::SimplexId ttk::PeriodicImplicitTriangulation<card>::getCellEdgeNumberInternal(
  const ttk::SimplexId &ttkNotUsed(cellId)) const {
  if(this->dimensionality_ == 3)
    return 6;
  else if(this->dimensionality_ == 2)
    return 3;

  return 0;
}

template<size_t card>
int ttk::PeriodicImplicitTriangulation<card>::getCellEdgeInternal(
  const ttk::SimplexId &cellId, const int &localEdgeId, ttk::SimplexId &edgeId) const {
  if(this->dimensionality_ == 3)
    this->getTetrahedronEdge(cellId, localEdgeId, edgeId);
  else if(this->dimensionality_ == 2)
    this->getTriangleEdgeInternal(cellId, localEdgeId, edgeId);
  else if(this->dimensionality_ == 1)
    getCellNeighbor(cellId, localEdgeId, edgeId);

  return 0;
}

template<size_t card>
const std::vector<std::vector<ttk::SimplexId>> *
  ttk::PeriodicImplicitTriangulation<card>::getCellEdgesInternal() {
  if(this->cellEdgeVector_.empty()) {
    Timer t;

    if(this->dimensionality_ == 3)
      this->getTetrahedronEdges(this->cellEdgeVector_);
    else if(this->dimensionality_ == 2)
      getTriangleEdgesInternal(this->cellEdgeVector_);

    this->printMsg("Built " + std::to_string(cellNumber_) + " cell edges.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->cellEdgeVector_;
}

template<size_t card>
int ttk::PeriodicImplicitTriangulation<card>::getCellTriangleInternal(
  const ttk::SimplexId &cellId,
  const int &localTriangleId,
  ttk::SimplexId &triangleId) const {
  if(this->dimensionality_ == 3)
    this->getTetrahedronTriangle(cellId, localTriangleId, triangleId);

  return 0;
}

template<size_t card>
const std::vector<std::vector<ttk::SimplexId>> *
  ttk::PeriodicImplicitTriangulation<card>::getCellTrianglesInternal() {
  if(this->cellTriangleVector_.empty()) {
    Timer t;

    if(this->dimensionality_ == 3)
      this->getTetrahedronTriangles(this->cellTriangleVector_);

    this->printMsg("Built " + std::to_string(cellNumber_) + " cell triangles.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->cellTriangleVector_;
}

template<size_t card>
ttk::SimplexId ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(
  getCellNeighborNumber)(const ttk::SimplexId &cellId) const {

  if(this->dimensionality_ == 3)
    return this->getTetrahedronNeighborNumber(cellId);
  else if(this->dimensionality_ == 2)
    return getTriangleNeighborNumber(cellId);
  else if(this->dimensionality_ == 1) {
    this->printErr("getCellNeighborNumber() not implemented in 1D! (TODO)");
    return -1;
  }

  return 0;
}

template<size_t card>
int ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getCellNeighbor)(
  const ttk::SimplexId &cellId,
  const int &localNeighborId,
  ttk::SimplexId &neighborId) const {

  if(this->dimensionality_ == 3)
    this->getTetrahedronNeighbor(cellId, localNeighborId, neighborId);
  else if(this->dimensionality_ == 2)
    getTriangleNeighbor(cellId, localNeighborId, neighborId);
  else if(this->dimensionality_ == 1) {
    this->printErr("getCellNeighbor() not implemented in 1D! (TODO)");
    return -1;
  }

  return 0;
}

template<size_t card>
const std::vector<std::vector<ttk::SimplexId>> *
  ttk::PeriodicImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(
    getCellNeighbors)() {

  if(this->cellNeighborList_.empty()) {
    Timer t;

    if(this->dimensionality_ == 3)
      this->getTetrahedronNeighbors(this->cellNeighborList_);
    else if(this->dimensionality_ == 2)
      getTriangleNeighbors(this->cellNeighborList_);
    else if(this->dimensionality_ == 1) {
      this->printErr("getCellNeighbors() not implemented in 1D! (TODO)");
      return nullptr;
    }

    this->printMsg("Built " + std::to_string(cellNumber_) + " cell neighbors.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->cellNeighborList_;
}

// explicit instantiations

template class ttk::PeriodicImplicitTriangulation<0>;
template class ttk::PeriodicImplicitTriangulation<1>;
template class ttk::PeriodicImplicitTriangulation<2>;
template class ttk::PeriodicImplicitTriangulation<3>;
/*
*/
template class ttk::PeriodicImplicitTriangulationCRTP<0, ttk::PeriodicWithPreconditions<0>>;
template class ttk::PeriodicImplicitTriangulationCRTP<0, ttk::PeriodicNoPreconditions<0>>;
template class ttk::PeriodicImplicitTriangulationCRTP<1, ttk::PeriodicWithPreconditions<1>>;
template class ttk::PeriodicImplicitTriangulationCRTP<1, ttk::PeriodicNoPreconditions<1>>;
template class ttk::PeriodicImplicitTriangulationCRTP<2, ttk::PeriodicWithPreconditions<2>>;
template class ttk::PeriodicImplicitTriangulationCRTP<2, ttk::PeriodicNoPreconditions<2>>;
template class ttk::PeriodicImplicitTriangulationCRTP<3, ttk::PeriodicWithPreconditions<3>>;
template class ttk::PeriodicImplicitTriangulationCRTP<3, ttk::PeriodicNoPreconditions<3>>;

/// @endcond
