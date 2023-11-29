#include <MultiresTriangulation.h>

using ttk::MultiresTriangulation;
using ttk::SimplexId;

template <size_t card>
MultiresTriangulation<card>::MultiresTriangulation() {
  decimation_ = 1;
  decimationLevel_ = 0;
  debugLevel_ = 0;
  gridDimensions_[0] = gridDimensions_[1] = gridDimensions_[2] = -1;
  triangulation_ = nullptr;
  this->setDebugMsgPrefix("MultiresTriangulation");
}

template <size_t card>
MultiresTriangulation<card>::~MultiresTriangulation() = default;

template <size_t card>
int ttk::MultiresTriangulation<card>::preconditionVerticesInternal() {
  vertexPositions_.resize(vertexNumber_);
  vertexCoords_.resize(vertexNumber_);

  if(dimensionality_ == 1) {
    vertexPositions_[0] = VertexPosition::LEFT_CORNER_1D;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 1; i < vertexNumber_ - 1; ++i) {
      vertexPositions_[i] = VertexPosition::CENTER_1D;
    }
    vertexPositions_[vertexNumber_ - 1] = VertexPosition::RIGHT_CORNER_1D;

  } else if(dimensionality_ == 2) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      std::array<SimplexId, 3> p{};
      vertexToPosition2d(i, p);

      if(0 < p[0] and p[0] < nbvoxels_[Di_]) {
        if(0 < p[1] and p[1] < nbvoxels_[Dj_])
          vertexPositions_[i] = VertexPosition::CENTER_2D;
        else if(p[1] == 0)
          vertexPositions_[i] = VertexPosition::TOP_EDGE_2D; // ab
        else
          vertexPositions_[i] = VertexPosition::BOTTOM_EDGE_2D; // cd
      } else if(p[0] == 0) {
        if(0 < p[1] and p[1] < nbvoxels_[Dj_])
          vertexPositions_[i] = VertexPosition::LEFT_EDGE_2D; // ac
        else if(p[1] == 0)
          vertexPositions_[i] = VertexPosition::TOP_LEFT_CORNER_2D; // a
        else
          vertexPositions_[i] = VertexPosition::BOTTOM_LEFT_CORNER_2D; // c
      } else {
        if(0 < p[1] and p[1] < nbvoxels_[Dj_])
          vertexPositions_[i] = VertexPosition::RIGHT_EDGE_2D; // bd
        else if(p[1] == 0)
          vertexPositions_[i] = VertexPosition::TOP_RIGHT_CORNER_2D; // b
        else
          vertexPositions_[i] = VertexPosition::BOTTOM_RIGHT_CORNER_2D; // d
      }
      vertexCoords_[i] = p;
    }

  } else if(dimensionality_ == 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      std::array<SimplexId, 3> p{};
      vertexToPosition(i, p);

      if(0 < p[0] and p[0] < nbvoxels_[0]) {
        if(0 < p[1] and p[1] < nbvoxels_[1]) {
          if(0 < p[2] and p[2] < nbvoxels_[2])
            vertexPositions_[i] = VertexPosition::CENTER_3D;
          else if(p[2] == 0)
            vertexPositions_[i] = VertexPosition::FRONT_FACE_3D; // abcd
          else
            vertexPositions_[i] = VertexPosition::BACK_FACE_3D; // efgh
        } else if(p[1] == 0) {
          if(0 < p[2] and p[2] < nbvoxels_[2])
            vertexPositions_[i] = VertexPosition::TOP_FACE_3D; // abef
          else if(p[2] == 0)
            vertexPositions_[i] = VertexPosition::TOP_FRONT_EDGE_3D; // ab
          else
            vertexPositions_[i] = VertexPosition::TOP_BACK_EDGE_3D; // ef
        } else {
          if(0 < p[2] and p[2] < nbvoxels_[2])
            vertexPositions_[i] = VertexPosition::BOTTOM_FACE_3D; // cdgh
          else if(p[2] == 0)
            vertexPositions_[i] = VertexPosition::BOTTOM_FRONT_EDGE_3D; // cd
          else
            vertexPositions_[i] = VertexPosition::BOTTOM_BACK_EDGE_3D; // gh
        }
      } else if(p[0] == 0) {
        if(0 < p[1] and p[1] < nbvoxels_[1]) {
          if(0 < p[2] and p[2] < nbvoxels_[2])
            vertexPositions_[i] = VertexPosition::LEFT_FACE_3D; // aceg
          else if(p[2] == 0)
            vertexPositions_[i] = VertexPosition::LEFT_FRONT_EDGE_3D; // ac
          else
            vertexPositions_[i] = VertexPosition::LEFT_BACK_EDGE_3D; // eg
        } else if(p[1] == 0) {
          if(0 < p[2] and p[2] < nbvoxels_[2])
            vertexPositions_[i] = VertexPosition::TOP_LEFT_EDGE_3D; // ae
          else if(p[2] == 0)
            vertexPositions_[i] = VertexPosition::TOP_LEFT_FRONT_CORNER_3D; // a
          else
            vertexPositions_[i] = VertexPosition::TOP_LEFT_BACK_CORNER_3D; // e
        } else {
          if(0 < p[2] and p[2] < nbvoxels_[2])
            vertexPositions_[i] = VertexPosition::BOTTOM_LEFT_EDGE_3D; // cg
          else if(p[2] == 0)
            vertexPositions_[i]
              = VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D; // c
          else
            vertexPositions_[i]
              = VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D; // g
        }
      } else {
        if(0 < p[1] and p[1] < nbvoxels_[1]) {
          if(0 < p[2] and p[2] < nbvoxels_[2])
            vertexPositions_[i] = VertexPosition::RIGHT_FACE_3D; // bdfh
          else if(p[2] == 0)
            vertexPositions_[i] = VertexPosition::RIGHT_FRONT_EDGE_3D; // bd
          else
            vertexPositions_[i] = VertexPosition::RIGHT_BACK_EDGE_3D; // fh
        } else if(p[1] == 0) {
          if(0 < p[2] and p[2] < nbvoxels_[2])
            vertexPositions_[i] = VertexPosition::TOP_RIGHT_EDGE_3D; // bf
          else if(p[2] == 0)
            vertexPositions_[i]
              = VertexPosition::TOP_RIGHT_FRONT_CORNER_3D; // b
          else
            vertexPositions_[i] = VertexPosition::TOP_RIGHT_BACK_CORNER_3D; // f
        } else {
          if(0 < p[2] and p[2] < nbvoxels_[2])
            vertexPositions_[i] = VertexPosition::BOTTOM_RIGHT_EDGE_3D; // dh
          else if(p[2] == 0)
            vertexPositions_[i]
              = VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D; // d
          else
            vertexPositions_[i]
              = VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D; // h
        }
      }
      vertexCoords_[i] = p;
    }
  }
  return 0;
}

template<size_t card>
int ttk::MultiresTriangulation<card>::computeVerticesShifts() {
  vertexShifts_.resize(vertexNumber_);

  if(dimensionality_ == 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < decimatedVertexNumber_; ++i) {
      const auto v = localToGlobalVertexId(i);
      const auto &p = this->vertexCoords_[v];
      std::array<SimplexId, 3> shifts{decimation_, decimation_, decimation_};

      if((nbvoxels_[0] % decimation_) and (p[0] + decimation_ > nbvoxels_[0])) {
        shifts[0] = nbvoxels_[0] % decimation_;
      }
      if((nbvoxels_[1] % decimation_) and (p[1] + decimation_ > nbvoxels_[1])) {
        shifts[1] = nbvoxels_[1] % decimation_;
      }
      if((nbvoxels_[2] % decimation_) and (p[2] + decimation_ > nbvoxels_[2])) {
        shifts[2] = nbvoxels_[2] % decimation_;
      }
      this->vertexShifts_[v] = shifts;
    }
  } else if(dimensionality_ == 2) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < decimatedVertexNumber_; ++i) {
      const auto v = localToGlobalVertexId(i);
      const auto &p = this->vertexCoords_[v];
      std::array<SimplexId, 3> shifts{decimation_, decimation_, decimation_};

      if((nbvoxels_[Di_] % decimation_)
         and (p[0] + decimation_ > nbvoxels_[Di_])) {
        shifts[0] = nbvoxels_[0] % decimation_;
      }
      if((nbvoxels_[Dj_] % decimation_)
         and (p[1] + decimation_ > nbvoxels_[Dj_])) {
        shifts[1] = nbvoxels_[1] % decimation_;
      }
      this->vertexShifts_[v] = shifts;
    }
  }

  return 0;
}

template<size_t card>
int MultiresTriangulation<card>::getInvertVertexNeighbor(
  const SimplexId &vertexId,
  const SimplexId &neighborId,
  SimplexId &localNeighborId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  // if(triangulation_->isEmptyCheck()) {
  //   return -1;
  // }
  // if(localNeighborId < 0
  //    or localNeighborId >= getInvertVertexNeighborNumber(vertexId))
  //   return -1;
#endif

  const auto &p = this->vertexCoords_[vertexId];

  SimplexId shiftX = decimation_;
  SimplexId shiftY = decimation_;
  SimplexId shiftZ = decimation_;

  if(dimensionality_ == 3) {
    if((nbvoxels_[0] % decimation_) and (p[0] + decimation_ > nbvoxels_[0])) {
      shiftX = nbvoxels_[0] % decimation_;
    }
    if((nbvoxels_[1] % decimation_) and (p[1] + decimation_ > nbvoxels_[1])) {
      shiftY = nbvoxels_[1] % decimation_;
    }
    if((nbvoxels_[2] % decimation_) and (p[2] + decimation_ > nbvoxels_[2])) {
      shiftZ = nbvoxels_[2] % decimation_;
    }
  } else if(dimensionality_ == 2) {
    if((nbvoxels_[Di_] % decimation_)
       and (p[0] + decimation_ > nbvoxels_[Di_])) {
      shiftX = nbvoxels_[0] % decimation_;
    }
    if((nbvoxels_[Dj_] % decimation_)
       and (p[1] + decimation_ > nbvoxels_[Dj_])) {
      shiftY = nbvoxels_[1] % decimation_;
    }
  }

  switch(vertexPositions_[vertexId]) {
    case VertexPosition::CENTER_3D:
      localNeighborId = getInvertVertexNeighborABCDEFGH(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::FRONT_FACE_3D:
      localNeighborId = getInvertVertexNeighborABDC(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::BACK_FACE_3D:
      localNeighborId = getInvertVertexNeighborEFHG(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::TOP_FACE_3D:
      localNeighborId = getInvertVertexNeighborAEFB(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::TOP_FRONT_EDGE_3D:
      localNeighborId = getInvertVertexNeighborAB(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::TOP_BACK_EDGE_3D:
      localNeighborId = getInvertVertexNeighborEF(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::BOTTOM_FACE_3D:
      localNeighborId = getInvertVertexNeighborGHDC(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
      localNeighborId = getInvertVertexNeighborCD(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
      localNeighborId = getInvertVertexNeighborGH(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::LEFT_FACE_3D:
      localNeighborId = getInvertVertexNeighborAEGC(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
      localNeighborId = getInvertVertexNeighborAC(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::LEFT_BACK_EDGE_3D: // eg
      localNeighborId = getInvertVertexNeighborEG(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::TOP_LEFT_EDGE_3D: // ae
      localNeighborId = getInvertVertexNeighborAE(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
      localNeighborId = getInvertVertexNeighborA(
        vertexId, neighborId, shiftX, shiftY, shiftZ); // e
      break;
    case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
      localNeighborId = getInvertVertexNeighborE(
        vertexId, neighborId, shiftX, shiftY, shiftZ); // a
      break;
    case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
      localNeighborId = getInvertVertexNeighborCG(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
      localNeighborId = getInvertVertexNeighborC(
        vertexId, neighborId, shiftX, shiftY, shiftZ); // g
      break;
    case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
      localNeighborId = getInvertVertexNeighborG(
        vertexId, neighborId, shiftX, shiftY, shiftZ); // c
      break;
    case VertexPosition::RIGHT_FACE_3D:
      localNeighborId = getInvertVertexNeighborBFHD(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
      localNeighborId = getInvertVertexNeighborBD(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
      localNeighborId = getInvertVertexNeighborFH(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
      localNeighborId = getInvertVertexNeighborBF(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
      localNeighborId = getInvertVertexNeighborB(
        vertexId, neighborId, shiftX, shiftY, shiftZ); // f
      break;
    case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
      localNeighborId = getInvertVertexNeighborF(
        vertexId, neighborId, shiftX, shiftY, shiftZ); // b
      break;
    case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
      localNeighborId = getInvertVertexNeighborDH(
        vertexId, neighborId, shiftX, shiftY, shiftZ); //
      break;
    case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
      localNeighborId = getInvertVertexNeighborD(
        vertexId, neighborId, shiftX, shiftY, shiftZ); // h
      break;
    case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
      localNeighborId = getInvertVertexNeighborH(
        vertexId, neighborId, shiftX, shiftY, shiftZ); // d
      break;
    case VertexPosition::CENTER_2D:
      localNeighborId = getInvertVertexNeighbor2dABCD(
        vertexId, neighborId, shiftX, shiftY); // abcd
      break;
    case VertexPosition::TOP_EDGE_2D:
      localNeighborId = getInvertVertexNeighbor2dAB(
        vertexId, neighborId, shiftX, shiftY); // ab
      break;
    case VertexPosition::BOTTOM_EDGE_2D:
      localNeighborId = getInvertVertexNeighbor2dCD(
        vertexId, neighborId, shiftX, shiftY); // cd
      break;
    case VertexPosition::LEFT_EDGE_2D:
      localNeighborId = getInvertVertexNeighbor2dAC(
        vertexId, neighborId, shiftX, shiftY); // ac
      break;
    case VertexPosition::TOP_LEFT_CORNER_2D: // a
      localNeighborId
        = getInvertVertexNeighbor2dA(vertexId, neighborId, shiftX, shiftY); // a
      break;
    case VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
      localNeighborId
        = getInvertVertexNeighbor2dC(vertexId, neighborId, shiftX, shiftY); // c
      break;
    case VertexPosition::RIGHT_EDGE_2D:
      localNeighborId = getInvertVertexNeighbor2dBD(
        vertexId, neighborId, shiftX, shiftY); // bd
      break;
    case VertexPosition::TOP_RIGHT_CORNER_2D: // b
      localNeighborId
        = getInvertVertexNeighbor2dB(vertexId, neighborId, shiftX, shiftY); // b
      break;
    case VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
      localNeighborId
        = getInvertVertexNeighbor2dD(vertexId, neighborId, shiftX, shiftY); // d
      break;
    case VertexPosition::CENTER_1D:
      localNeighborId = localNeighborId == 0 ? vertexId + decimation_
                                             : vertexId - decimation_;
      break;
    case VertexPosition::LEFT_CORNER_1D:
      localNeighborId = vertexId + decimation_;
      break;
    case VertexPosition::RIGHT_CORNER_1D:
      localNeighborId = vertexId - decimation_;
      break;
    default:
      localNeighborId = -1;
      break;
  }

  return 0;
}

template<size_t card>
int MultiresTriangulation<card>::getVertexNeighbor(const SimplexId &vertexId,
                                             const int &localNeighborId,
                                             SimplexId &neighborId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  // if(triangulation_->isEmptyCheck()) {
  //   return -1;
  // }
  // if(localNeighborId < 0
  //    or localNeighborId >= getVertexNeighborNumber(vertexId))
  //   return -1;
#endif

  const auto &s = vertexShifts_[vertexId];

  switch(vertexPositions_[vertexId]) {
    case VertexPosition::CENTER_3D:
      neighborId = getVertexNeighborABCDEFGH(
        vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::FRONT_FACE_3D:
      neighborId
        = getVertexNeighborABDC(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::BACK_FACE_3D:
      neighborId
        = getVertexNeighborEFHG(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::TOP_FACE_3D:
      neighborId
        = getVertexNeighborAEFB(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::BOTTOM_FACE_3D:
      neighborId
        = getVertexNeighborGHDC(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::LEFT_FACE_3D:
      neighborId
        = getVertexNeighborAEGC(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::RIGHT_FACE_3D:
      neighborId
        = getVertexNeighborBFHD(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::TOP_FRONT_EDGE_3D: // ab
      neighborId
        = getVertexNeighborAB(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
      neighborId
        = getVertexNeighborCD(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
      neighborId
        = getVertexNeighborAC(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
      neighborId
        = getVertexNeighborBD(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::TOP_BACK_EDGE_3D: // ef
      neighborId
        = getVertexNeighborEF(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
      neighborId
        = getVertexNeighborGH(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::LEFT_BACK_EDGE_3D: // eg
      neighborId
        = getVertexNeighborEG(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
      neighborId
        = getVertexNeighborFH(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::TOP_LEFT_EDGE_3D: // ae
      neighborId
        = getVertexNeighborAE(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
      neighborId
        = getVertexNeighborBF(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
      neighborId
        = getVertexNeighborCG(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
      neighborId
        = getVertexNeighborDH(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
      neighborId
        = getVertexNeighborA(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
      neighborId
        = getVertexNeighborB(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
      neighborId
        = getVertexNeighborC(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
      neighborId
        = getVertexNeighborD(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
      neighborId
        = getVertexNeighborE(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
      neighborId
        = getVertexNeighborF(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
      neighborId
        = getVertexNeighborG(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
      neighborId
        = getVertexNeighborH(vertexId, localNeighborId, s[0], s[1], s[2]);
      break;
    case VertexPosition::CENTER_2D:
      neighborId
        = getVertexNeighbor2dABCD(vertexId, localNeighborId, s[0], s[1]);
      break;
    case VertexPosition::TOP_EDGE_2D:
      neighborId = getVertexNeighbor2dAB(vertexId, localNeighborId, s[0], s[1]);
      break;
    case VertexPosition::BOTTOM_EDGE_2D:
      neighborId = getVertexNeighbor2dCD(vertexId, localNeighborId, s[0], s[1]);
      break;
    case VertexPosition::LEFT_EDGE_2D:
      neighborId = getVertexNeighbor2dAC(vertexId, localNeighborId, s[0], s[1]);
      break;
    case VertexPosition::RIGHT_EDGE_2D:
      neighborId = getVertexNeighbor2dBD(vertexId, localNeighborId, s[0], s[1]);
      break;
    case VertexPosition::TOP_LEFT_CORNER_2D: // a
      neighborId = getVertexNeighbor2dA(vertexId, localNeighborId, s[0], s[1]);
      break;
    case VertexPosition::TOP_RIGHT_CORNER_2D: // b
      neighborId = getVertexNeighbor2dB(vertexId, localNeighborId, s[0], s[1]);
      break;
    case VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
      neighborId = getVertexNeighbor2dC(vertexId, localNeighborId, s[0], s[1]);
      break;
    case VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
      neighborId = getVertexNeighbor2dD(vertexId, localNeighborId, s[0], s[1]);
      break;
    case VertexPosition::CENTER_1D:
      neighborId = (localNeighborId == 0 ? vertexId + decimation_
                                         : vertexId - decimation_);
      break;
    case VertexPosition::LEFT_CORNER_1D:
      neighborId = vertexId + decimation_;
      break;
    case VertexPosition::RIGHT_CORNER_1D:
      neighborId = vertexId - decimation_;
      break;
    default:
      neighborId = -1;
      break;
  }

  return 0;
}

template<size_t card>
SimplexId MultiresTriangulation<card>::getInvertVertexNeighbor2dA(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY) const {
  // V(a)={b,c}
  if(neighborId == v + shiftX)
    return 0; // b
  else if(neighborId == v + gridDimensions_[Di_] * shiftY)
    return 1; // c
  return -1;
}

template<size_t card>
SimplexId
  MultiresTriangulation<card>::getVertexNeighbor2dA(const SimplexId v,
                                              const int id,
                                              const SimplexId shiftX,
                                              const SimplexId shiftY) const {
  // V(a)={b,c}
  switch(id) {
    case 0:
      return v + shiftX; // b
    case 1:
      return v + gridDimensions_[Di_] * shiftY; // c
  }
  return -1;
}

template<size_t card>
SimplexId MultiresTriangulation<card>::getInvertVertexNeighbor2dB(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY) const {
  // V(b)={a,c,d}
  if(neighborId == v - shiftX)
    return 0; // a
  else if(neighborId == v + gridDimensions_[Di_] * shiftY)
    return 1; // d
  else if(neighborId == v + (gridDimensions_[Di_] * shiftY - shiftX))
    return 2; // c
  return -1;
}

template<size_t card>
SimplexId
  MultiresTriangulation<card>::getVertexNeighbor2dB(const SimplexId v,
                                              const int id,
                                              const SimplexId shiftX,
                                              const SimplexId shiftY) const {
  // V(b)={a,c,d}
  switch(id) {
    case 0:
      return v - shiftX; // a
    case 1:
      return v + gridDimensions_[Di_] * shiftY; // d
    case 2:
      return v + (gridDimensions_[Di_] * shiftY - shiftX); // c
  }
  return -1;
}

template<size_t card>
SimplexId MultiresTriangulation<card>::getInvertVertexNeighbor2dC(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY) const {
  // V(c)={a,b,d}
  if(neighborId == v + shiftX)
    return 0;
  else if(neighborId == v - gridDimensions_[Di_] * shiftY)
    return 1; // a
  else if(neighborId == v + (shiftX - gridDimensions_[Di_] * shiftY))
    return 2; // b
  return -1;
}

template<size_t card>
SimplexId
  MultiresTriangulation<card>::getVertexNeighbor2dC(const SimplexId v,
                                              const int id,
                                              const SimplexId shiftX,
                                              const SimplexId shiftY) const {
  // V(c)={a,b,d}
  switch(id) {
    case 0:
      return v + shiftX;
    case 1:
      return v - gridDimensions_[Di_] * shiftY; // a
    case 2:
      return v + (shiftX - gridDimensions_[Di_] * shiftY); // b
  }
  return -1;
}

template<size_t card>
SimplexId MultiresTriangulation<card>::getInvertVertexNeighbor2dD(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY) const {
  // V(d)={c,b}
  if(neighborId == v - shiftX)
    return 0; // c
  else if(neighborId == v - gridDimensions_[Di_] * shiftY)
    return 1; // b
  return -1;
}

template<size_t card>
SimplexId
  MultiresTriangulation<card>::getVertexNeighbor2dD(const SimplexId v,
                                              const int id,
                                              const SimplexId shiftX,
                                              const SimplexId shiftY) const {
  // V(d)={c,b}
  switch(id) {
    case 0:
      return v - shiftX; // c
    case 1:
      return v - gridDimensions_[Di_] * shiftY; // b
  }
  return -1;
}

template<size_t card>
SimplexId MultiresTriangulation<card>::getInvertVertexNeighbor2dAB(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY) const {
  // V(ab)=V(b)::{a,c,d}+V(a)::{b}
  if(neighborId == v - decimation_)
    return 0; // V(b)::a
  else if(neighborId == v + (gridDimensions_[Di_] * shiftY - decimation_))
    return 1; // V(b)::c
  else if(neighborId == v + gridDimensions_[Di_] * shiftY)
    return 2; // V(b)::d
  else if(neighborId == v + shiftX)
    return 3; // V(a)::b
  return -1;
}

template<size_t card>
SimplexId
  MultiresTriangulation<card>::getVertexNeighbor2dAB(const SimplexId v,
                                               const int id,
                                               const SimplexId shiftX,
                                               const SimplexId shiftY) const {
  // V(ab)=V(b)::{a,c,d}+V(a)::{b}
  switch(id) {
    case 0:
      return v - decimation_; // V(b)::a
    case 1:
      return v + (gridDimensions_[Di_] * shiftY - decimation_); // V(b)::c
    case 2:
      return v + gridDimensions_[Di_] * shiftY; // V(b)::d
    case 3:
      return v + shiftX; // V(a)::b
  }
  return -1;
}

template<size_t card>
SimplexId MultiresTriangulation<card>::getInvertVertexNeighbor2dCD(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY) const {
  // V(cd)=V(c)::{a,b,d}+V(d)::{c}
  if(neighborId == v - decimation_)
    return 0; // V(d)::c
  else if(neighborId == v - gridDimensions_[Di_] * shiftY)
    return 1; // V(c)::a
  else if(neighborId == v + (shiftX - gridDimensions_[Di_] * shiftY))
    return 2; // V(c)::b
  else if(neighborId == v + shiftX)
    return 3; // V(c)::d
  return -1;
}

template<size_t card>
SimplexId
  MultiresTriangulation<card>::getVertexNeighbor2dCD(const SimplexId v,
                                               const int id,
                                               const SimplexId shiftX,
                                               const SimplexId shiftY) const {
  // V(cd)=V(c)::{a,b,d}+V(d)::{c}
  switch(id) {
    case 0:
      return v - decimation_; // V(d)::c
    case 1:
      return v - gridDimensions_[Di_] * shiftY; // V(c)::a
    case 2:
      return v + (shiftX - gridDimensions_[Di_] * shiftY); // V(c)::b
    case 3:
      return v + shiftX; // V(c)::d
  }
  return -1;
}

template<size_t card>
SimplexId MultiresTriangulation<card>::getInvertVertexNeighbor2dAC(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY) const {
  // V(ac)=V(c)::{a,b,d}+V(a)::{c}
  if(neighborId == v - gridDimensions_[Di_] * decimation_)
    return 0; // V(c)::{a}
  else if(neighborId == v + (shiftX - gridDimensions_[Di_] * decimation_))
    return 1; // V(c)::{b}
  else if(neighborId == v + shiftX)
    return 2; // V(c)::{d}
  else if(neighborId == v + gridDimensions_[Di_] * shiftY)
    return 3; // V(a)::{c}
  return -1;
}

template<size_t card>
SimplexId
  MultiresTriangulation<card>::getVertexNeighbor2dAC(const SimplexId v,
                                               const int id,
                                               const SimplexId shiftX,
                                               const SimplexId shiftY) const {
  // V(ac)=V(c)::{a,b,d}+V(a)::{c}
  switch(id) {
    case 0:
      return v - gridDimensions_[Di_] * decimation_; // V(c)::{a}
    case 1:
      return v + (shiftX - gridDimensions_[Di_] * decimation_); // V(c)::{b}
    case 2:
      return v + shiftX; // V(c)::{d}
    case 3:
      return v + gridDimensions_[Di_] * shiftY; // V(a)::{c}
  }
  return -1;
}

template<size_t card>
SimplexId MultiresTriangulation<card>::getInvertVertexNeighbor2dBD(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY) const {
  // V(bd)=V(b)::{c,d}+V(d)::{b,c}
  if(neighborId == v + (gridDimensions_[Di_] * shiftY - shiftX))
    return 0; // V(b)::{c}
  else if(neighborId == v + gridDimensions_[Di_] * shiftY)
    return 1; // V(b)::{d}
  else if(neighborId == v - gridDimensions_[Di_] * decimation_)
    return 2; // V(d)::{b}
  else if(neighborId == v - shiftX)
    return 3; // V(d)::{c}
  return -1;
}

template<size_t card>
SimplexId
  MultiresTriangulation<card>::getVertexNeighbor2dBD(const SimplexId v,
                                               const int id,
                                               const SimplexId shiftX,
                                               const SimplexId shiftY) const {
  // V(bd)=V(b)::{c,d}+V(d)::{b,c}
  switch(id) {
    case 0:
      return v + (gridDimensions_[Di_] * shiftY - shiftX); // V(b)::{c}
    case 1:
      return v + gridDimensions_[Di_] * shiftY; // V(b)::{d}
    case 2:
      return v - gridDimensions_[Di_] * decimation_; // V(d)::{b}
    case 3:
      return v - shiftX; // V(d)::{c}
  }
  return -1;
}

template<size_t card>
SimplexId MultiresTriangulation<card>::getInvertVertexNeighbor2dABCD(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY) const {
  // V(abcd)=V(d)::{b,c}+V(c)::{b,d}+V(a)::{b}+V(b)::{c}
  if(neighborId == v - decimation_)
    return 0;
  else if(neighborId == v - gridDimensions_[Di_] * decimation_)
    return 1;
  else if(neighborId == v + (shiftX - gridDimensions_[Di_] * decimation_))
    return 2;
  else if(neighborId == v + shiftX)
    return 3;
  else if(neighborId == v + gridDimensions_[Di_] * shiftY)
    return 4;
  else if(neighborId == v + (gridDimensions_[Di_] * shiftY - decimation_))
    return 5;
  return -1;
}

template<size_t card>
SimplexId
  MultiresTriangulation<card>::getVertexNeighbor2dABCD(const SimplexId v,
                                                 const int id,
                                                 const SimplexId shiftX,
                                                 const SimplexId shiftY) const {
  // V(abcd)=V(d)::{b,c}+V(c)::{b,d}+V(a)::{b}+V(b)::{c}
  switch(id) {
    case 0:
      return v - decimation_;
    case 1:
      return v - gridDimensions_[Di_] * decimation_;
    case 2:
      return v + (shiftX - gridDimensions_[Di_] * decimation_);
    case 3:
      return v + shiftX;
    case 4:
      return v + gridDimensions_[Di_] * shiftY;
    case 5:
      return v + (gridDimensions_[Di_] * shiftY - decimation_);
  }
  return -1;
}

template<size_t card>
SimplexId MultiresTriangulation<card>::getVertexNeighborNumber(
  const SimplexId &vertexId) const {

#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#endif

  switch(vertexPositions_[vertexId]) {
    case VertexPosition::CENTER_3D:
      return 14;
    case VertexPosition::FRONT_FACE_3D:
    case VertexPosition::BACK_FACE_3D:
    case VertexPosition::TOP_FACE_3D:
    case VertexPosition::BOTTOM_FACE_3D:
    case VertexPosition::LEFT_FACE_3D:
    case VertexPosition::RIGHT_FACE_3D:
      return 10;
    case VertexPosition::TOP_FRONT_EDGE_3D: // ab
    case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
    case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
    case VertexPosition::LEFT_BACK_EDGE_3D: // eg
    case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
    case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
      return 8;
    case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
    case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
      return 7;
    case VertexPosition::TOP_BACK_EDGE_3D: // ef
    case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
    case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
    case VertexPosition::TOP_LEFT_EDGE_3D: // ae
    case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
    case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
    case VertexPosition::CENTER_2D:
      return 6;
    case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
    case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
    case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
    case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
    case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
    case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
    case VertexPosition::TOP_EDGE_2D:
    case VertexPosition::BOTTOM_EDGE_2D:
    case VertexPosition::LEFT_EDGE_2D:
    case VertexPosition::RIGHT_EDGE_2D:
      return 4;
    case VertexPosition::TOP_RIGHT_CORNER_2D: // b
    case VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
      return 3;
    case VertexPosition::TOP_LEFT_CORNER_2D: // a
    case VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
    case VertexPosition::CENTER_1D:
      return 2;
    case VertexPosition::LEFT_CORNER_1D:
    case VertexPosition::RIGHT_CORNER_1D:
      return 1;
  }

  return -1;
}

template<size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborA(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(a)={b,c,e,g}
  if(neighborId == v + shiftX)
    return 0; // b
  else if(neighborId == v + vshift_[0] * shiftY)
    return 1; // c
  else if(neighborId == v + vshift_[1] * shiftZ)
    return 2; // e
  else if(neighborId == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ))
    return 3; // g
  return -1;
}

template<size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborA(const SimplexId v,
                                            const int id,
                                            const SimplexId shiftX,
                                            const SimplexId shiftY,
                                            const SimplexId shiftZ) const {
  // V(a)={b,c,e,g}
  switch(id) {
    case 0:
      return v + shiftX; // b
    case 1:
      return v + vshift_[0] * shiftY; // c
    case 2:
      return v + vshift_[1] * shiftZ; // e
    case 3:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // g
  }
  return -1;
}

template<size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborB(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(b)={a,c,d,e,f,g,h}
  if(neighborId == v - shiftX)
    return 0; // a
  else if(neighborId == v + (vshift_[0] * shiftY - shiftX))
    return 1; // c
  else if(neighborId == v + vshift_[0] * shiftY)
    return 2; // d
  else if(neighborId == v + (vshift_[1] * shiftZ - shiftX))
    return 3; // e
  else if(neighborId == v + vshift_[1] * shiftZ)
    return 4; // f
  else if(neighborId
          == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - shiftX))
    return 5; // g
  else if(neighborId == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ))
    return 6; // h
  return -1;
}

template<size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborB(const SimplexId v,
                                            const int id,
                                            const SimplexId shiftX,
                                            const SimplexId shiftY,
                                            const SimplexId shiftZ) const {
  // V(b)={a,c,d,e,f,g,h}
  switch(id) {
    case 0:
      return v - shiftX; // a
    case 1:
      return v + (vshift_[0] * shiftY - shiftX); // c
    case 2:
      return v + vshift_[0] * shiftY; // d
    case 3:
      return v + (vshift_[1] * shiftZ - shiftX); // e
    case 4:
      return v + vshift_[1] * shiftZ; // f
    case 5:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - shiftX); // g
    case 6:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // h
  }
  return -1;
}

template<size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborC(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(c)={a,b,d,g}
  if(neighborId == v - vshift_[0] * shiftY)
    return 0; // a
  else if(neighborId == v + (shiftX - vshift_[0] * shiftY))
    return 1; // b
  else if(neighborId == v + shiftX)
    return 2; // d
  else if(neighborId == v + vshift_[1] * shiftZ)
    return 3; // g
  return -1;
}

template<size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborC(const SimplexId v,
                                            const int id,
                                            const SimplexId shiftX,
                                            const SimplexId shiftY,
                                            const SimplexId shiftZ) const {
  // V(c)={a,b,d,g}
  switch(id) {
    case 0:
      return v - vshift_[0] * shiftY; // a
    case 1:
      return v + (shiftX - vshift_[0] * shiftY); // b
    case 2:
      return v + shiftX; // d
    case 3:
      return v + vshift_[1] * shiftZ; // g
  }
  return -1;
}

template<size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborD(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(d)={b,c,g,h}
  if(neighborId == v - vshift_[0] * shiftY)
    return 0; // b
  else if(neighborId == v - shiftX)
    return 1; // c
  else if(neighborId == v + (-shiftX + vshift_[1] * shiftZ))
    return 2; // g
  else if(neighborId == v + vshift_[1] * shiftZ)
    return 3; // h
  return -1;
}

template<size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborD(const SimplexId v,
                                            const int id,
                                            const SimplexId shiftX,
                                            const SimplexId shiftY,
                                            const SimplexId shiftZ) const {
  // V(d)={b,c,g,h}
  switch(id) {
    case 0:
      return v - vshift_[0] * shiftY; // b
    case 1:
      return v - shiftX; // c
    case 2:
      return v + (-shiftX + vshift_[1] * shiftZ); // g
    case 3:
      return v + vshift_[1] * shiftZ; // h
  }
  return -1;
}

template<size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborE(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(e)={a,b,f,g}
  if(neighborId == v - vshift_[1] * shiftZ)
    return 0; // a
  else if(neighborId == v + (shiftX - vshift_[1] * shiftZ))
    return 1; // b
  else if(neighborId == v + shiftX)
    return 2; // f
  else if(neighborId == v + vshift_[0] * shiftY)
    return 3; // g
  return -1;
}

template<size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborE(const SimplexId v,
                                            const int id,
                                            const SimplexId shiftX,
                                            const SimplexId shiftY,
                                            const SimplexId shiftZ) const {
  // V(e)={a,b,f,g}
  switch(id) {
    case 0:
      return v - vshift_[1] * shiftZ; // a
    case 1:
      return v + (shiftX - vshift_[1] * shiftZ); // b
    case 2:
      return v + shiftX; // f
    case 3:
      return v + vshift_[0] * shiftY; // g
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborF(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(f)={b,e,g,h}
  if(neighborId == v - vshift_[1] * shiftZ)
    return 0; // b
  else if(neighborId == v - shiftX)
    return 1; // e
  else if(neighborId == v + (vshift_[0] * shiftY - shiftX))
    return 2; // g
  else if(neighborId == v + vshift_[0] * shiftY)
    return 3; // h
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborF(const SimplexId v,
                                            const int id,
                                            const SimplexId shiftX,
                                            const SimplexId shiftY,
                                            const SimplexId shiftZ) const {
  // V(f)={b,e,g,h}
  switch(id) {
    case 0:
      return v - vshift_[1] * shiftZ; // b
    case 1:
      return v - shiftX; // e
    case 2:
      return v + (vshift_[0] * shiftY - shiftX); // g
    case 3:
      return v + vshift_[0] * shiftY; // h
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborG(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(g)={a,b,c,d,e,f,h}
  if(neighborId == v - (vshift_[0] * shiftY + vshift_[1] * shiftZ))
    return 0; // a
  else if(neighborId
          == v + (shiftX - vshift_[0] * shiftY - vshift_[1] * shiftZ))
    return 1; // b
  else if(neighborId == v - vshift_[1] * shiftZ)
    return 2; // c
  else if(neighborId == v + (shiftX - vshift_[1] * shiftZ))
    return 3; // d
  else if(neighborId == v - vshift_[0] * shiftY)
    return 4; // e
  else if(neighborId == v + (shiftX - vshift_[0] * shiftY))
    return 5; // f
  else if(neighborId == v + shiftX)
    return 6; // h
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborG(const SimplexId v,
                                            const int id,
                                            const SimplexId shiftX,
                                            const SimplexId shiftY,
                                            const SimplexId shiftZ) const {
  // V(g)={a,b,c,d,e,f,h}
  switch(id) {
    case 0:
      return v - (vshift_[0] * shiftY + vshift_[1] * shiftZ); // a
    case 1:
      return v + (shiftX - vshift_[0] * shiftY - vshift_[1] * shiftZ); // b
    case 2:
      return v - vshift_[1] * shiftZ; // c
    case 3:
      return v + (shiftX - vshift_[1] * shiftZ); // d
    case 4:
      return v - vshift_[0] * shiftY; // e
    case 5:
      return v + (shiftX - vshift_[0] * shiftY); // f
    case 6:
      return v + shiftX; // h
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborH(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(h)={b,d,f,g}
  if(neighborId == v - (vshift_[0] * shiftY + vshift_[1] * shiftZ))
    return 0; // b
  else if(neighborId == v - vshift_[1] * shiftZ)
    return 1; // d
  else if(neighborId == v - vshift_[0] * shiftY)
    return 2; // f
  else if(neighborId == v - shiftX)
    return 3; // g
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborH(const SimplexId v,
                                            const int id,
                                            const SimplexId shiftX,
                                            const SimplexId shiftY,
                                            const SimplexId shiftZ) const {
  // V(h)={b,d,f,g}
  switch(id) {
    case 0:
      return v - (vshift_[0] * shiftY + vshift_[1] * shiftZ); // b
    case 1:
      return v - vshift_[1] * shiftZ; // d
    case 2:
      return v - vshift_[0] * shiftY; // f
    case 3:
      return v - shiftX; // g
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborAB(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(ab)=V(b)+V(a)::{b}
  if(neighborId == v - decimation_)
    return 0; // a
  else if(neighborId == v + (vshift_[0] * shiftY - decimation_))
    return 1; // c
  else if(neighborId == v + vshift_[0] * shiftY)
    return 2; // d
  else if(neighborId == v + (vshift_[1] * shiftZ - decimation_))
    return 3; // e
  else if(neighborId == v + vshift_[1] * shiftZ)
    return 4; // f
  else if(neighborId
          == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - decimation_))
    return 5; // g
  else if(neighborId == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ))
    return 6; // h
  else if(neighborId == v + shiftX)
    return 7; // V(a)::{b}
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborAB(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ) const {
  // V(ab)=V(b)+V(a)::{b}
  switch(id) {
    case 0:
      return v - decimation_; // a
    case 1:
      return v + (vshift_[0] * shiftY - decimation_); // c
    case 2:
      return v + vshift_[0] * shiftY; // d
    case 3:
      return v + (vshift_[1] * shiftZ - decimation_); // e
    case 4:
      return v + vshift_[1] * shiftZ; // f
    case 5:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - decimation_); // g
    case 6:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // h
    case 7:
      return v + shiftX; // V(a)::{b}
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborCD(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(cd)=V(d)+V(c)::{b,d}
  if(neighborId == v - vshift_[0] * shiftY)
    return 0; // b
  else if(neighborId == v - decimation_)
    return 1; // c
  else if(neighborId == v + (vshift_[1] * shiftZ - decimation_))
    return 2; // g
  else if(neighborId == v + vshift_[1] * shiftZ)
    return 3; // h
  else if(neighborId == v + (shiftX - vshift_[0] * shiftY))
    return 4; // V(c)::{b}
  else if(neighborId == v + shiftX)
    return 5; // V(c)::{d}
  return -1;
}
template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborCD(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ) const {
  // V(cd)=V(d)+V(c)::{b,d}
  switch(id) {
    case 0:
      return v - vshift_[0] * shiftY; // b
    case 1:
      return v - decimation_; // c
    case 2:
      return v + (vshift_[1] * shiftZ - decimation_); // g
    case 3:
      return v + vshift_[1] * shiftZ; // h
    case 4:
      return v + (shiftX - vshift_[0] * shiftY); // V(c)::{b}
    case 5:
      return v + shiftX; // V(c)::{d}
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborEF(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(ef)=V(f)+V(e)::{b,f}
  if(neighborId == v - vshift_[1] * shiftZ)
    return 0; // b
  else if(neighborId == v - decimation_)
    return 1; // e
  else if(neighborId == v + (vshift_[0] * shiftY - decimation_))
    return 2; // g
  else if(neighborId == v + vshift_[0] * shiftY)
    return 3; // h
  else if(neighborId == v + (shiftX - vshift_[1] * shiftZ))
    return 4; // V(e)::{b}
  else if(neighborId == v + shiftX)
    return 5; // V(e)::{f}
  return -1;
}
template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborEF(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ) const {
  // V(ef)=V(f)+V(e)::{b,f}
  switch(id) {
    case 0:
      return v - vshift_[1] * shiftZ; // b
    case 1:
      return v - decimation_; // e
    case 2:
      return v + (vshift_[0] * shiftY - decimation_); // g
    case 3:
      return v + vshift_[0] * shiftY; // h
    case 4:
      return v + (shiftX - vshift_[1] * shiftZ); // V(e)::{b}
    case 5:
      return v + shiftX; // V(e)::{f}
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborGH(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(gh)=V(g)+V(h)::{g}
  if(neighborId == v - (vshift_[0] * shiftY + vshift_[1] * shiftZ))
    return 0; // a
  else if(neighborId
          == v + (shiftX - vshift_[0] * shiftY - vshift_[1] * shiftZ))
    return 1; // b
  else if(neighborId == v - vshift_[1] * shiftZ)
    return 2; // c
  else if(neighborId == v + (shiftX - vshift_[1] * shiftZ))
    return 3; // d
  else if(neighborId == v - vshift_[0] * shiftY)
    return 4; // e
  else if(neighborId == v + (shiftX - vshift_[0] * shiftY))
    return 5; // f
  else if(neighborId == v + shiftX)
    return 6; // h
  else if(neighborId == v - decimation_)
    return 7; // V(h)::{g}

  return -1;
}
template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborGH(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ) const {
  // V(gh)=V(g)+V(h)::{g}
  switch(id) {
    case 0:
      return v - (vshift_[0] * shiftY + vshift_[1] * shiftZ); // a
    case 1:
      return v + (shiftX - vshift_[0] * shiftY - vshift_[1] * shiftZ); // b
    case 2:
      return v - vshift_[1] * shiftZ; // c
    case 3:
      return v + (shiftX - vshift_[1] * shiftZ); // d
    case 4:
      return v - vshift_[0] * shiftY; // e
    case 5:
      return v + (shiftX - vshift_[0] * shiftY); // f
    case 6:
      return v + shiftX; // h
    case 7:
      return v - decimation_; // V(h)::{g}
  }

  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborAC(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(ac)=V(c)+V(a)::{c,g}
  if(neighborId == v - vshift_[0] * decimation_)
    return 0; // a
  else if(neighborId == v + (shiftX - vshift_[0] * decimation_))
    return 1; // b
  else if(neighborId == v + shiftX)
    return 2; // d
  else if(neighborId == v + vshift_[1] * shiftZ)
    return 3; // g
  else if(neighborId == v + vshift_[0] * shiftY)
    return 4; // V(a)::{c}
  else if(neighborId == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ))
    return 5; // V(a)::{c}
  return -1;
}
template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborAC(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ) const {
  // V(ac)=V(c)+V(a)::{c,g}
  switch(id) {
    case 0:
      return v - vshift_[0] * decimation_; // a
    case 1:
      return v + (shiftX - vshift_[0] * decimation_); // b
    case 2:
      return v + shiftX; // d
    case 3:
      return v + vshift_[1] * shiftZ; // g
    case 4:
      return v + vshift_[0] * shiftY; // V(a)::{c}
    case 5:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // V(a)::{c}
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborBD(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(bd)=V(b)+V(d)::{b}
  if(neighborId == v - shiftX)
    return 0; // a
  else if(neighborId == v + (vshift_[0] * shiftY - shiftX))
    return 1; // c
  else if(neighborId == v + vshift_[0] * shiftY)
    return 2; // d
  else if(neighborId == v + (vshift_[1] * shiftZ - shiftX))
    return 3; // e
  else if(neighborId == v + vshift_[1] * shiftZ)
    return 4; // f
  else if(neighborId
          == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - shiftX))
    return 5; // g
  else if(neighborId == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ))
    return 6; // h
  else if(neighborId == v - vshift_[0] * decimation_)
    return 7; // V(d)::{b}
  return -1;
}
template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborBD(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ) const {
  // V(bd)=V(b)+V(d)::{b}
  switch(id) {
    case 0:
      return v - shiftX; // a
    case 1:
      return v + (vshift_[0] * shiftY - shiftX); // c
    case 2:
      return v + vshift_[0] * shiftY; // d
    case 3:
      return v + (vshift_[1] * shiftZ - shiftX); // e
    case 4:
      return v + vshift_[1] * shiftZ; // f
    case 5:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - shiftX); // g
    case 6:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // h
    case 7:
      return v - vshift_[0] * decimation_; // V(d)::{b}
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborEG(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(eg)=V(g)+V(e)::{g}
  if(neighborId == v - (vshift_[0] * decimation_ + vshift_[1] * shiftZ))
    return 0; // a
  else if(neighborId
          == v + (shiftX - vshift_[0] * decimation_ - vshift_[1] * shiftZ))
    return 1; // b
  else if(neighborId == v - vshift_[1] * shiftZ)
    return 2; // c
  else if(neighborId == v + (shiftX - vshift_[1] * shiftZ))
    return 3; // d
  else if(neighborId == v - vshift_[0] * decimation_)
    return 4; // e
  else if(neighborId == v + (shiftX - vshift_[0] * decimation_))
    return 5; // f
  else if(neighborId == v + shiftX)
    return 6; // h
  else if(neighborId == v + vshift_[0] * shiftY)
    return 7; // V(e)::{g}
  return -1;
}
template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborEG(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ) const {
  // V(eg)=V(g)+V(e)::{g}
  switch(id) {
    case 0:
      return v - (vshift_[0] * decimation_ + vshift_[1] * shiftZ); // a
    case 1:
      return v + (shiftX - vshift_[0] * decimation_ - vshift_[1] * shiftZ); // b
    case 2:
      return v - vshift_[1] * shiftZ; // c
    case 3:
      return v + (shiftX - vshift_[1] * shiftZ); // d
    case 4:
      return v - vshift_[0] * decimation_; // e
    case 5:
      return v + (shiftX - vshift_[0] * decimation_); // f
    case 6:
      return v + shiftX; // h
    case 7:
      return v + vshift_[0] * shiftY; // V(e)::{g}
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborFH(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(fh)=V(f)+V(h)::{b,f}
  if(neighborId == v - vshift_[1] * shiftZ)
    return 0; // b
  else if(neighborId == v - shiftX)
    return 1; // e
  else if(neighborId == v + (vshift_[0] * shiftY - shiftX))
    return 2; // g
  else if(neighborId == v + vshift_[0] * shiftY)
    return 3; // h
  else if(neighborId == v - (vshift_[0] * decimation_ + vshift_[1] * shiftZ))
    return 4; // V(h)::{b}
  else if(neighborId == v - vshift_[0] * decimation_)
    return 5; // V(h)::{f}
  return -1;
}
template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborFH(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ) const {
  // V(fh)=V(f)+V(h)::{b,f}
  switch(id) {
    case 0:
      return v - vshift_[1] * shiftZ; // b
    case 1:
      return v - shiftX; // e
    case 2:
      return v + (vshift_[0] * shiftY - shiftX); // g
    case 3:
      return v + vshift_[0] * shiftY; // h
    case 4:
      return v - (vshift_[0] * decimation_ + vshift_[1] * shiftZ); // V(h)::{b}
    case 5:
      return v - vshift_[0] * decimation_; // V(h)::{f}
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborAE(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(ae)=V(a)+V(e)::{a,b}
  if(neighborId == v + shiftX)
    return 0; // b
  else if(neighborId == v + vshift_[0] * shiftY)
    return 1; // c
  else if(neighborId == v + vshift_[1] * shiftZ)
    return 2; // e
  else if(neighborId == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ))
    return 3; // g
  else if(neighborId == v - vshift_[1] * decimation_)
    return 4; // V(e)::{a}
  else if(neighborId == v + (shiftX - vshift_[1] * decimation_))
    return 5; // V(e)::{b}
  return -1;
}
template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborAE(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ) const {
  // V(ae)=V(a)+V(e)::{a,b}
  switch(id) {
    case 0:
      return v + shiftX; // b
    case 1:
      return v + vshift_[0] * shiftY; // c
    case 2:
      return v + vshift_[1] * shiftZ; // e
    case 3:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // g
    case 4:
      return v - vshift_[1] * decimation_; // V(e)::{a}
    case 5:
      return v + (shiftX - vshift_[1] * decimation_); // V(e)::{b}
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborBF(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(bf)=V(b)+V(f)::{b}
  if(neighborId == v - shiftX)
    return 0; // a
  else if(neighborId == v + (vshift_[0] * shiftY - shiftX))
    return 1; // c
  else if(neighborId == v + vshift_[0] * shiftY)
    return 2; // d
  else if(neighborId == v + (vshift_[1] * shiftZ - shiftX))
    return 3; // e
  else if(neighborId == v + vshift_[1] * shiftZ)
    return 4; // f
  else if(neighborId
          == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - shiftX))
    return 5; // g
  else if(neighborId == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ))
    return 6; // h
  else if(neighborId == v - vshift_[1] * decimation_)
    return 7; // V(f)::{b}
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborBF(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ) const {
  // V(bf)=V(b)+V(f)::{b}
  switch(id) {
    case 0:
      return v - shiftX; // a
    case 1:
      return v + (vshift_[0] * shiftY - shiftX); // c
    case 2:
      return v + vshift_[0] * shiftY; // d
    case 3:
      return v + (vshift_[1] * shiftZ - shiftX); // e
    case 4:
      return v + vshift_[1] * shiftZ; // f
    case 5:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - shiftX); // g
    case 6:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // h
    case 7:
      return v - vshift_[1] * decimation_; // V(f)::{b}
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborCG(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(cg)=V(g)+V(c)::{g}
  if(neighborId == v - (vshift_[0] * shiftY + vshift_[1] * decimation_))
    return 0; // a
  else if(neighborId
          == v + (shiftX - vshift_[0] * shiftY - vshift_[1] * decimation_))
    return 1; // b
  else if(neighborId == v - vshift_[1] * decimation_)
    return 2; // c
  else if(neighborId == v + (shiftX - vshift_[1] * decimation_))
    return 3; // d
  else if(neighborId == v - vshift_[0] * shiftY)
    return 4; // e
  else if(neighborId == v + (shiftX - vshift_[0] * shiftY))
    return 5; // f
  else if(neighborId == v + shiftX)
    return 6; // h
  else if(neighborId == v + vshift_[1] * shiftZ)
    return 7; // V(c)::{g}
  return -1;
}
template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborCG(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ) const {
  // V(cg)=V(g)+V(c)::{g}
  switch(id) {
    case 0:
      return v - (vshift_[0] * shiftY + vshift_[1] * decimation_); // a
    case 1:
      return v + (shiftX - vshift_[0] * shiftY - vshift_[1] * decimation_); // b
    case 2:
      return v - vshift_[1] * decimation_; // c
    case 3:
      return v + (shiftX - vshift_[1] * decimation_); // d
    case 4:
      return v - vshift_[0] * shiftY; // e
    case 5:
      return v + (shiftX - vshift_[0] * shiftY); // f
    case 6:
      return v + shiftX; // h
    case 7:
      return v + vshift_[1] * shiftZ; // V(c)::{g}
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborDH(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(dh)=V(d)+V(h)::{b,d}
  if(neighborId == v - vshift_[0] * shiftY)
    return 0; // b
  else if(neighborId == v - shiftX)
    return 1; // c
  else if(neighborId == v + (vshift_[1] * shiftZ - shiftX))
    return 2; // g
  else if(neighborId == v + vshift_[1] * shiftZ)
    return 3; // h
  else if(neighborId == v - (vshift_[0] * shiftY + vshift_[1] * decimation_))
    return 4; // V(h)::{b}
  else if(neighborId == v - vshift_[1] * decimation_)
    return 5; // V(h)::{d}
  return -1;
}
template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborDH(const SimplexId v,
                                             const int id,
                                             const SimplexId shiftX,
                                             const SimplexId shiftY,
                                             const SimplexId shiftZ) const {
  // V(dh)=V(d)+V(h)::{b,d}
  switch(id) {
    case 0:
      return v - vshift_[0] * shiftY; // b
    case 1:
      return v - shiftX; // c
    case 2:
      return v + (vshift_[1] * shiftZ - shiftX); // g
    case 3:
      return v + vshift_[1] * shiftZ; // h
    case 4:
      return v - (vshift_[0] * shiftY + vshift_[1] * decimation_); // V(h)::{b}
    case 5:
      return v - vshift_[1] * decimation_; // V(h)::{d}
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborABDC(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(abdc)=V(b)+V(d)::{b}+V(c)::{b}+V(a)::{b}
  if(neighborId == v - decimation_)
    return 0;
  else if(neighborId == v + (vshift_[0] * shiftY - decimation_)) // c
    return 1;
  else if(neighborId == v + vshift_[0] * shiftY) // d
    return 2;
  else if(neighborId == v + (vshift_[1] * shiftZ - decimation_)) // e
    return 3;
  else if(neighborId == v + vshift_[1] * shiftZ) // f
    return 4;
  else if(neighborId
          == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - decimation_)) // g
    return 5;
  else if(neighborId == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ)) // h
    return 6;
  else if(neighborId == v - vshift_[0] * decimation_) // V(d)::{b}
    return 7;
  else if(neighborId == v + (shiftX - vshift_[0] * decimation_)) // V(c)::{b}
    return 8;
  else if(neighborId == v + shiftX) // V(a)::{b}
    return 9;
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborABDC(const SimplexId v,
                                               const int id,
                                               const SimplexId shiftX,
                                               const SimplexId shiftY,
                                               const SimplexId shiftZ) const {
  // V(abdc)=V(b)+V(d)::{b}+V(c)::{b}+V(a)::{b}
  switch(id) {
    case 0:
      return v - decimation_; // a
    case 1:
      return v + (vshift_[0] * shiftY - decimation_); // c
    case 2:
      return v + vshift_[0] * shiftY; // d
    case 3:
      return v + (vshift_[1] * shiftZ - decimation_); // e
    case 4:
      return v + vshift_[1] * shiftZ; // f
    case 5:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - decimation_); // g
    case 6:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // h
    case 7:
      return v - vshift_[0] * decimation_; // V(d)::{b}
    case 8:
      return v + (shiftX - vshift_[0] * decimation_); // V(c)::{b}
    case 9:
      return v + shiftX; // V(a)::{b}
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborEFHG(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(efhg)=V(g)+V(h)::{g}+V(f)::{g,h}
  if(neighborId == v - (vshift_[0] * decimation_ + vshift_[1] * shiftZ))
    return 0; // a
  else if(neighborId
          == v + (shiftX - vshift_[0] * decimation_ - vshift_[1] * shiftZ))
    return 1; // b
  else if(neighborId == v - vshift_[1] * shiftZ)
    return 2; // c
  else if(neighborId == v + (shiftX - vshift_[1] * shiftZ))
    return 3; // d
  else if(neighborId == v - vshift_[0] * decimation_)
    return 4; // e
  else if(neighborId == v + (shiftX - vshift_[0] * decimation_))
    return 5; // f
  else if(neighborId == v + shiftX)
    return 6; // h
  else if(neighborId == v - decimation_)
    return 7; // V(h)::{g}
  else if(neighborId == v + (vshift_[0] * shiftY - decimation_))
    return 8; // V(f)::{g}
  else if(neighborId == v + vshift_[0] * shiftY)
    return 9; // V(f)::{h}
  return -1;
}
template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborEFHG(const SimplexId v,
                                               const int id,
                                               const SimplexId shiftX,
                                               const SimplexId shiftY,
                                               const SimplexId shiftZ) const {
  // V(efhg)=V(g)+V(h)::{g}+V(f)::{g,h}
  switch(id) {
    case 0:
      return v - (vshift_[0] * decimation_ + vshift_[1] * shiftZ); // a
    case 1:
      return v + (shiftX - vshift_[0] * decimation_ - vshift_[1] * shiftZ); // b
    case 2:
      return v - vshift_[1] * shiftZ; // c
    case 3:
      return v + (shiftX - vshift_[1] * shiftZ); // d
    case 4:
      return v - vshift_[0] * decimation_; // e
    case 5:
      return v + (shiftX - vshift_[0] * decimation_); // f
    case 6:
      return v + shiftX; // h
    case 7:
      return v - decimation_; // V(h)::{g}
    case 8:
      return v + (vshift_[0] * shiftY - decimation_); // V(f)::{g}
    case 9:
      return v + vshift_[0] * shiftY; // V(f)::{h}
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborAEGC(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(aegc)=V(g)+V(a)::{c,g}+V(c)::{g}
  if(neighborId == v - (vshift_[0] * decimation_ + vshift_[1] * decimation_))
    return 0; // a
  else if(neighborId
          == v + (shiftX - vshift_[0] * decimation_ - vshift_[1] * decimation_))
    return 1; // b
  else if(neighborId == v - vshift_[1] * decimation_)
    return 2; // c
  else if(neighborId == v + (shiftX - vshift_[1] * decimation_))
    return 3; // d
  else if(neighborId == v - vshift_[0] * decimation_)
    return 4; // e
  else if(neighborId == v + (shiftX - vshift_[0] * decimation_))
    return 5; // f
  else if(neighborId == v + shiftX)
    return 6; // h
  else if(neighborId == v + vshift_[0] * shiftY)
    return 7; // V(a)::{c}
  else if(neighborId == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ))
    return 8; // V(a)::{g}
  else if(neighborId == v + vshift_[1] * shiftZ)
    return 9; // V(c)::{g}
  return -1;
}
template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborAEGC(const SimplexId v,
                                               const int id,
                                               const SimplexId shiftX,
                                               const SimplexId shiftY,
                                               const SimplexId shiftZ) const {
  // V(aegc)=V(g)+V(a)::{c,g}+V(c)::{g}
  switch(id) {
    case 0:
      return v - (vshift_[0] * decimation_ + vshift_[1] * decimation_); // a
    case 1:
      return v
             + (shiftX - vshift_[0] * decimation_
                - vshift_[1] * decimation_); // b
    case 2:
      return v - vshift_[1] * decimation_; // c
    case 3:
      return v + (shiftX - vshift_[1] * decimation_); // d
    case 4:
      return v - vshift_[0] * decimation_; // e
    case 5:
      return v + (shiftX - vshift_[0] * decimation_); // f
    case 6:
      return v + shiftX; // h
    case 7:
      return v + vshift_[0] * shiftY; // V(a)::{c}
    case 8:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // V(a)::{g}
    case 9:
      return v + vshift_[1] * shiftZ; // V(c)::{g}
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborBFHD(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(bfhd)=V(b)+V(f)::{b}+V(h)::{b}+V(d)::{b}
  if(neighborId == v - shiftX)
    return 0; // a
  else if(neighborId == v + (vshift_[0] * shiftY - shiftX))
    return 1; // c
  else if(neighborId == v + vshift_[0] * shiftY)
    return 2; // d
  else if(neighborId == v + (vshift_[1] * shiftZ - shiftX))
    return 3; // e
  else if(neighborId == v + vshift_[1] * shiftZ)
    return 4; // f
  else if(neighborId
          == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - shiftX))
    return 5; // g
  else if(neighborId == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ))
    return 6; // h
  else if(neighborId == v - vshift_[1] * decimation_)
    return 7; // V(f)::{b}
  else if(neighborId
          == v - (vshift_[0] * decimation_ + vshift_[1] * decimation_))
    return 8; // V(h)::{b}
  else if(neighborId == v - vshift_[0] * decimation_)
    return 9; // V(d)::{b}
  return -1;
}
template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborBFHD(const SimplexId v,
                                               const int id,
                                               const SimplexId shiftX,
                                               const SimplexId shiftY,
                                               const SimplexId shiftZ) const {
  // V(bfhd)=V(b)+V(f)::{b}+V(h)::{b}+V(d)::{b}
  switch(id) {
    case 0:
      return v - shiftX; // a
    case 1:
      return v + (vshift_[0] * shiftY - shiftX); // c
    case 2:
      return v + vshift_[0] * shiftY; // d
    case 3:
      return v + (vshift_[1] * shiftZ - shiftX); // e
    case 4:
      return v + vshift_[1] * shiftZ; // f
    case 5:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - shiftX); // g
    case 6:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // h
    case 7:
      return v - vshift_[1] * decimation_; // V(f)::{b}
    case 8:
      return v
             - (vshift_[0] * decimation_
                + vshift_[1] * decimation_); // V(h)::{b}
    case 9:
      return v - vshift_[0] * decimation_; // V(d)::{b}
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborAEFB(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(aefb)=V(b)+V(a)::{b}+V(e)::{b}+V(f)::{b}
  if(neighborId == v - decimation_)
    return 0; // a
  else if(neighborId == v + (vshift_[0] * shiftY - decimation_))
    return 1; // c
  else if(neighborId == v + vshift_[0] * shiftY)
    return 2; // d
  else if(neighborId == v + (vshift_[1] * shiftZ - decimation_))
    return 3; // e
  else if(neighborId == v + vshift_[1] * shiftZ)
    return 4; // f
  else if(neighborId
          == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - decimation_))
    return 5; // g
  else if(neighborId == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ))
    return 6; // h
  else if(neighborId == v + shiftX)
    return 7; // V(a)::{b}
  else if(neighborId == v + (shiftX - vshift_[1] * decimation_))
    return 8; // V(e)::{b}
  else if(neighborId == v - vshift_[1] * decimation_)
    return 9; // V(f)::{b}
  return -1;
}
template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborAEFB(const SimplexId v,
                                               const int id,
                                               const SimplexId shiftX,
                                               const SimplexId shiftY,
                                               const SimplexId shiftZ) const {
  // V(aefb)=V(b)+V(a)::{b}+V(e)::{b}+V(f)::{b}
  switch(id) {
    case 0:
      return v - decimation_; // a
    case 1:
      return v + (vshift_[0] * shiftY - decimation_); // c
    case 2:
      return v + vshift_[0] * shiftY; // d
    case 3:
      return v + (vshift_[1] * shiftZ - decimation_); // e
    case 4:
      return v + vshift_[1] * shiftZ; // f
    case 5:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - decimation_); // g
    case 6:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // h
    case 7:
      return v + shiftX; // V(a)::{b}
    case 8:
      return v + (shiftX - vshift_[1] * decimation_); // V(e)::{b}
    case 9:
      return v - vshift_[1] * decimation_; // V(f)::{b}
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborGHDC(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(ghdc)=V(g)+V(h)::{g}+V(d)::{g,h}
  if(neighborId == v - (vshift_[0] * shiftY + vshift_[1] * decimation_))
    return 0; // a
  else if(neighborId
          == v + (shiftX - vshift_[0] * shiftY - vshift_[1] * decimation_))
    return 1; // b
  else if(neighborId == v - vshift_[1] * decimation_)
    return 2; // c
  else if(neighborId == v + (shiftX - vshift_[1] * decimation_))
    return 3; // d
  else if(neighborId == v - vshift_[0] * shiftY)
    return 4; // e
  else if(neighborId == v + (shiftX - vshift_[0] * shiftY))
    return 5; // f
  else if(neighborId == v + shiftX)
    return 6; // h
  else if(neighborId == v - decimation_)
    return 7; // V(h)::{g}
  else if(neighborId == v + (vshift_[1] * shiftZ - decimation_))
    return 8; // V(d)::{g}
  else if(neighborId == v + vshift_[1] * shiftZ)
    return 9; // V(d)::{h}
  return -1;
}
template <size_t card>
inline ttk::SimplexId
  MultiresTriangulation<card>::getVertexNeighborGHDC(const SimplexId v,
                                               const int id,
                                               const SimplexId shiftX,
                                               const SimplexId shiftY,
                                               const SimplexId shiftZ) const {
  // V(ghdc)=V(g)+V(h)::{g}+V(d)::{g,h}
  switch(id) {
    case 0:
      return v - (vshift_[0] * shiftY + vshift_[1] * decimation_); // a
    case 1:
      return v + (shiftX - vshift_[0] * shiftY - vshift_[1] * decimation_); // b
    case 2:
      return v - vshift_[1] * decimation_; // c
    case 3:
      return v + (shiftX - vshift_[1] * decimation_); // d
    case 4:
      return v - vshift_[0] * shiftY; // e
    case 5:
      return v + (shiftX - vshift_[0] * shiftY); // f
    case 6:
      return v + shiftX; // h
    case 7:
      return v - decimation_; // V(h)::{g}
    case 8:
      return v + (vshift_[1] * shiftZ - decimation_); // V(d)::{g}
    case 9:
      return v + vshift_[1] * shiftZ; // V(d)::{h}
  }
  return -1;
}

template<size_t card>
ttk::SimplexId MultiresTriangulation<card>::getInvertVertexNeighborABCDEFGH(
  const SimplexId v,
  const SimplexId neighborId,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(abcdefgh)=V(g)+V(d)::{g,h}+V(h)::{g}+V(b)::{c,d,g,h}
  if(neighborId == v - (vshift_[0] * decimation_ + vshift_[1] * decimation_))
    return 0; // a
  else if(neighborId
          == v + (shiftX - vshift_[0] * decimation_ - vshift_[1] * decimation_))
    return 1; // b
  else if(neighborId == v - vshift_[1] * decimation_)
    return 2; // c
  else if(neighborId == v + (shiftX - vshift_[1] * decimation_))
    return 3; // d
  else if(neighborId == v - vshift_[0] * decimation_)
    return 4; // e
  else if(neighborId == v + (shiftX - vshift_[0] * decimation_))
    return 5; // f
  else if(neighborId == v + shiftX)
    return 6; // h
  else if(neighborId == v + (vshift_[1] * shiftZ - decimation_))
    return 7; // V(d)::{g}
  else if(neighborId == v + vshift_[1] * shiftZ)
    return 8; // V(d)::{h}
  else if(neighborId == v - decimation_)
    return 9; // V(h)::{g}
  else if(neighborId == v + (vshift_[0] * shiftY - decimation_))
    return 10; // V(b)::{c}
  else if(neighborId == v + vshift_[0] * shiftY)
    return 11; // V(b)::{d}
  else if(neighborId
          == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ - decimation_))
    return 12;
  else if(neighborId == v + (vshift_[0] * shiftY + vshift_[1] * shiftZ))
    return 13;
  return -1;
}

template <size_t card>
ttk::SimplexId MultiresTriangulation<card>::getVertexNeighborABCDEFGH(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ) const {
  // V(abcdefgh)=V(g)+V(d)::{g,h}+V(h)::{g}+V(b)::{c,d,g,h}
  switch(id) {
    case 0:
      return v - (vshift_[0] * decimation_ + vshift_[1] * decimation_); // a
    case 1:
      return v
             + (shiftX - vshift_[0] * decimation_
                - vshift_[1] * decimation_); // b
    case 2:
      return v - vshift_[1] * decimation_; // c
    case 3:
      return v + (shiftX - vshift_[1] * decimation_); // d
    case 4:
      return v - vshift_[0] * decimation_; // e
    case 5:
      return v + (shiftX - vshift_[0] * decimation_); // f
    case 6:
      return v + shiftX; // h
    case 7:
      return v + (vshift_[1] * shiftZ - decimation_); // V(d)::{g}
    case 8:
      return v + vshift_[1] * shiftZ; // V(d)::{h}
    case 9:
      return v - decimation_; // V(h)::{g}
    case 10:
      return v + (vshift_[0] * shiftY - decimation_); // V(b)::{c}
    case 11:
      return v + vshift_[0] * shiftY; // V(b)::{d}
    case 12:
      return v
             + (vshift_[0] * shiftY + vshift_[1] * shiftZ
                - decimation_); // V(b)::{g}
    case 13:
      return v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // V(b)::{h}
  }
  return -1;
}
template <size_t card>
ttk::SimplexId MultiresTriangulation<card>::getInvertedVertexNeighborABCDEFGH(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  const SimplexId shiftZ,
  SimplexId &invertedLocalNeighbor) const {
  // V(abcdefgh)=V(g)+V(d)::{g,h}+V(h)::{g}+V(b)::{c,d,g,h}
  SimplexId invertedVertexId = -1;
  switch(id) {
    case 0:
      invertedLocalNeighbor = 13;
      invertedVertexId
        = v - (vshift_[0] * decimation_ + vshift_[1] * decimation_); // a
      break;
    case 1:
      invertedLocalNeighbor = 12;
      invertedVertexId
        = v
          + (shiftX - vshift_[0] * decimation_ - vshift_[1] * decimation_); // b
      break;
    case 2:
      invertedLocalNeighbor = 8;
      invertedVertexId = v - vshift_[1] * decimation_; // c
      break;
    case 3:
      invertedLocalNeighbor = 7;
      invertedVertexId = v + (shiftX - vshift_[1] * decimation_); // d
      break;
    case 4:
      invertedLocalNeighbor = 11;
      invertedVertexId = v - vshift_[0] * decimation_; // e
      break;
    case 5:
      invertedLocalNeighbor = 10;
      invertedVertexId = v + (shiftX - vshift_[0] * decimation_); // f
      break;
    case 6:
      invertedLocalNeighbor = 9;
      invertedVertexId = v + shiftX; // h
      break;
    case 7:
      invertedLocalNeighbor = 3;
      invertedVertexId = v + (vshift_[1] * shiftZ - decimation_); // V(d)::{g}
      break;
    case 8:
      invertedLocalNeighbor = 2;
      invertedVertexId = v + vshift_[1] * shiftZ; // V(d)::{h}
      break;
    case 9:
      invertedLocalNeighbor = 6;
      invertedVertexId = v - decimation_; // V(h)::{g}
      break;
    case 10:
      invertedLocalNeighbor = 5;
      invertedVertexId = v + (vshift_[0] * shiftY - decimation_); // V(b)::{c}
      break;
    case 11:
      invertedLocalNeighbor = 4;
      invertedVertexId = v + vshift_[0] * shiftY; // V(b)::{d}
      break;
    case 12:
      invertedLocalNeighbor = 1;
      invertedVertexId = v
                         + (vshift_[0] * shiftY + vshift_[1] * shiftZ
                            - decimation_); // V(b)::{g}
      break;
    case 13:
      invertedLocalNeighbor = 0;
      invertedVertexId
        = v + (vshift_[0] * shiftY + vshift_[1] * shiftZ); // V(b)::{h}
      break;
  }
  return invertedVertexId;
}
template<size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborABDC(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 0:
      invertedLocalNeighbor = 6;
      break;
    case 1:
      invertedLocalNeighbor = 5;
      break;
    case 2:
      invertedLocalNeighbor = 4;
      break;
    case 3:
      invertedLocalNeighbor = 3;
      break;
    case 4:
      invertedLocalNeighbor = 2;
      break;
    case 5:
      invertedLocalNeighbor = 1;
      break;
    case 6:
      invertedLocalNeighbor = 0;
      break;
    case 9:
      invertedLocalNeighbor = 9;
      break;
    case 10:
      invertedLocalNeighbor = 8;
      break;
    case 11:
      invertedLocalNeighbor = 7;
      break;
  }
}

template<size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborEFHG(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 13:
      invertedLocalNeighbor = 0;
      break;
    case 12:
      invertedLocalNeighbor = 1;
      break;
    case 8:
      invertedLocalNeighbor = 2;
      break;
    case 7:
      invertedLocalNeighbor = 3;
      break;
    case 11:
      invertedLocalNeighbor = 4;
      break;
    case 10:
      invertedLocalNeighbor = 5;
      break;
    case 9:
      invertedLocalNeighbor = 6;
      break;
    case 6:
      invertedLocalNeighbor = 7;
      break;
    case 5:
      invertedLocalNeighbor = 8;
      break;
    case 4:
      invertedLocalNeighbor = 9;
      break;
  }
}

template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborAEFB(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 6:
      invertedLocalNeighbor = 0;
      break;
    case 5:
      invertedLocalNeighbor = 1;
      break;
    case 4:
      invertedLocalNeighbor = 2;
      break;
    case 3:
      invertedLocalNeighbor = 3;
      break;
    case 2:
      invertedLocalNeighbor = 4;
      break;
    case 1:
      invertedLocalNeighbor = 5;
      break;
    case 0:
      invertedLocalNeighbor = 6;
      break;
    case 9:
      invertedLocalNeighbor = 7;
      break;
    case 7:
      invertedLocalNeighbor = 8;
      break;
    case 8:
      invertedLocalNeighbor = 9;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborGHDC(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 13:
      invertedLocalNeighbor = 0;
      break;
    case 12:
      invertedLocalNeighbor = 1;
      break;
    case 8:
      invertedLocalNeighbor = 2;
      break;
    case 7:
      invertedLocalNeighbor = 3;
      break;
    case 11:
      invertedLocalNeighbor = 4;
      break;
    case 10:
      invertedLocalNeighbor = 5;
      break;
    case 9:
      invertedLocalNeighbor = 6;
      break;
    case 6:
      invertedLocalNeighbor = 7;
      break;
    case 3:
      invertedLocalNeighbor = 8;
      break;
    case 2:
      invertedLocalNeighbor = 9;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborAEGC(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 13:
      invertedLocalNeighbor = 0;
      break;
    case 12:
      invertedLocalNeighbor = 1;
      break;
    case 8:
      invertedLocalNeighbor = 2;
      break;
    case 7:
      invertedLocalNeighbor = 3;
      break;
    case 11:
      invertedLocalNeighbor = 4;
      break;
    case 10:
      invertedLocalNeighbor = 5;
      break;
    case 9:
      invertedLocalNeighbor = 6;
      break;
    case 4:
      invertedLocalNeighbor = 7;
      break;
    case 0:
      invertedLocalNeighbor = 8;
      break;
    case 2:
      invertedLocalNeighbor = 9;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborBFHD(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 6:
      invertedLocalNeighbor = 0;
      break;
    case 5:
      invertedLocalNeighbor = 1;
      break;
    case 4:
      invertedLocalNeighbor = 2;
      break;
    case 3:
      invertedLocalNeighbor = 3;
      break;
    case 2:
      invertedLocalNeighbor = 4;
      break;
    case 1:
      invertedLocalNeighbor = 5;
      break;
    case 0:
      invertedLocalNeighbor = 6;
      break;
    case 8:
      invertedLocalNeighbor = 7;
      break;
    case 13:
      invertedLocalNeighbor = 8;
      break;
    case 11:
      invertedLocalNeighbor = 9;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborAB(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 6:
      invertedLocalNeighbor = 0;
      break;
    case 5:
      invertedLocalNeighbor = 1;
      break;
    case 4:
      invertedLocalNeighbor = 2;
      break;
    case 3:
      invertedLocalNeighbor = 3;
      break;
    case 2:
      invertedLocalNeighbor = 4;
      break;
    case 1:
      invertedLocalNeighbor = 5;
      break;
    case 0:
      invertedLocalNeighbor = 6;
      break;
    case 9:
      invertedLocalNeighbor = 7;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborEF(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 8:
      invertedLocalNeighbor = 0;
      break;
    case 6:
      invertedLocalNeighbor = 1;
      break;
    case 5:
      invertedLocalNeighbor = 2;
      break;
    case 4:
      invertedLocalNeighbor = 3;
      break;
    case 7:
      invertedLocalNeighbor = 4;
      break;
    case 9:
      invertedLocalNeighbor = 5;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborCD(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 11:
      invertedLocalNeighbor = 0;
      break;
    case 6:
      invertedLocalNeighbor = 1;
      break;
    case 3:
      invertedLocalNeighbor = 2;
      break;
    case 2:
      invertedLocalNeighbor = 3;
      break;
    case 10:
      invertedLocalNeighbor = 4;
      break;
    case 9:
      invertedLocalNeighbor = 5;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborGH(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 13:
      invertedLocalNeighbor = 0;
      break;
    case 12:
      invertedLocalNeighbor = 1;
      break;
    case 8:
      invertedLocalNeighbor = 2;
      break;
    case 7:
      invertedLocalNeighbor = 3;
      break;
    case 11:
      invertedLocalNeighbor = 4;
      break;
    case 10:
      invertedLocalNeighbor = 5;
      break;
    case 9:
      invertedLocalNeighbor = 6;
      break;
    case 6:
      invertedLocalNeighbor = 7;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborAC(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 11:
      invertedLocalNeighbor = 0;
      break;
    case 10:
      invertedLocalNeighbor = 1;
      break;
    case 9:
      invertedLocalNeighbor = 2;
      break;
    case 2:
      invertedLocalNeighbor = 3;
      break;
    case 4:
      invertedLocalNeighbor = 4;
      break;
    case 0:
      invertedLocalNeighbor = 5;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborEG(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 13:
      invertedLocalNeighbor = 0;
      break;
    case 12:
      invertedLocalNeighbor = 1;
      break;
    case 8:
      invertedLocalNeighbor = 2;
      break;
    case 7:
      invertedLocalNeighbor = 3;
      break;
    case 11:
      invertedLocalNeighbor = 4;
      break;
    case 10:
      invertedLocalNeighbor = 5;
      break;
    case 9:
      invertedLocalNeighbor = 6;
      break;
    case 4:
      invertedLocalNeighbor = 7;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborAE(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 9:
      invertedLocalNeighbor = 0;
      break;
    case 4:
      invertedLocalNeighbor = 1;
      break;
    case 2:
      invertedLocalNeighbor = 2;
      break;
    case 0:
      invertedLocalNeighbor = 3;
      break;
    case 8:
      invertedLocalNeighbor = 4;
      break;
    case 7:
      invertedLocalNeighbor = 5;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborCG(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 13:
      invertedLocalNeighbor = 0;
      break;
    case 12:
      invertedLocalNeighbor = 1;
      break;
    case 8:
      invertedLocalNeighbor = 2;
      break;
    case 7:
      invertedLocalNeighbor = 3;
      break;
    case 11:
      invertedLocalNeighbor = 4;
      break;
    case 10:
      invertedLocalNeighbor = 5;
      break;
    case 9:
      invertedLocalNeighbor = 6;
      break;
    case 2:
      invertedLocalNeighbor = 7;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborBD(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 6:
      invertedLocalNeighbor = 0;
      break;
    case 5:
      invertedLocalNeighbor = 1;
      break;
    case 4:
      invertedLocalNeighbor = 2;
      break;
    case 3:
      invertedLocalNeighbor = 3;
      break;
    case 2:
      invertedLocalNeighbor = 4;
      break;
    case 1:
      invertedLocalNeighbor = 5;
      break;
    case 0:
      invertedLocalNeighbor = 6;
      break;
    case 11:
      invertedLocalNeighbor = 7;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborFH(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 8:
      invertedLocalNeighbor = 0;
      break;
    case 6:
      invertedLocalNeighbor = 1;
      break;
    case 5:
      invertedLocalNeighbor = 2;
      break;
    case 4:
      invertedLocalNeighbor = 3;
      break;
    case 13:
      invertedLocalNeighbor = 4;
      break;
    case 11:
      invertedLocalNeighbor = 5;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborBF(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 6:
      invertedLocalNeighbor = 0;
      break;
    case 5:
      invertedLocalNeighbor = 1;
      break;
    case 4:
      invertedLocalNeighbor = 2;
      break;
    case 3:
      invertedLocalNeighbor = 3;
      break;
    case 2:
      invertedLocalNeighbor = 4;
      break;
    case 1:
      invertedLocalNeighbor = 5;
      break;
    case 0:
      invertedLocalNeighbor = 6;
      break;
    case 8:
      invertedLocalNeighbor = 7;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborDH(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 11:
      invertedLocalNeighbor = 0;
      break;
    case 6:
      invertedLocalNeighbor = 1;
      break;
    case 3:
      invertedLocalNeighbor = 2;
      break;
    case 2:
      invertedLocalNeighbor = 3;
      break;
    case 13:
      invertedLocalNeighbor = 4;
      break;
    case 8:
      invertedLocalNeighbor = 5;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborA(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 9:
      invertedLocalNeighbor = 0;
      break;
    case 4:
      invertedLocalNeighbor = 1;
      break;
    case 2:
      invertedLocalNeighbor = 2;
      break;
    case 0:
      invertedLocalNeighbor = 3;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborE(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 8:
      invertedLocalNeighbor = 0;
      break;
    case 7:
      invertedLocalNeighbor = 1;
      break;
    case 9:
      invertedLocalNeighbor = 2;
      break;
    case 4:
      invertedLocalNeighbor = 3;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborC(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 11:
      invertedLocalNeighbor = 0;
      break;
    case 10:
      invertedLocalNeighbor = 1;
      break;
    case 9:
      invertedLocalNeighbor = 2;
      break;
    case 2:
      invertedLocalNeighbor = 3;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborG(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 13:
      invertedLocalNeighbor = 0;
      break;
    case 12:
      invertedLocalNeighbor = 1;
      break;
    case 8:
      invertedLocalNeighbor = 2;
      break;
    case 7:
      invertedLocalNeighbor = 3;
      break;
    case 11:
      invertedLocalNeighbor = 4;
      break;
    case 10:
      invertedLocalNeighbor = 5;
      break;
    case 9:
      invertedLocalNeighbor = 6;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborB(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 6:
      invertedLocalNeighbor = 0;
      break;
    case 5:
      invertedLocalNeighbor = 1;
      break;
    case 4:
      invertedLocalNeighbor = 2;
      break;
    case 3:
      invertedLocalNeighbor = 3;
      break;
    case 2:
      invertedLocalNeighbor = 4;
      break;
    case 1:
      invertedLocalNeighbor = 5;
      break;
    case 0:
      invertedLocalNeighbor = 6;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborF(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 8:
      invertedLocalNeighbor = 0;
      break;
    case 6:
      invertedLocalNeighbor = 1;
      break;
    case 5:
      invertedLocalNeighbor = 2;
      break;
    case 4:
      invertedLocalNeighbor = 3;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborD(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 11:
      invertedLocalNeighbor = 0;
      break;
    case 6:
      invertedLocalNeighbor = 1;
      break;
    case 3:
      invertedLocalNeighbor = 2;
      break;
    case 2:
      invertedLocalNeighbor = 3;
      break;
  }
}
template <size_t card>
void MultiresTriangulation<card>::getInvertedLocalNeighborH(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 13:
      invertedLocalNeighbor = 0;
      break;
    case 8:
      invertedLocalNeighbor = 1;
      break;
    case 11:
      invertedLocalNeighbor = 2;
      break;
    case 6:
      invertedLocalNeighbor = 3;
      break;
  }
}

template <size_t card>
void MultiresTriangulation<card>::vertexToPosition2d(
  const SimplexId vertex, std::array<SimplexId, 3> &p) const {
  p[0] = vertex % vshift_[0];
  p[1] = vertex / vshift_[0];
}

template <size_t card>
void MultiresTriangulation<card>::vertexToPosition(
  const SimplexId vertex, std::array<SimplexId, 3> &p) const {
  p[0] = vertex % vshift_[0];
  p[1] = (vertex % vshift_[1]) / vshift_[0];
  p[2] = vertex / vshift_[1];
}

template <size_t card>
bool MultiresTriangulation<card>::isInTriangulation(const SimplexId vertexId) const {
  bool is_in_triangulation = false;
  if(dimensionality_ == 1) {
    is_in_triangulation = ((vertexId % decimation_) == 0);
  } else if(dimensionality_ == 2) {
    std::array<SimplexId, 3> p{};
    vertexToPosition2d(vertexId, p);
    is_in_triangulation
      = ((p[0] % decimation_) == 0) and ((p[1] % decimation_) == 0);
  } else if(dimensionality_ == 3) {
    std::array<SimplexId, 3> p{};
    vertexToPosition(vertexId, p);
    is_in_triangulation = ((p[0] % decimation_) == 0)
                          and ((p[1] % decimation_) == 0)
                          and ((p[2] % decimation_) == 0);
  } else {
    this->printErr("Unknown dimension " + std::to_string(dimensionality_));
  }
  return is_in_triangulation;
}

template <size_t card>
SimplexId
  MultiresTriangulation<card>::localToGlobalVertexId(const SimplexId localId) const {
  // ID in the decimated grid TO ID in the highest resolution grid
  if(decimation_ == 1)
    return localId;
  SimplexId pLocal[3];
  pLocal[0] = localId % gridDecimatedDimensions_[Di_];
  pLocal[1]
    = (localId
       % (gridDecimatedDimensions_[Di_] * gridDecimatedDimensions_[Dj_]))
      / gridDecimatedDimensions_[Di_];
  pLocal[2]
    = localId / (gridDecimatedDimensions_[Di_] * gridDecimatedDimensions_[Dj_]);

  SimplexId p0 = pLocal[0] * decimation_;
  SimplexId p1 = pLocal[1] * decimation_;
  SimplexId p2 = pLocal[2] * decimation_;

  if(nbvoxels_[0] % decimation_
     and pLocal[0] == gridDecimatedDimensions_[Di_] - 1) {
    p0 -= decimation_ - nbvoxels_[0] % decimation_;
  }
  if(nbvoxels_[1] % decimation_
     and pLocal[1] == gridDecimatedDimensions_[Dj_] - 1) {
    p1 -= decimation_ - nbvoxels_[1] % decimation_;
  }
  if(nbvoxels_[2] % decimation_
     and pLocal[2] == gridDecimatedDimensions_[2] - 1) {
    p2 -= decimation_ - nbvoxels_[2] % decimation_;
  }
  SimplexId const globalId = p0 + p1 * gridDimensions_[0]
                             + p2 * gridDimensions_[0] * gridDimensions_[1];

  return globalId;
}

// Gets impacted vertices v0 and v1 when a new point pops up on an edge.
// Fills the v0 and v1 array with 3 values:
//      v[0] : the local neighbor index of v0 for vertexId
//      v[1] : the global index of v0
//      v[2] : the local neighbor index of vertexId for v0
template <size_t card>
void MultiresTriangulation<card>::getImpactedVertices(SimplexId vertexId,
                                                SimplexId v0[3],
                                                SimplexId v1[3]) {

  SimplexId localNeighborId0 = -1, localNeighborId1 = -1;
  if(dimensionality_ == 3) {
    std::array<SimplexId, 3> p{};
    vertexToPosition(vertexId, p);

    if(0 < p[0] and p[0] < nbvoxels_[0]) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          getImpactedVerticesABCDEFGH(p, localNeighborId0,
                                      localNeighborId1); //
        else if(p[2] == 0)
          getImpactedVerticesABDC(p, localNeighborId0, localNeighborId1); //
        else
          getImpactedVerticesEFHG(p, localNeighborId0, localNeighborId1); //
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          getImpactedVerticesAEFB(p, localNeighborId0, localNeighborId1); //
        else if(p[2] == 0)
          getImpactedVerticesAB(p, localNeighborId0, localNeighborId1); //
        else
          getImpactedVerticesEF(p, localNeighborId0, localNeighborId1); //
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          getImpactedVerticesGHDC(p, localNeighborId0, localNeighborId1); //
        else if(p[2] == 0)
          getImpactedVerticesCD(p, localNeighborId0, localNeighborId1); //
        else
          getImpactedVerticesGH(p, localNeighborId0, localNeighborId1); //
      }
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          getImpactedVerticesAEGC(p, localNeighborId0, localNeighborId1); //
        else if(p[2] == 0)
          getImpactedVerticesAC(p, localNeighborId0,
                                localNeighborId1); //
        else
          getImpactedVerticesEG(p, localNeighborId0,
                                localNeighborId1); //
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          getImpactedVerticesAE(p, localNeighborId0,
                                localNeighborId1); //
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          getImpactedVerticesCG(p, localNeighborId0,
                                localNeighborId1); //
      }
    } else {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          getImpactedVerticesBFHD(p, localNeighborId0,
                                  localNeighborId1); //
        else if(p[2] == 0)
          getImpactedVerticesBD(p, localNeighborId0,
                                localNeighborId1); //
        else
          getImpactedVerticesFH(p, localNeighborId0,
                                localNeighborId1); //
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          getImpactedVerticesBF(p, localNeighborId0,
                                localNeighborId1); //
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          getImpactedVerticesDH(p, localNeighborId0,
                                localNeighborId1); //
      }
    }

  } else if(dimensionality_ == 2) {
    std::array<SimplexId, 3> p{};
    vertexToPosition2d(vertexId, p);

    if(0 < p[0] and p[0] < nbvoxels_[Di_]) {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_]) {
        getImpactedVertices2dABCD(p, localNeighborId0,
                                  localNeighborId1); // abcd
      } else if(p[1] == 0) {
        getImpactedVertices2dAB(p, localNeighborId0,
                                localNeighborId1); // ab
      } else {
        getImpactedVertices2dCD(p, localNeighborId0,
                                localNeighborId1); // cd
      }
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_]) {
        getImpactedVertices2dAC(p, localNeighborId0,
                                localNeighborId1); // ac
      }
    } else {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_]) {
        getImpactedVertices2dBD(p, localNeighborId0,
                                localNeighborId1); // bd
      }
    }
  } else {
    this->printWrn("getImpactedVertices not implemented for 1D yet");
  }
  v0[0] = localNeighborId0;
  v1[0] = localNeighborId1;
  getVertexNeighbor(vertexId, v0[0], v0[1]);
  getVertexNeighbor(vertexId, v1[0], v1[1]);
  getInvertVertexNeighbor(v0[1], vertexId, v0[2]);
  getInvertVertexNeighbor(v1[1], vertexId, v1[2]);
}

template <size_t card>
char MultiresTriangulation<card>::localNeighborId(SimplexId neighborId,
                                            SimplexId vertexId) {
  this->printMsg("BOUNDARY CASES TO TAKE CARE OF");

  SimplexId pLocalNeighbor[3];
  pLocalNeighbor[0] = neighborId % gridDecimatedDimensions_[Di_];
  pLocalNeighbor[1]
    = (neighborId
       % (gridDecimatedDimensions_[Di_] * gridDecimatedDimensions_[Dj_]))
      / gridDecimatedDimensions_[Di_];
  pLocalNeighbor[2]
    = neighborId
      / (gridDecimatedDimensions_[Di_] * gridDecimatedDimensions_[Dj_]);

  SimplexId pLocalVertex[3];
  pLocalVertex[0] = vertexId % gridDecimatedDimensions_[Di_];
  pLocalVertex[1]
    = (vertexId
       % (gridDecimatedDimensions_[Di_] * gridDecimatedDimensions_[Dj_]))
      / gridDecimatedDimensions_[Di_];
  pLocalVertex[2]
    = vertexId
      / (gridDecimatedDimensions_[Di_] * gridDecimatedDimensions_[Dj_]);

  int localNeighbor = -1;

  if(dimensionality_ == 2) {
    if(pLocalNeighbor[0] == pLocalVertex[0]) {
      if(pLocalNeighbor[1] == pLocalVertex[1]) {
        this->printErr("localNeighbor: not a neighbor!");
      } else if(pLocalNeighbor[1] == pLocalVertex[1] + 1) {
        localNeighbor = 4;
      } else if(pLocalNeighbor[1] == pLocalVertex[1] - 1) {
        localNeighbor = 1;
      }
    } else if(pLocalNeighbor[0] == pLocalVertex[0] + 1) {
      if(pLocalNeighbor[1] == pLocalVertex[1]) {
        localNeighbor = 3;
      } else if(pLocalNeighbor[1] == pLocalVertex[1] - 1) {
        localNeighbor = 2;
      } else {
        this->printErr("localNeighbor: not a neighbor!");
      }
    } else if(pLocalNeighbor[0] == pLocalVertex[0] - 1) {
      if(pLocalNeighbor[1] == pLocalVertex[1]) {
        localNeighbor = 0;
      } else if(pLocalNeighbor[1] == pLocalVertex[1] + 1) {
        localNeighbor = 5;
      } else {
        this->printErr("localNeighbor: not a neighbor!");
      }
    } else {
      this->printErr("localNeighbor: not a neighbor!");
    }
  } else if(dimensionality_ == 3) {
    this->printWrn("Local neighbors ids not implemented yet for 3D");
  }

  return localNeighbor;
}

template <size_t card>
int MultiresTriangulation<card>::getInteriorInvertedVertexNeighbor(
  SimplexId vertexId,
  SimplexId localNeighborId,
  SimplexId &invertedVertexId,
  SimplexId &invertedLocalNeighborId) const {

  if(dimensionality_ == 3) {
    std::array<SimplexId, 3> p{};
    vertexToPosition(vertexId, p);
    SimplexId shiftX = decimation_;
    SimplexId shiftY = decimation_;
    SimplexId shiftZ = decimation_;

    if((nbvoxels_[0] % decimation_) and (p[0] + decimation_ > nbvoxels_[0])) {
      shiftX = nbvoxels_[0] % decimation_;
    }
    if((nbvoxels_[1] % decimation_) and (p[1] + decimation_ > nbvoxels_[1])) {
      shiftY = nbvoxels_[1] % decimation_;
    }
    if((nbvoxels_[2] % decimation_) and (p[2] + decimation_ > nbvoxels_[2])) {
      shiftZ = nbvoxels_[2] % decimation_;
    }
    invertedVertexId = getInvertedVertexNeighborABCDEFGH(
      vertexId, localNeighborId, shiftX, shiftY, shiftZ,
      invertedLocalNeighborId); //

    if(decimation_ < p[0] and p[0] < nbvoxels_[0] - shiftX) {
      if(decimation_ < p[1] and p[1] < nbvoxels_[1] - shiftY) {
        if(decimation_ < p[2] and p[2] < nbvoxels_[2] - shiftZ) {
          // do nothing
          return 0;
        } else if(p[2] == decimation_)
          getInvertedLocalNeighborABDC(vertexId, invertedLocalNeighborId);
        // boundary = 0;
        else
          getInvertedLocalNeighborEFHG(vertexId, invertedLocalNeighborId);
        // boundary = 1;
      } else if(p[1] == decimation_) {
        if(decimation_ < p[2] and p[2] < nbvoxels_[2] - shiftZ)
          getInvertedLocalNeighborAEFB(vertexId, invertedLocalNeighborId);
        // boundary = 2;
        else if(p[2] == decimation_)
          getInvertedLocalNeighborAB(vertexId, invertedLocalNeighborId);
        // boundary = 6;
        else
          getInvertedLocalNeighborEF(vertexId, invertedLocalNeighborId);
        // boundary = 7;
      } else {
        if(decimation_ < p[2] and p[2] < nbvoxels_[2] - shiftZ)
          getInvertedLocalNeighborGHDC(vertexId, invertedLocalNeighborId);
        // boundary=3;
        else if(p[2] == decimation_)
          getInvertedLocalNeighborCD(vertexId, invertedLocalNeighborId);
        // boundary = 8;
        else
          getInvertedLocalNeighborGH(vertexId, invertedLocalNeighborId);
        // boundary = 9;
      }
    } else if(p[0] == decimation_) {
      if(decimation_ < p[1] and p[1] < nbvoxels_[1] - shiftY) {
        if(decimation_ < p[2] and p[2] < nbvoxels_[2] - shiftZ)
          getInvertedLocalNeighborAEGC(vertexId, invertedLocalNeighborId);
        // boundary = 4;
        else if(p[2] == decimation_)
          getInvertedLocalNeighborAC(vertexId, invertedLocalNeighborId);
        // boundary = 10;
        else
          getInvertedLocalNeighborEG(vertexId, invertedLocalNeighborId);
        // boundary = 11;
      } else if(p[1] == decimation_) {
        if(decimation_ < p[2] and p[2] < nbvoxels_[2] - shiftZ)
          getInvertedLocalNeighborAE(vertexId, invertedLocalNeighborId);
        // boundary = 12;
        else if(p[2] == decimation_)
          getInvertedLocalNeighborA(vertexId, invertedLocalNeighborId);
        // boundary = 19;
        else
          getInvertedLocalNeighborE(vertexId, invertedLocalNeighborId);
        // boundary = 18;
      } else {
        if(decimation_ < p[2] and p[2] < nbvoxels_[2] - shiftZ)
          // boundary = 13;
          getInvertedLocalNeighborCG(vertexId, invertedLocalNeighborId);
        else if(p[2] == decimation_)
          getInvertedLocalNeighborC(vertexId, invertedLocalNeighborId);
        // boundary = 20;
        else
          getInvertedLocalNeighborG(vertexId, invertedLocalNeighborId);
        // boundary = 21;
      }
    } else {
      if(decimation_ < p[1] and p[1] < nbvoxels_[1] - shiftY) {
        if(decimation_ < p[2] and p[2] < nbvoxels_[2] - shiftZ)
          getInvertedLocalNeighborBFHD(vertexId, invertedLocalNeighborId);
        // boundary = 5;
        else if(p[2] == decimation_)
          getInvertedLocalNeighborBD(vertexId, invertedLocalNeighborId);
        // boundary = 14;
        else
          getInvertedLocalNeighborFH(vertexId, invertedLocalNeighborId);
        // boundary = 15;
      } else if(p[1] == decimation_) {
        if(0 < p[2] and p[2] < nbvoxels_[2] - shiftZ)
          getInvertedLocalNeighborBF(vertexId, invertedLocalNeighborId);
        // boundary = 16;
        else if(p[2] == decimation_)
          getInvertedLocalNeighborB(vertexId, invertedLocalNeighborId);
        // boundary = 22;
        else
          getInvertedLocalNeighborF(vertexId, invertedLocalNeighborId);
        // boundary = 23;
      } else {
        if(decimation_ < p[2] and p[2] < nbvoxels_[2] - shiftZ)
          getInvertedLocalNeighborDH(vertexId, invertedLocalNeighborId);
        // boundary = 17;
        else if(p[2] == decimation_)
          getInvertedLocalNeighborD(vertexId, invertedLocalNeighborId);
        // boundary = 24;
        else
          getInvertedLocalNeighborH(vertexId, invertedLocalNeighborId);
        // boundary = 25;
      }
    }
    return -1;

  } else if(dimensionality_ == 2) {
    std::array<SimplexId, 3> p{};
    vertexToPosition2d(vertexId, p);
    SimplexId shiftX = decimation_;
    SimplexId shiftY = decimation_;

    if((nbvoxels_[Di_] % decimation_)
       and (p[0] + decimation_ > nbvoxels_[Di_])) {
      shiftX = nbvoxels_[0] % decimation_;
    }
    if((nbvoxels_[Dj_] % decimation_)
       and (p[1] + decimation_ > nbvoxels_[Dj_])) {
      shiftY = nbvoxels_[1] % decimation_;
    }

    invertedVertexId = getInvertedVertexNeighbor2dABCD(
      vertexId, localNeighborId, shiftX, shiftY,
      invertedLocalNeighborId); // abcd

    if(0 < p[0] and p[0] < nbvoxels_[Di_] - shiftX) {
      if(decimation_ < p[1] and p[1] < nbvoxels_[Dj_]) {
        // nothing to do here
      } else if(p[1] == decimation_) {
        getInvertedLocalNeighbor2dAB(vertexId, invertedLocalNeighborId);
      } else {
        getInvertedLocalNeighbor2dCD(vertexId, invertedLocalNeighborId); // cd
      }
    } else if(p[0] == decimation_) {
      if(decimation_ < p[1] and p[1] < nbvoxels_[Dj_] - shiftY) {
        getInvertedLocalNeighbor2dAC(vertexId, invertedLocalNeighborId); // ac
      } else if(p[1] == decimation_) {
        getInvertedLocalNeighbor2dA(vertexId, invertedLocalNeighborId); // a
      } else {
        getInvertedLocalNeighbor2dC(vertexId, invertedLocalNeighborId); // c
      }
    } else {
      if(decimation_ < p[1] and p[1] < nbvoxels_[Dj_] - shiftY) {
        getInvertedLocalNeighbor2dBD(vertexId, invertedLocalNeighborId); // bd
      } else if(p[1] == decimation_) {
        getInvertedLocalNeighbor2dB(vertexId, invertedLocalNeighborId); // b
      } else {
        getInvertedLocalNeighbor2dD(vertexId, invertedLocalNeighborId); // d
      }
    }
  } else if(dimensionality_ == 1) {
    // ab
    this->printWrn("NOT TESTED IN 1D");
    if(vertexId > decimation_ and vertexId < nbvoxels_[Di_]) {
      if(localNeighborId == decimation_)
        invertedVertexId = vertexId + decimation_;
      else
        invertedVertexId = vertexId - decimation_;
    } else if(vertexId == decimation_)
      invertedVertexId = vertexId + decimation_;
    else
      invertedVertexId = vertexId - decimation_;
  }
  return 0;
}

template <size_t card>
SimplexId MultiresTriangulation<card>::getInvertedVertexNeighbor2dABCD(
  const SimplexId v,
  const int id,
  const SimplexId shiftX,
  const SimplexId shiftY,
  SimplexId &invertedLocalNeighbor) const {
  SimplexId invertedVertexId = -1;
  switch(id) {
    case 0:
      invertedLocalNeighbor = 3;
      invertedVertexId = v - decimation_;
      break;
    case 1:
      invertedLocalNeighbor = 4;
      invertedVertexId = v - gridDimensions_[Di_] * decimation_;
      break;
    case 2:
      invertedLocalNeighbor = 5;
      invertedVertexId = v + (shiftX - gridDimensions_[Di_] * decimation_);
      break;
    case 3:
      invertedLocalNeighbor = 0;
      invertedVertexId = v + shiftX;
      break;
    case 4:
      invertedLocalNeighbor = 1;
      invertedVertexId = v + gridDimensions_[Di_] * shiftY;
      break;
    case 5:
      invertedLocalNeighbor = 2;
      invertedVertexId = v + (gridDimensions_[Di_] * shiftY - decimation_);
      break;
  }
  return invertedVertexId;
}

template <size_t card>
SimplexId MultiresTriangulation<card>::getInvertedLocalNeighbor2dA(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 0:
      invertedLocalNeighbor = 0;
      break;
    case 1:
      invertedLocalNeighbor = 1;
      break;
  }
  return 0;
}

template <size_t card>
SimplexId MultiresTriangulation<card>::getInvertedLocalNeighbor2dB(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 3:
      invertedLocalNeighbor = 0;
      break;
    case 1:
      invertedLocalNeighbor = 1;
      break;
    case 2:
      invertedLocalNeighbor = 2;
      break;
  }
  return 0;
}

template <size_t card>
SimplexId MultiresTriangulation<card>::getInvertedLocalNeighbor2dC(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 0:
      invertedLocalNeighbor = 0;
      break;
    case 4:
      invertedLocalNeighbor = 1;
      break;
    case 5:
      invertedLocalNeighbor = 2;
      break;
  }
  return 0;
}

template <size_t card>
SimplexId MultiresTriangulation<card>::getInvertedLocalNeighbor2dD(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 3:
      invertedLocalNeighbor = 0;
      break;
    case 4:
      invertedLocalNeighbor = 1;
      break;
  }
  return 0;
}

template <size_t card>
SimplexId MultiresTriangulation<card>::getInvertedLocalNeighbor2dAB(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 3:
      invertedLocalNeighbor = 0;
      break;
    case 2:
      invertedLocalNeighbor = 1;
      break;
    case 1:
      invertedLocalNeighbor = 2;
      break;
    case 0:
      invertedLocalNeighbor = 3;
      break;
  }
  return 0;
}

template <size_t card>
SimplexId MultiresTriangulation<card>::getInvertedLocalNeighbor2dAC(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 4:
      invertedLocalNeighbor = 0;
      break;
    case 5:
      invertedLocalNeighbor = 1;
      break;
    case 0:
      invertedLocalNeighbor = 2;
      break;
    case 1:
      invertedLocalNeighbor = 3;
      break;
  }
  return 0;
}

template <size_t card>
SimplexId MultiresTriangulation<card>::getInvertedLocalNeighbor2dBD(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 2:
      invertedLocalNeighbor = 0;
      break;
    case 1:
      invertedLocalNeighbor = 1;
      break;
    case 4:
      invertedLocalNeighbor = 2;
      break;
    case 3:
      invertedLocalNeighbor = 3;
      break;
  }
  return 0;
}

template <size_t card>
SimplexId MultiresTriangulation<card>::getInvertedLocalNeighbor2dCD(
  SimplexId id, SimplexId &invertedLocalNeighbor) const {
  switch(id) {
    case 3:
      invertedLocalNeighbor = 0;
      break;
    case 2:
      invertedLocalNeighbor = 1;
      break;
    case 1:
      invertedLocalNeighbor = 2;
      break;
    case 0:
      invertedLocalNeighbor = 3;
      break;
  }
  return 0;
}

template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesError(
  const int prev_decim, const std::array<SimplexId, 3> &p) const {
  this->printErr("THIS SHOULDNT HAPPEN");
  this->printErr("previous decimation: " + std::to_string(prev_decim));
  this->printErr("position: " + std::to_string(p[0]) + " "
                 + std::to_string(p[1]) + " " + std::to_string(p[2]));
  this->printErr("grid size: " + std::to_string(gridDimensions_[0]) + " "
                 + std::to_string(gridDimensions_[1]) + " "
                 + std::to_string(gridDimensions_[2]));
}

template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesABCDEFGH(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {

  int const previous_decimation = pow(2, (decimationLevel_ + 1));
  if((p[0] % previous_decimation) and (p[1] % previous_decimation)
     and (p[2] % previous_decimation)) {
    localNeighborId0 = 1;
    localNeighborId1 = 12;

  } else if((p[0] % previous_decimation) and (p[1] % previous_decimation)) {
    localNeighborId0 = 5;
    localNeighborId1 = 10;
  } else if((p[0] % previous_decimation) and (p[2] % previous_decimation)) {
    localNeighborId0 = 3;
    localNeighborId1 = 7;

  } else if((p[1] % previous_decimation) and (p[2] % previous_decimation)) {
    localNeighborId0 = 0;
    localNeighborId1 = 13;

  } else if(p[0] % previous_decimation) {
    localNeighborId0 = 6;
    localNeighborId1 = 9;

  } else if(p[1] % previous_decimation) {
    localNeighborId0 = 4;
    localNeighborId1 = 11;

  } else if(p[2] % previous_decimation) {
    localNeighborId0 = 2;
    localNeighborId1 = 8;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}

template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesABDC(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if((p[0] % previous_decimation) and (p[1] % previous_decimation)) {
    localNeighborId0 = 1;
    localNeighborId1 = 8;

  } else if(p[0] % previous_decimation) {
    localNeighborId0 = 0;
    localNeighborId1 = 9;

  } else if(p[1] % previous_decimation) {
    localNeighborId0 = 2;
    localNeighborId1 = 7;

  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}
template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesEFHG(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if((p[0] % previous_decimation) and (p[1] % previous_decimation)) {
    localNeighborId0 = 5;
    localNeighborId1 = 8;

  } else if(p[0] % previous_decimation) {
    localNeighborId0 = 7;
    localNeighborId1 = 6;

  } else if(p[1] % previous_decimation) {
    localNeighborId0 = 4;
    localNeighborId1 = 9;

  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}
template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesAEGC(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if((p[1] % previous_decimation) and (p[2] % previous_decimation)) {
    localNeighborId0 = 0;
    localNeighborId1 = 8;

  } else if(p[1] % previous_decimation) {
    localNeighborId0 = 4;
    localNeighborId1 = 7;

  } else if(p[2] % previous_decimation) {
    localNeighborId0 = 2;
    localNeighborId1 = 9;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}
template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesBFHD(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if((p[1] % previous_decimation) and (p[2] % previous_decimation)) {
    localNeighborId0 = 8;
    localNeighborId1 = 6;

  } else if(p[1] % previous_decimation) {
    localNeighborId0 = 9;
    localNeighborId1 = 2;

  } else if(p[2] % previous_decimation) {
    localNeighborId0 = 7;
    localNeighborId1 = 4;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}
template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesAEFB(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if((p[0] % previous_decimation) and (p[2] % previous_decimation)) {
    localNeighborId0 = 3;
    localNeighborId1 = 8;

  } else if(p[0] % previous_decimation) {
    localNeighborId0 = 0;
    localNeighborId1 = 7;

  } else if(p[2] % previous_decimation) {
    localNeighborId0 = 9;
    localNeighborId1 = 4;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}
template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesGHDC(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if((p[0] % previous_decimation) and (p[2] % previous_decimation)) {
    localNeighborId0 = 3;
    localNeighborId1 = 8;

  } else if(p[0] % previous_decimation) {
    localNeighborId0 = 7;
    localNeighborId1 = 6;

  } else if(p[2] % previous_decimation) {
    localNeighborId0 = 2;
    localNeighborId1 = 9;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}
template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesAB(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if(p[0] % previous_decimation) {
    localNeighborId0 = 0;
    localNeighborId1 = 7;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}
template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesCD(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if(p[0] % previous_decimation) {
    localNeighborId0 = 1;
    localNeighborId1 = 5;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}
template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesEF(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if(p[0] % previous_decimation) {
    localNeighborId0 = 1;
    localNeighborId1 = 5;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}
template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesGH(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if(p[0] % previous_decimation) {
    localNeighborId0 = 7;
    localNeighborId1 = 6;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}
template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesAC(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if(p[1] % previous_decimation) {
    localNeighborId0 = 0;
    localNeighborId1 = 4;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}
template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesBD(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if(p[1] % previous_decimation) {
    localNeighborId0 = 2;
    localNeighborId1 = 7;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}
template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesEG(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if(p[1] % previous_decimation) {
    localNeighborId0 = 4;
    localNeighborId1 = 7;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}
template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesFH(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if(p[1] % previous_decimation) {
    localNeighborId0 = 5;
    localNeighborId1 = 3;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}
template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesAE(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if(p[2] % previous_decimation) {
    localNeighborId0 = 4;
    localNeighborId1 = 2;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}

template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesBF(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if(p[2] % previous_decimation) {
    localNeighborId0 = 7;
    localNeighborId1 = 4;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}

template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesCG(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if(p[2] % previous_decimation) {
    localNeighborId0 = 2;
    localNeighborId1 = 7;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}

template <size_t card>
void MultiresTriangulation<card>::getImpactedVerticesDH(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if(p[2] % previous_decimation) {
    localNeighborId0 = 5;
    localNeighborId1 = 3;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}

template <size_t card>
void MultiresTriangulation<card>::getImpactedVertices2dABCD(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if((p[0] % previous_decimation) and (p[1] % previous_decimation)) {
    localNeighborId0 = 2;
    localNeighborId1 = 5;

  } else if(p[0] % previous_decimation) {
    localNeighborId0 = 0;
    localNeighborId1 = 3;

  } else if(p[1] % previous_decimation) {
    localNeighborId0 = 1;
    localNeighborId1 = 4;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}

template <size_t card>
void MultiresTriangulation<card>::getImpactedVertices2dAB(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if(p[0] % previous_decimation) {
    localNeighborId0 = 0;
    localNeighborId1 = 3;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}

template <size_t card>
void MultiresTriangulation<card>::getImpactedVertices2dCD(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if(p[0] % previous_decimation) {
    localNeighborId0 = 0;
    localNeighborId1 = 3;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}

template <size_t card>
void MultiresTriangulation<card>::getImpactedVertices2dAC(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if(p[1] % previous_decimation) {
    localNeighborId0 = 0;
    localNeighborId1 = 3;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}
template <size_t card>
void MultiresTriangulation<card>::getImpactedVertices2dBD(
  std::array<SimplexId, 3> &p,
  SimplexId &localNeighborId0,
  SimplexId &localNeighborId1) const {
  int const previous_decimation = pow(2, (decimationLevel_ + 1));

  if(p[1] % previous_decimation) {
    localNeighborId0 = 2;
    localNeighborId1 = 1;
  } else {
    this->getImpactedVerticesError(previous_decimation, p);
  }
}

template <size_t card>
std::vector<SimplexId>
  MultiresTriangulation<card>::getExtendedStar(const SimplexId &vertexId) const {

  std::vector<SimplexId> result;
  if(dimensionality_ == 3) {
    std::array<SimplexId, 3> p{};
    vertexToPosition(vertexId, p);
    SimplexId shiftX = decimation_;
    SimplexId shiftY = decimation_;
    SimplexId shiftZ = decimation_;

    if((nbvoxels_[0] % decimation_) and (p[0] + decimation_ > nbvoxels_[0])) {
      shiftX = nbvoxels_[0] % decimation_;
    }
    if((nbvoxels_[1] % decimation_) and (p[1] + decimation_ > nbvoxels_[1])) {
      shiftY = nbvoxels_[1] % decimation_;
    }
    if((nbvoxels_[2] % decimation_) and (p[2] + decimation_ > nbvoxels_[2])) {
      shiftZ = nbvoxels_[2] % decimation_;
    }

    std::vector<SimplexId> vi;
    vi.push_back(-decimation_);
    vi.push_back(0);
    vi.push_back(shiftX);

    std::vector<SimplexId> vj;
    vj.push_back(-decimation_);
    vj.push_back(0);
    vj.push_back(shiftY);

    std::vector<SimplexId> vk;
    vk.push_back(-decimation_);
    vk.push_back(0);
    vk.push_back(shiftZ);

    if(p[0] == 0)
      vi.erase(vi.begin() + 0);
    if(p[0] == nbvoxels_[0])
      vi.erase(vi.begin() + 2);
    if(p[1] == 0)
      vj.erase(vj.begin() + 0);
    if(p[1] == nbvoxels_[1])
      vj.erase(vj.begin() + 2);
    if(p[2] == 0)
      vk.erase(vk.begin() + 0);
    if(p[2] == nbvoxels_[2])
      vk.erase(vk.begin() + 2);

    for(SimplexId const i : vi) {
      for(SimplexId const j : vj) {
        for(SimplexId const k : vk) {
          if(i != 0 or j != 0 or k != 0)
            result.push_back(vertexId + i + j * vshift_[0] + k * vshift_[1]);
        }
      }
    }

  } else if(dimensionality_ == 2) {
    std::array<SimplexId, 3> p{};
    vertexToPosition2d(vertexId, p);
    SimplexId shiftX = decimation_;
    SimplexId shiftY = decimation_;

    if((nbvoxels_[Di_] % decimation_)
       and (p[0] + decimation_ > nbvoxels_[Di_])) {
      shiftX = nbvoxels_[0] % decimation_;
    }
    if((nbvoxels_[Dj_] % decimation_)
       and (p[1] + decimation_ > nbvoxels_[Dj_])) {
      shiftY = nbvoxels_[1] % decimation_;
    }

    std::vector<SimplexId> vi;
    vi.push_back(-decimation_);
    vi.push_back(0);
    vi.push_back(shiftX);

    std::vector<SimplexId> vj;
    vj.push_back(-decimation_);
    vj.push_back(0);
    vj.push_back(shiftY);

    if(p[0] == 0)
      vi.erase(vi.begin() + 0);
    if(p[0] == nbvoxels_[0])
      vi.erase(vi.begin() + 2);
    if(p[1] == 0)
      vj.erase(vj.begin() + 0);
    if(p[1] == nbvoxels_[1])
      vj.erase(vj.begin() + 2);

    for(SimplexId const i : vi) {
      for(SimplexId const j : vj) {
        if(i != 0 or j != 0) {
          result.push_back(vertexId + i + j * vshift_[Di_]);
        }
      }
    }
  }
  return result;
}

template <size_t card>
int MultiresTriangulation<card>::getVertexBoundaryIndex(
  const SimplexId vertexId) const {

  if(dimensionality_ == 3) {
    std::array<SimplexId, 3> p{};
    vertexToPosition(vertexId, p);

    if(0 < p[0] and p[0] < nbvoxels_[0]) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 0;
        else if(p[2] == 0)
          return 1;
        else
          return 2;
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 3;
        else if(p[2] == 0)
          return 4;
        else
          return 5;
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 6;
        else if(p[2] == 0)
          return 7;
        else
          return 8;
      }
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 9;
        else if(p[2] == 0)
          return 10;
        else
          return 11;
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 12;
        else if(p[2] == 0)
          return 13;
        else
          return 14;
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 15;
        else if(p[2] == 0)
          return 16;
        else
          return 17;
      }
    } else {
      if(0 < p[1] and p[1] < nbvoxels_[1]) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 18;
        else if(p[2] == 0)
          return 19;
        else
          return 20;
      } else if(p[1] == 0) {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 21;
        else if(p[2] == 0)
          return 22;
        else
          return 23;
      } else {
        if(0 < p[2] and p[2] < nbvoxels_[2])
          return 24;
        else if(p[2] == 0)
          return 25;
        else
          return 26;
      }
    }

  } else if(dimensionality_ == 2) {
    std::array<SimplexId, 3> p{};
    vertexToPosition2d(vertexId, p);

    if(0 < p[0] and p[0] < nbvoxels_[Di_]) {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_]) {
        return 0;
      } else if(p[1] == 0) {
        return 1;
      } else {
        return 2;
      }
    } else if(p[0] == 0) {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_]) {
        return 3;
      } else if(p[1] == 0) {
        return 4;
      } else {
        return 5;
      }
    } else {
      if(0 < p[1] and p[1] < nbvoxels_[Dj_]) {
        return 6;
      } else if(p[1] == 0) {
        return 7;
      } else {
        return 8;
      }
    }
  }
  return -1;
}

template <size_t card>
void MultiresTriangulation<card>::findBoundaryRepresentatives(
  std::vector<SimplexId> &boundaryRepresentatives) {
  int const currentDecimationLevel = decimationLevel_;
  setDecimationLevel(0);
  if(dimensionality_ == 3) { // here gridDimension[i] > 1 for all i =0,1,2
    boundaryRepresentatives.resize(27, -1);
    // ADD THE CORNERS
    boundaryRepresentatives[13] = 0; // 0 0 0
    boundaryRepresentatives[22] = gridDimensions_[0] - 1; // Nx 0 0
    boundaryRepresentatives[16]
      = vshift_[0] * (gridDimensions_[1] - 1); // 0 Ny 0
    boundaryRepresentatives[14]
      = vshift_[1] * (gridDimensions_[2] - 1); // 0 0 Nz
    boundaryRepresentatives[25]
      = gridDimensions_[0] - 1
        + vshift_[0] * (gridDimensions_[1] - 1); // Nx Ny 0
    boundaryRepresentatives[23]
      = gridDimensions_[0] - 1
        + vshift_[1] * (gridDimensions_[2] - 1); // Nx 0 Nz
    boundaryRepresentatives[17]
      = vshift_[0] * (gridDimensions_[1] - 1)
        + vshift_[1] * (gridDimensions_[2] - 1); // 0 Ny Nz
    boundaryRepresentatives[26]
      = gridDimensions_[0] - 1 + vshift_[0] * (gridDimensions_[1] - 1)
        + vshift_[1] * (gridDimensions_[2] - 1); // Nx Ny Nz

    if(gridDimensions_[0] > 2) {
      boundaryRepresentatives[4] = 1; // 1 0 0
      boundaryRepresentatives[7]
        = 1 + vshift_[0] * (gridDimensions_[1] - 1); // 1 Ny 0
      boundaryRepresentatives[5]
        = 1 + vshift_[1] * (gridDimensions_[2] - 1); // 1 0 Nz
      boundaryRepresentatives[8]
        = 1 + vshift_[0] * (gridDimensions_[1] - 1)
          + vshift_[1] * (gridDimensions_[2] - 1); // 1 Ny Nz
    }
    if(gridDimensions_[1] > 2) {
      boundaryRepresentatives[10] = vshift_[0]; // 0 1 0
      boundaryRepresentatives[19]
        = gridDimensions_[0] - 1 + vshift_[0]; // Nx 1 0
      boundaryRepresentatives[11]
        = vshift_[0] + vshift_[1] * (gridDimensions_[2] - 1); // 0 1 Nz
      boundaryRepresentatives[20]
        = gridDimensions_[0] - 1 + vshift_[0]
          + vshift_[1] * (gridDimensions_[2] - 1); // Nx 1 Nz
    }
    if(gridDimensions_[2] > 2) {
      boundaryRepresentatives[12] = vshift_[1]; // 0 0 1
      boundaryRepresentatives[15]
        = vshift_[0] * (gridDimensions_[1] - 1) + vshift_[1]; // 0 Ny 1
      boundaryRepresentatives[21]
        = (gridDimensions_[0] - 1) + vshift_[1]; // Nx 0 1
      boundaryRepresentatives[24] = (gridDimensions_[0] - 1)
                                    + vshift_[0] * (gridDimensions_[1] - 1)
                                    + vshift_[1]; // Nx Ny 1
    }
    if(gridDimensions_[0] > 2 and gridDimensions_[1] > 2) {
      boundaryRepresentatives[1] = 1 + vshift_[0]; // 1 1 0
      boundaryRepresentatives[2]
        = 1 + vshift_[0] + vshift_[1] * (gridDimensions_[2] - 1); // 1 1 Nz
    }
    if(gridDimensions_[0] > 2 and gridDimensions_[2] > 2) {
      boundaryRepresentatives[3] = 1 + vshift_[1]; // 1 0 1
      boundaryRepresentatives[6]
        = 1 + vshift_[0] * (gridDimensions_[1] - 1) + vshift_[1]; // 1 Ny 1
    }
    if(gridDimensions_[1] > 2 and gridDimensions_[2] > 2) {
      boundaryRepresentatives[9] = vshift_[0] + vshift_[1]; // 0 1 1
      boundaryRepresentatives[18]
        = (gridDimensions_[0] - 1) + vshift_[0] + vshift_[1]; // Nx 1 1
    }
    if(gridDimensions_[0] > 2 and gridDimensions_[1] > 2
       and gridDimensions_[2] > 2) {
      boundaryRepresentatives[0] = 1 + vshift_[0] + vshift_[1];
    }
  } else if(dimensionality_ == 2) {
    boundaryRepresentatives.resize(9, -1);
    boundaryRepresentatives[4] = 0; // 0 0
    boundaryRepresentatives[5] = vshift_[0] * (gridDimensions_[1] - 1); // 0 Ny
    boundaryRepresentatives[7] = gridDimensions_[0] - 1; // Nx 0
    boundaryRepresentatives[8]
      = (gridDimensions_[0] - 1)
        + vshift_[0] * (gridDimensions_[1] - 1); // Nx Ny

    boundaryRepresentatives[1] = 1; // 1 0
    boundaryRepresentatives[2]
      = 1 + vshift_[0] * (gridDimensions_[1] - 1); // 1 Ny
    boundaryRepresentatives[3] = vshift_[0]; // 0 1
    boundaryRepresentatives[6] = (gridDimensions_[0] - 1) + vshift_[0]; // Nx 1

    boundaryRepresentatives[0] = 1 + vshift_[0]; // 1 1
  }
  setDecimationLevel(currentDecimationLevel);
}

template<size_t card>
bool ttk::MultiresTriangulation<card>::isBoundaryImpacted(SimplexId v) const {
  bool ret = false;
  if(dimensionality_ == 3) {
    std::array<SimplexId, 3> p{};
    vertexToPosition(v, p);
    if((nbvoxels_[0] % decimation_) and (p[0] + decimation_ > nbvoxels_[0])) {
      ret = true;
    }
    if((nbvoxels_[1] % decimation_) and (p[1] + decimation_ > nbvoxels_[1])) {
      ret = true;
    }
    if((nbvoxels_[2] % decimation_) and (p[2] + decimation_ > nbvoxels_[2])) {
      ret = true;
    }
  }
  return ret;
}

template<size_t card>
void ttk::MultiresTriangulation<card>::computeCoarsestDecimationLevel() {
  int maxDim = std::max(
    gridDimensions_[0], std::max(gridDimensions_[1], gridDimensions_[2]));
  int dl = 0;
  while(maxDim > 2) {
    maxDim
      = ((maxDim - 1) % 2) ? ((maxDim - 1) / 2 + 2) : ((maxDim - 1) / 2 + 1);
    dl++;
  }
  coarsestDL_ = dl;
}

template<size_t card>
int ttk::MultiresTriangulation<card>::RL_to_DL(int rl) {
  if(rl < 0) {
    return 0;
  } else {
    return std::max(coarsestDL_ - rl, 0);
  }
}

template<size_t card>
int ttk::MultiresTriangulation<card>::DL_to_RL(int dl) {
  return coarsestDL_ - dl;
}

// explicit instanciations
template class ttk::MultiresTriangulation<0>;
template class ttk::MultiresTriangulation<1>;
template class ttk::MultiresTriangulation<2>;
template class ttk::MultiresTriangulation<3>;
