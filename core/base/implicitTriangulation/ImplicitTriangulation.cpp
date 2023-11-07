#include <ImplicitTriangulation.h>


#include <numeric>

using namespace std;
using namespace ttk;

#define CASE_EDGE_POSITION_L_3D \
  case CRTPClass::EdgePosition::L_xnn_3D:  \
  case CRTPClass::EdgePosition::L_xn0_3D:  \
  case CRTPClass::EdgePosition::L_xnN_3D:  \
  case CRTPClass::EdgePosition::L_x0n_3D:  \
  case CRTPClass::EdgePosition::L_x00_3D:  \
  case CRTPClass::EdgePosition::L_x0N_3D:  \
  case CRTPClass::EdgePosition::L_xNn_3D:  \
  case CRTPClass::EdgePosition::L_xN0_3D:  \
  case CRTPClass::EdgePosition::L_xNN_3D
#define CASE_EDGE_POSITION_H_3D \
  case CRTPClass::EdgePosition::H_nyn_3D:  \
  case CRTPClass::EdgePosition::H_ny0_3D:  \
  case CRTPClass::EdgePosition::H_nyN_3D:  \
  case CRTPClass::EdgePosition::H_0yn_3D:  \
  case CRTPClass::EdgePosition::H_0y0_3D:  \
  case CRTPClass::EdgePosition::H_0yN_3D:  \
  case CRTPClass::EdgePosition::H_Nyn_3D:  \
  case CRTPClass::EdgePosition::H_Ny0_3D:  \
  case CRTPClass::EdgePosition::H_NyN_3D
#define CASE_EDGE_POSITION_P_3D \
  case CRTPClass::EdgePosition::P_nnz_3D:  \
  case CRTPClass::EdgePosition::P_n0z_3D:  \
  case CRTPClass::EdgePosition::P_nNz_3D:  \
  case CRTPClass::EdgePosition::P_0nz_3D:  \
  case CRTPClass::EdgePosition::P_00z_3D:  \
  case CRTPClass::EdgePosition::P_0Nz_3D:  \
  case CRTPClass::EdgePosition::P_Nnz_3D:  \
  case CRTPClass::EdgePosition::P_N0z_3D:  \
  case CRTPClass::EdgePosition::P_NNz_3D
#define CASE_EDGE_POSITION_D1_3D \
  case CRTPClass::EdgePosition::D1_xyn_3D:  \
  case CRTPClass::EdgePosition::D1_xy0_3D:  \
  case CRTPClass::EdgePosition::D1_xyN_3D
#define CASE_EDGE_POSITION_D2_3D \
  case CRTPClass::EdgePosition::D2_nyz_3D:  \
  case CRTPClass::EdgePosition::D2_0yz_3D:  \
  case CRTPClass::EdgePosition::D2_Nyz_3D
#define CASE_EDGE_POSITION_D3_3D \
  case CRTPClass::EdgePosition::D3_xnz_3D:  \
  case CRTPClass::EdgePosition::D3_x0z_3D:  \
  case CRTPClass::EdgePosition::D3_xNz_3D
#define CASE_EDGE_POSITION_L_2D \
  case CRTPClass::EdgePosition::L_xn_2D:   \
  case CRTPClass::EdgePosition::L_x0_2D:   \
  case CRTPClass::EdgePosition::L_xN_2D
#define CASE_EDGE_POSITION_H_2D \
  case CRTPClass::EdgePosition::H_ny_2D:   \
  case CRTPClass::EdgePosition::H_0y_2D:   \
  case CRTPClass::EdgePosition::H_Ny_2D

template<size_t card>
ImplicitTriangulation<card>::ImplicitTriangulation()
  : cellNumber_{}, vertexNumber_{}, edgeNumber_{}, triangleNumber_{},
    tetrahedronNumber_{}, isAccelerated_{} {
  this->setDebugMsgPrefix("ImplicitTriangulation");
#ifdef TTK_ENABLE_MPI
  this->hasPreconditionedDistributedEdges_ = true;
  this->hasPreconditionedDistributedTriangles_ = true;
#endif // TTK_ENABLE_MPI
}

template<size_t card>
ImplicitTriangulation<card>::~ImplicitTriangulation() = default;

template<size_t card>
int ImplicitTriangulation<card>::setInputGrid(const float &xOrigin,
                                        const float &yOrigin,
                                        const float &zOrigin,
                                        const float &xSpacing,
                                        const float &ySpacing,
                                        const float &zSpacing,
                                        const SimplexId &xDim,
                                        const SimplexId &yDim,
                                        const SimplexId &zDim) {

  // Dimensionality //
  if(xDim < 1 or yDim < 1 or zDim < 1)
    this->dimensionality_ = -1;
  else if(xDim > 1 and yDim > 1 and zDim > 1)
    this->dimensionality_ = 3;
  else if((xDim > 1 and yDim > 1) or (yDim > 1 and zDim > 1)
          or (xDim > 1 and zDim > 1))
    this->dimensionality_ = 2;
  else if(xDim > 1 or yDim > 1 or zDim > 1)
    this->dimensionality_ = 1;
  else
    this->dimensionality_ = 0;

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
  this->nbvoxels_[0] = xDim - 1;
  this->nbvoxels_[1] = yDim - 1;
  this->nbvoxels_[2] = zDim - 1;

  if(this->dimensionality_ == 3) {
    // VertexShift
    this->vshift_[0] = xDim;
    this->vshift_[1] = xDim * yDim;
    // EdgeSetDimensions
    esetdims_[0] = (xDim - 1) * yDim * zDim;
    esetdims_[1] = xDim * (yDim - 1) * zDim;
    esetdims_[2] = xDim * yDim * (zDim - 1);
    esetdims_[3] = (xDim - 1) * (yDim - 1) * zDim;
    esetdims_[4] = xDim * (yDim - 1) * (zDim - 1);
    esetdims_[5] = (xDim - 1) * yDim * (zDim - 1);
    esetdims_[6] = (xDim - 1) * (yDim - 1) * (zDim - 1);
    // EdgeSetShift
    this->esetshift_[0] = esetdims_[0];
    for(int k = 1; k < 7; ++k)
      this->esetshift_[k] = esetshift_[k - 1] + esetdims_[k];
    // EdgeShift
    this->eshift_[0] = xDim - 1;
    this->eshift_[1] = (xDim - 1) * yDim;
    this->eshift_[2] = xDim;
    this->eshift_[3] = xDim * (yDim - 1);
    this->eshift_[4] = xDim;
    this->eshift_[5] = xDim * yDim;
    this->eshift_[6] = xDim - 1;
    this->eshift_[7] = (xDim - 1) * (yDim - 1);
    this->eshift_[8] = xDim;
    this->eshift_[9] = xDim * (yDim - 1);
    this->eshift_[10] = xDim - 1;
    this->eshift_[11] = (xDim - 1) * yDim;
    this->eshift_[12] = xDim - 1;
    this->eshift_[13] = (xDim - 1) * (yDim - 1);
    // TriangleSetDimensions
    tsetdims_[0] = (xDim - 1) * (yDim - 1) * zDim * 2;
    tsetdims_[1] = (xDim - 1) * yDim * (zDim - 1) * 2;
    tsetdims_[2] = xDim * (yDim - 1) * (zDim - 1) * 2;
    tsetdims_[3] = (xDim - 1) * (yDim - 1) * (zDim - 1) * 2;
    tsetdims_[4] = (xDim - 1) * (yDim - 1) * (zDim - 1) * 2;
    tsetdims_[5] = (xDim - 1) * (yDim - 1) * (zDim - 1) * 2;
    // TriangleSetShift
    this->tsetshift_[0] = tsetdims_[0];
    for(int k = 1; k < 6; ++k)
      this->tsetshift_[k] = this->tsetshift_[k - 1] + tsetdims_[k];
    // TriangleShift
    this->tshift_[0] = (xDim - 1) * 2;
    this->tshift_[1] = (xDim - 1) * (yDim - 1) * 2;
    this->tshift_[2] = (xDim - 1) * 2;
    this->tshift_[3] = (xDim - 1) * yDim * 2;
    this->tshift_[4] = xDim * 2;
    this->tshift_[5] = xDim * (yDim - 1) * 2;
    this->tshift_[6] = (xDim - 1) * 2;
    this->tshift_[7] = (xDim - 1) * (yDim - 1) * 2;
    this->tshift_[8] = (xDim - 1) * 2;
    this->tshift_[9] = (xDim - 1) * (yDim - 1) * 2;
    this->tshift_[10] = (xDim - 1) * 2;
    this->tshift_[11] = (xDim - 1) * (yDim - 1) * 2;
    // TetrahedronShift
    this->tetshift_[0] = (xDim - 1) * 6;
    this->tetshift_[1] = (xDim - 1) * (yDim - 1) * 6;

    // Numbers
    vertexNumber_ = xDim * yDim * zDim;
    this->edgeNumber_ = 0;
    for(int k = 0; k < 7; ++k)
      this->edgeNumber_ += esetdims_[k];
    triangleNumber_ = 0;
    for(int k = 0; k < 6; ++k)
      triangleNumber_ += tsetdims_[k];
    tetrahedronNumber_ = (xDim - 1) * (yDim - 1) * (zDim - 1) * 6;
    cellNumber_ = tetrahedronNumber_;

    checkAcceleration();
  } else if(this->dimensionality_ == 2) {
    // dimensions selectors
    if(xDim == 1) {
      this->Di_ = 1;
      this->Dj_ = 2;
    } else if(yDim == 1) {
      this->Di_ = 0;
      this->Dj_ = 2;
    } else {
      this->Di_ = 0;
      this->Dj_ = 1;
    }
    // VertexShift
    this->vshift_[0] = this->dimensions_[this->Di_];
    // EdgeSetDimensions
    esetdims_[0] = (this->dimensions_[this->Di_] - 1) * this->dimensions_[this->Dj_];
    esetdims_[1] = this->dimensions_[this->Di_] * (this->dimensions_[this->Dj_] - 1);
    esetdims_[2] = (this->dimensions_[this->Di_] - 1) * (this->dimensions_[this->Dj_] - 1);
    // EdgeSetShift
    this->esetshift_[0] = esetdims_[0];
    for(int k = 1; k < 3; ++k)
      this->esetshift_[k] = esetshift_[k - 1] + esetdims_[k];
    // EdgeShift
    this->eshift_[0] = this->dimensions_[this->Di_] - 1;
    this->eshift_[2] = this->dimensions_[this->Di_];
    this->eshift_[4] = this->dimensions_[this->Di_] - 1;
    // TriangleShift
    this->tshift_[0] = (this->dimensions_[this->Di_] - 1) * 2;

    // Numbers
    vertexNumber_ = this->dimensions_[this->Di_] * this->dimensions_[this->Dj_];
    this->edgeNumber_ = 0;
    for(int k = 0; k < 3; ++k)
      this->edgeNumber_ += esetdims_[k];
    triangleNumber_ = (this->dimensions_[this->Di_] - 1) * (this->dimensions_[this->Dj_] - 1) * 2;
    cellNumber_ = triangleNumber_;

    checkAcceleration();
  } else if(this->dimensionality_ == 1) {
    // dimensions selectors
    for(int k = 0; k < 3; ++k) {
      if(this->dimensions_[k] > 1) {
        this->Di_ = k;
        break;
      }
    }

    // Numbers
    vertexNumber_ = this->dimensions_[this->Di_];
    this->edgeNumber_ = vertexNumber_ - 1;
    cellNumber_ = this->edgeNumber_;
  }

  return 0;
}

template<size_t card>
int ImplicitTriangulation<card>::checkAcceleration() {
  isAccelerated_ = false;

  unsigned long long int msb[3];
  if(this->dimensionality_ == 3) {
    bool allDimensionsArePowerOfTwo = true;
    for(int k = 0; k < 3; ++k)
      if(!isPowerOfTwo(this->dimensions_[k], msb[k]))
        allDimensionsArePowerOfTwo = false;

    if(allDimensionsArePowerOfTwo) {
      mod_[0] = this->dimensions_[0] - 1;
      mod_[1] = this->dimensions_[0] * this->dimensions_[1] - 1;
      this->div_[0] = msb[0];
      this->div_[1] = msb[0] + msb[1];
      isAccelerated_ = true;
    }
  } else if(this->dimensionality_ == 2) {
    bool const isDi = isPowerOfTwo(this->dimensions_[this->Di_], msb[this->Di_]);
    bool const isDj = isPowerOfTwo(this->dimensions_[this->Dj_], msb[this->Dj_]);
    bool const allDimensionsArePowerOfTwo = (isDi and isDj);

    if(allDimensionsArePowerOfTwo) {
      mod_[0] = this->dimensions_[this->Di_] - 1;
      this->div_[0] = msb[this->Di_];
      isAccelerated_ = true;
    }
  }

  if(isAccelerated_) {
    this->printMsg("Accelerated getVertex*() requests.", debug::Priority::INFO);
  }

  return 0;
}

template<size_t card>
bool ImplicitTriangulation<card>::isPowerOfTwo(unsigned long long int v,
                                         unsigned long long int &r) {
  if(v && !(v & (v - 1))) {
    r = 0;
    while(v >>= 1)
      r++;
    return true;
  }
  return false;
}

template <size_t card, typename Derived>
bool ImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  isVertexOnBoundary)(const SimplexId &vertexId) const {

  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifdef TTK_ENABLE_MPI
  if(this->metaGrid_ != nullptr) {
    return this->isVertexOnGlobalBoundaryInternal(vertexId);
  }
#endif // TTK_ENABLE_MPI

#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return false;
#endif // !TTK_ENABLE_KAMIKAZE

  switch(this->underlying().getVertexPosition(vertexId)) {
    case CRTPClass::VertexPosition::CENTER_3D:
    case CRTPClass::VertexPosition::CENTER_2D:
    case CRTPClass::VertexPosition::CENTER_1D:
      return false;
    default:
      return true;
  }
}

template <size_t card, typename Derived>
bool ImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  isEdgeOnBoundary)(const SimplexId &edgeId) const {

  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifdef TTK_ENABLE_MPI
  if(this->metaGrid_ != nullptr) {
    return this->isEdgeOnGlobalBoundaryInternal(edgeId);
  }
#endif // TTK_ENABLE_MPI

#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= this->edgeNumber_)
    return false;
#endif // !TTK_ENABLE_KAMIKAZE

  switch(this->underlying().getEdgePosition(edgeId)) {
    case CRTPClass::EdgePosition::L_xnn_3D:
    case CRTPClass::EdgePosition::H_nyn_3D:
    case CRTPClass::EdgePosition::P_nnz_3D:
    case CRTPClass::EdgePosition::D1_xyn_3D:
    case CRTPClass::EdgePosition::D2_nyz_3D:
    case CRTPClass::EdgePosition::D3_xnz_3D:
    case CRTPClass::EdgePosition::D4_3D:
    case CRTPClass::EdgePosition::L_xn_2D:
    case CRTPClass::EdgePosition::H_ny_2D:
    case CRTPClass::EdgePosition::D1_2D:
      return false;
    default:
      break;
  }
  return true;
}

template<size_t card>
bool ImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(isTriangleOnBoundary)(
  const SimplexId &triangleId) const {

#ifdef TTK_ENABLE_MPI
  if(this->metaGrid_ != nullptr) {
    return this->isTriangleOnGlobalBoundaryInternal(triangleId);
  }
#endif // TTK_ENABLE_MPI

#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return false;
#endif // !TTK_ENABLE_KAMIKAZE

  if(this->dimensionality_ == 3)
    return (TTK_TRIANGULATION_INTERNAL(this->getTriangleStarNumber)(triangleId) == 1);

  return false;
}

template <size_t card, typename Derived>
int ImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  getVertexNeighbor)(const SimplexId &vertexId,
                     const int &localNeighborId,
                     SimplexId &neighborId) const {

  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(localNeighborId < 0
     or localNeighborId >= this->getVertexNeighborNumber(vertexId))
    return -1;
#endif // !TTK_ENABLE_KAMIKAZE

  switch(this->underlying().getVertexPosition(vertexId)) {
    case CRTPClass::VertexPosition::CENTER_3D:
      neighborId = vertexId + this->vertexNeighborABCDEFGH_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::FRONT_FACE_3D:
      neighborId = vertexId + this->vertexNeighborABCD_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::BACK_FACE_3D:
      neighborId = vertexId + this->vertexNeighborEFGH_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::TOP_FACE_3D:
      neighborId = vertexId + this->vertexNeighborAEFB_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::BOTTOM_FACE_3D:
      neighborId = vertexId + this->vertexNeighborGHDC_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::LEFT_FACE_3D:
      neighborId = vertexId + this->vertexNeighborAEGC_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::RIGHT_FACE_3D:
      neighborId = vertexId + this->vertexNeighborBFHD_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::TOP_FRONT_EDGE_3D: // ab
      neighborId = vertexId + this->vertexNeighborAB_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
      neighborId = vertexId + this->vertexNeighborCD_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::LEFT_FRONT_EDGE_3D: // ac
      neighborId = vertexId + this->vertexNeighborAC_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
      neighborId = vertexId + this->vertexNeighborBD_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::TOP_BACK_EDGE_3D: // ef
      neighborId = vertexId + this->vertexNeighborEF_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
      neighborId = vertexId + this->vertexNeighborGH_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::LEFT_BACK_EDGE_3D: // eg
      neighborId = vertexId + this->vertexNeighborEG_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::RIGHT_BACK_EDGE_3D: // fh
      neighborId = vertexId + this->vertexNeighborFH_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_EDGE_3D: // ae
      neighborId = vertexId + this->vertexNeighborAE_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_EDGE_3D: // bf
      neighborId = vertexId + this->vertexNeighborBF_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
      neighborId = vertexId + this->vertexNeighborCG_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
      neighborId = vertexId + this->vertexNeighborDH_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
      neighborId = vertexId + this->vertexNeighborA_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
      neighborId = vertexId + this->vertexNeighborB_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
      neighborId = vertexId + this->vertexNeighborC_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
      neighborId = vertexId + this->vertexNeighborD_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
      neighborId = vertexId + this->vertexNeighborE_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
      neighborId = vertexId + this->vertexNeighborF_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
      neighborId = vertexId + this->vertexNeighborG_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
      neighborId = vertexId + this->vertexNeighborH_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::CENTER_2D:
      neighborId = vertexId + this->vertexNeighbor2dABCD_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::TOP_EDGE_2D:
      neighborId = vertexId + this->vertexNeighbor2dAB_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::BOTTOM_EDGE_2D:
      neighborId = vertexId + this->vertexNeighbor2dCD_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::LEFT_EDGE_2D:
      neighborId = vertexId + this->vertexNeighbor2dAC_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::RIGHT_EDGE_2D:
      neighborId = vertexId + this->vertexNeighbor2dBD_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_CORNER_2D: // a
      neighborId = vertexId + this->vertexNeighbor2dA_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_CORNER_2D: // b
      neighborId = vertexId + this->vertexNeighbor2dB_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
      neighborId = vertexId + this->vertexNeighbor2dC_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
      neighborId = vertexId + this->vertexNeighbor2dD_[localNeighborId];
      break;
    case CRTPClass::VertexPosition::CENTER_1D:
      neighborId = (localNeighborId == 0 ? vertexId + 1 : vertexId - 1);
      break;
    case CRTPClass::VertexPosition::LEFT_CORNER_1D:
      neighborId = vertexId + 1;
      break;
    case CRTPClass::VertexPosition::RIGHT_CORNER_1D:
      neighborId = vertexId - 1;
      break;
    default:
      neighborId = -1;
      break;
  }

  return 0;
}

template<size_t card>
const vector<vector<SimplexId>> *
  ImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getVertexNeighbors)() {
  if(this->vertexNeighborList_.empty()) {
    Timer t;
    this->vertexNeighborList_.resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      this->vertexNeighborList_[i].resize(this->getVertexNeighborNumber(i));
      for(SimplexId j = 0; j < (SimplexId)this->vertexNeighborList_[i].size(); ++j)
        this->getVertexNeighbor(i, j, this->vertexNeighborList_[i][j]);
    }

    this->printMsg("Built " + to_string(vertexNumber_) + " vertex neighbors.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->vertexNeighborList_;
}

template<size_t card>
SimplexId ImplicitTriangulation<card>::getVertexEdgeNumberInternal(
  const SimplexId &vertexId) const {
  return TTK_TRIANGULATION_INTERNAL(this->getVertexNeighborNumber)(vertexId);
}

template <size_t card, typename Derived>
int ImplicitTriangulationCRTP<card, Derived>::getVertexEdgeInternal(
  const SimplexId &vertexId, const int &localEdgeId, SimplexId &edgeId) const {
  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(localEdgeId < 0 or localEdgeId >= this->getVertexEdgeNumberInternal(vertexId))
    return -1;
#endif
  //    e--------f
  //   /|       /|
  //  / |      / |
  // a--g-----b--h
  // | /      | /
  // |/       |/
  // c--------d
  //
  // Classement des "Edges" et dans cet ordre:
  // L: largeur (type ab)
  // H: hauteur (type ac)
  // P: profondeur (type ae)
  // D1: diagonale1 (type bc)
  // D2: diagonale2 (type ag)
  // D3: diagonale3 (type be)
  // D4: diagonale4 (type bg)

  const auto &p = this->underlying().getVertexCoords(vertexId);

  switch(this->underlying().getVertexPosition(vertexId)) {
    case CRTPClass::VertexPosition::CENTER_3D:
      edgeId = this->getVertexEdgeABCDEFGH(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::FRONT_FACE_3D:
      edgeId = this->getVertexEdgeABDC(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::BACK_FACE_3D:
      edgeId = this->getVertexEdgeEFHG(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::TOP_FACE_3D:
      edgeId = this->getVertexEdgeAEFB(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_FACE_3D:
      edgeId = this->getVertexEdgeGHDC(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::LEFT_FACE_3D:
      edgeId = this->getVertexEdgeAEGC(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::RIGHT_FACE_3D:
      edgeId = this->getVertexEdgeBFHD(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::TOP_FRONT_EDGE_3D: // ab
      edgeId = this->getVertexEdgeAB(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
      edgeId = this->getVertexEdgeCD(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::LEFT_FRONT_EDGE_3D: // ac
      edgeId = this->getVertexEdgeAC(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
      edgeId = this->getVertexEdgeBD(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::TOP_BACK_EDGE_3D: // ef
      edgeId = this->getVertexEdgeEF(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
      edgeId = this->getVertexEdgeGH(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::LEFT_BACK_EDGE_3D: // eg
      edgeId = this->getVertexEdgeEG(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::RIGHT_BACK_EDGE_3D: // fh
      edgeId = this->getVertexEdgeFH(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_EDGE_3D: // ae
      edgeId = this->getVertexEdgeAE(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_EDGE_3D: // bf
      edgeId = this->getVertexEdgeBF(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
      edgeId = this->getVertexEdgeCG(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
      edgeId = this->getVertexEdgeDH(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
      edgeId = this->getVertexEdgeA(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
      edgeId = this->getVertexEdgeB(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
      edgeId = this->getVertexEdgeC(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
      edgeId = this->getVertexEdgeD(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
      edgeId = this->getVertexEdgeE(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
      edgeId = this->getVertexEdgeF(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
      edgeId = this->getVertexEdgeG(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
      edgeId = this->getVertexEdgeH(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::CENTER_2D:
      edgeId = this->getVertexEdge2dABCD(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::TOP_EDGE_2D:
      edgeId = this->getVertexEdge2dAB(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_EDGE_2D:
      edgeId = this->getVertexEdge2dCD(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::LEFT_EDGE_2D:
      edgeId = this->getVertexEdge2dAC(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::RIGHT_EDGE_2D:
      edgeId = this->getVertexEdge2dBD(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_CORNER_2D: // a
      edgeId = this->getVertexEdge2dA(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_CORNER_2D: // b
      edgeId = this->getVertexEdge2dB(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
      edgeId = this->getVertexEdge2dC(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
      edgeId = this->getVertexEdge2dD(p.data(), localEdgeId);
      break;
    case CRTPClass::VertexPosition::CENTER_1D:
      edgeId = (localEdgeId == 0 ? vertexId : vertexId - 1);
      break;
    case CRTPClass::VertexPosition::LEFT_CORNER_1D:
      edgeId = vertexId;
      break;
    case CRTPClass::VertexPosition::RIGHT_CORNER_1D:
      edgeId = vertexId - 1;
      break;
    default:
      edgeId = -1;
  }

  return 0;
}

template<size_t card>
const vector<vector<SimplexId>> *
  ImplicitTriangulation<card>::getVertexEdgesInternal() {
  if(this->vertexEdgeList_.empty()) {
    Timer t;

    this->vertexEdgeList_.resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      this->vertexEdgeList_[i].resize(getVertexEdgeNumberInternal(i));
      for(SimplexId j = 0; j < (SimplexId)this->vertexEdgeList_[i].size(); ++j)
        this->getVertexEdgeInternal(i, j, this->vertexEdgeList_[i][j]);
    }

    this->printMsg("Built " + to_string(vertexNumber_) + " vertex edges.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->vertexEdgeList_;
}

template <size_t card, typename Derived>
SimplexId ImplicitTriangulationCRTP<card, Derived>::getVertexTriangleNumberInternal(
  const SimplexId &vertexId) const {
  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#endif

  switch(this->underlying().getVertexPosition(vertexId)) {
    case CRTPClass::VertexPosition::CENTER_3D:
      return 36;
    case CRTPClass::VertexPosition::FRONT_FACE_3D:
    case CRTPClass::VertexPosition::BACK_FACE_3D:
    case CRTPClass::VertexPosition::TOP_FACE_3D:
    case CRTPClass::VertexPosition::BOTTOM_FACE_3D:
    case CRTPClass::VertexPosition::LEFT_FACE_3D:
    case CRTPClass::VertexPosition::RIGHT_FACE_3D:
      return 21;
    case CRTPClass::VertexPosition::TOP_FRONT_EDGE_3D: // ab
    case CRTPClass::VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
    case CRTPClass::VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
    case CRTPClass::VertexPosition::LEFT_BACK_EDGE_3D: // eg
    case CRTPClass::VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
    case CRTPClass::VertexPosition::TOP_RIGHT_EDGE_3D: // bf
      return 15;
    case CRTPClass::VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
    case CRTPClass::VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
      return 12;
    case CRTPClass::VertexPosition::TOP_BACK_EDGE_3D: // ef
    case CRTPClass::VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
    case CRTPClass::VertexPosition::LEFT_FRONT_EDGE_3D: // ac
    case CRTPClass::VertexPosition::TOP_LEFT_EDGE_3D: // ae
    case CRTPClass::VertexPosition::RIGHT_BACK_EDGE_3D: // fh
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
      return 9;
    case CRTPClass::VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
    case CRTPClass::VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
    case CRTPClass::VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
    case CRTPClass::VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
      return 5;
    default: // 1D + 2D
      break;
  }

  return 0;
}

template <size_t card, typename Derived>
int ImplicitTriangulationCRTP<card, Derived>::getVertexTriangleInternal(
  const SimplexId &vertexId,
  const int &localTriangleId,
  SimplexId &triangleId) const {
  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(localTriangleId < 0
     or localTriangleId >= this->getVertexTriangleNumberInternal(vertexId))
    return -1;
#endif

  const auto &p = this->underlying().getVertexCoords(vertexId);

  switch(this->underlying().getVertexPosition(vertexId)) {
    case CRTPClass::VertexPosition::CENTER_3D:
      triangleId = this->getVertexTriangleABCDEFGH(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::FRONT_FACE_3D:
      triangleId = this->getVertexTriangleABDC(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::BACK_FACE_3D:
      triangleId = this->getVertexTriangleEFHG(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::TOP_FACE_3D:
      triangleId = this->getVertexTriangleAEFB(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_FACE_3D:
      triangleId = this->getVertexTriangleGHDC(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::LEFT_FACE_3D:
      triangleId = this->getVertexTriangleAEGC(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::RIGHT_FACE_3D:
      triangleId = this->getVertexTriangleBFHD(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::TOP_FRONT_EDGE_3D: // ab
      triangleId = this->getVertexTriangleAB(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
      triangleId = this->getVertexTriangleCD(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::LEFT_FRONT_EDGE_3D: // ac
      triangleId = this->getVertexTriangleAC(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
      triangleId = this->getVertexTriangleBD(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::TOP_BACK_EDGE_3D: // ef
      triangleId = this->getVertexTriangleEF(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
      triangleId = this->getVertexTriangleGH(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::LEFT_BACK_EDGE_3D: // eg
      triangleId = this->getVertexTriangleEG(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::RIGHT_BACK_EDGE_3D: // fh
      triangleId = this->getVertexTriangleFH(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_EDGE_3D: // ae
      triangleId = this->getVertexTriangleAE(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_EDGE_3D: // bf
      triangleId = this->getVertexTriangleBF(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
      triangleId = this->getVertexTriangleCG(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
      triangleId = this->getVertexTriangleDH(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
      triangleId = this->getVertexTriangleA(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
      triangleId = this->getVertexTriangleB(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
      triangleId = this->getVertexTriangleC(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
      triangleId = this->getVertexTriangleD(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
      triangleId = this->getVertexTriangleE(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
      triangleId = this->getVertexTriangleF(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
      triangleId = this->getVertexTriangleG(p.data(), localTriangleId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
      triangleId = this->getVertexTriangleH(p.data(), localTriangleId);
      break;
    default: // 1D + 2D
      triangleId = -1;
      break;
  }

  return 0;
}

template<size_t card>
const vector<vector<SimplexId>> *
  ImplicitTriangulation<card>::getVertexTrianglesInternal() {
  if(this->vertexTriangleList_.empty()) {
    Timer t;

    this->vertexTriangleList_.resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      this->vertexTriangleList_[i].resize(this->getVertexTriangleNumberInternal(i));
      for(SimplexId j = 0; j < (SimplexId)this->vertexTriangleList_[i].size(); ++j)
        this->getVertexTriangleInternal(i, j, this->vertexTriangleList_[i][j]);
    }

    this->printMsg("Built " + to_string(vertexNumber_) + " vertex triangles.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->vertexTriangleList_;
}

template<size_t card>
SimplexId ImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(
  getVertexLinkNumber)(const SimplexId &vertexId) const {
  return TTK_TRIANGULATION_INTERNAL(this->getVertexStarNumber)(vertexId);
}

template <size_t card, typename Derived>
int ImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  getVertexLink)(const SimplexId &vertexId,
                 const int &localLinkId,
                 SimplexId &linkId) const {

  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(localLinkId < 0 or localLinkId >= this->getVertexLinkNumber(vertexId))
    return -1;
#endif // !TTK_ENABLE_KAMIKAZE

  const auto &p = this->underlying().getVertexCoords(vertexId);

  switch(this->underlying().getVertexPosition(vertexId)) {
    case CRTPClass::VertexPosition::CENTER_3D:
      linkId = this->getVertexLinkABCDEFGH(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::FRONT_FACE_3D:
      linkId = this->getVertexLinkABDC(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::BACK_FACE_3D:
      linkId = this->getVertexLinkEFHG(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::TOP_FACE_3D:
      linkId = this->getVertexLinkAEFB(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_FACE_3D:
      linkId = this->getVertexLinkGHDC(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::LEFT_FACE_3D:
      linkId = this->getVertexLinkAEGC(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::RIGHT_FACE_3D:
      linkId = this->getVertexLinkBFHD(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::TOP_FRONT_EDGE_3D: // ab
      linkId = this->getVertexLinkAB(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
      linkId = this->getVertexLinkCD(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::LEFT_FRONT_EDGE_3D: // ac
      linkId = this->getVertexLinkAC(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
      linkId = this->getVertexLinkBD(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::TOP_BACK_EDGE_3D: // ef
      linkId = this->getVertexLinkEF(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
      linkId = this->getVertexLinkGH(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::LEFT_BACK_EDGE_3D: // eg
      linkId = this->getVertexLinkEG(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::RIGHT_BACK_EDGE_3D: // fh
      linkId = this->getVertexLinkFH(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_EDGE_3D: // ae
      linkId = this->getVertexLinkAE(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_EDGE_3D: // bf
      linkId = this->getVertexLinkBF(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
      linkId = this->getVertexLinkCG(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
      linkId = this->getVertexLinkDH(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
      linkId = this->getVertexLinkA(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
      linkId = this->getVertexLinkB(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
      linkId = this->getVertexLinkC(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
      linkId = this->getVertexLinkD(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
      linkId = this->getVertexLinkE(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
      linkId = this->getVertexLinkF(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
      linkId = this->getVertexLinkG(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
      linkId = this->getVertexLinkH(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::CENTER_2D:
      linkId = this->getVertexLink2dABCD(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::TOP_EDGE_2D:
      linkId = this->getVertexLink2dAB(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_EDGE_2D:
      linkId = this->getVertexLink2dCD(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::LEFT_EDGE_2D:
      linkId = this->getVertexLink2dAC(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::RIGHT_EDGE_2D:
      linkId = this->getVertexLink2dBD(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_CORNER_2D: // a
      linkId = this->getVertexLink2dA(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_CORNER_2D: // b
      linkId = this->getVertexLink2dB(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
      linkId = this->getVertexLink2dC(p.data(), localLinkId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
      linkId = this->getVertexLink2dD(p.data(), localLinkId);
      break;
    default: // 1D
      linkId = -1;
      break;
  };

  return 0;
}

template<size_t card>
const vector<vector<SimplexId>> *
  ImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getVertexLinks)() {
  if(this->vertexLinkList_.empty()) {
    Timer t;

    this->vertexLinkList_.resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      this->vertexLinkList_[i].resize(getVertexLinkNumber(i));
      for(SimplexId j = 0; j < (SimplexId)this->vertexLinkList_[i].size(); ++j)
        this->getVertexLink(i, j, this->vertexLinkList_[i][j]);
    }

    this->printMsg("Built " + to_string(vertexNumber_) + " vertex links.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->vertexLinkList_;
}

template <size_t card, typename Derived>
SimplexId ImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  getVertexStarNumber)(const SimplexId &vertexId) const {

  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(vertexId < 0 or vertexId >= vertexNumber_)
    return -1;
#endif // !TTK_ENABLE_KAMIKAZE

  switch(this->underlying().getVertexPosition(vertexId)) {
    case CRTPClass::VertexPosition::CENTER_3D:
      return 24;
    case CRTPClass::VertexPosition::FRONT_FACE_3D:
    case CRTPClass::VertexPosition::BACK_FACE_3D:
    case CRTPClass::VertexPosition::TOP_FACE_3D:
    case CRTPClass::VertexPosition::BOTTOM_FACE_3D:
    case CRTPClass::VertexPosition::LEFT_FACE_3D:
    case CRTPClass::VertexPosition::RIGHT_FACE_3D:
      return 12;
    case CRTPClass::VertexPosition::TOP_FRONT_EDGE_3D: // ab
    case CRTPClass::VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
    case CRTPClass::VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
    case CRTPClass::VertexPosition::LEFT_BACK_EDGE_3D: // eg
    case CRTPClass::VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
    case CRTPClass::VertexPosition::TOP_RIGHT_EDGE_3D: // bf
      return 8;
    case CRTPClass::VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
    case CRTPClass::VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
    case CRTPClass::VertexPosition::CENTER_2D:
      return 6;
    case CRTPClass::VertexPosition::TOP_BACK_EDGE_3D: // ef
    case CRTPClass::VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
    case CRTPClass::VertexPosition::LEFT_FRONT_EDGE_3D: // ac
    case CRTPClass::VertexPosition::TOP_LEFT_EDGE_3D: // ae
    case CRTPClass::VertexPosition::RIGHT_BACK_EDGE_3D: // fh
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
      return 4;
    case CRTPClass::VertexPosition::TOP_EDGE_2D: // ab
    case CRTPClass::VertexPosition::BOTTOM_EDGE_2D: // cd
    case CRTPClass::VertexPosition::LEFT_EDGE_2D: // ac
    case CRTPClass::VertexPosition::RIGHT_EDGE_2D: // bd
      return 3;
    case CRTPClass::VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
    case CRTPClass::VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
    case CRTPClass::VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
    case CRTPClass::VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
    case CRTPClass::VertexPosition::TOP_RIGHT_CORNER_2D: // b
    case CRTPClass::VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
      return 2;
    case CRTPClass::VertexPosition::TOP_LEFT_CORNER_2D: // a
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
      return 1;
    default: // 1D
      break;
  }

  return 0;
}

template <size_t card, typename Derived>
int ImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  getVertexStar)(const SimplexId &vertexId,
                 const int &localStarId,
                 SimplexId &starId) const {

  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(localStarId < 0 or localStarId >= this->getVertexStarNumber(vertexId))
    return -1;
#endif // !TTK_ENABLE_KAMIKAZE

  const auto &p = this->underlying().getVertexCoords(vertexId);

  switch(this->underlying().getVertexPosition(vertexId)) {
    case CRTPClass::VertexPosition::CENTER_3D:
      starId = this->getVertexStarABCDEFGH(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::FRONT_FACE_3D:
      starId = this->getVertexStarABDC(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::BACK_FACE_3D:
      starId = this->getVertexStarEFHG(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::TOP_FACE_3D:
      starId = this->getVertexStarAEFB(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_FACE_3D:
      starId = this->getVertexStarGHDC(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::LEFT_FACE_3D:
      starId = this->getVertexStarAEGC(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::RIGHT_FACE_3D:
      starId = this->getVertexStarBFHD(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::TOP_FRONT_EDGE_3D: // ab
      starId = this->getVertexStarAB(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
      starId = this->getVertexStarCD(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::LEFT_FRONT_EDGE_3D: // ac
      starId = this->getVertexStarAC(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
      starId = this->getVertexStarBD(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::TOP_BACK_EDGE_3D: // ef
      starId = this->getVertexStarEF(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
      starId = this->getVertexStarGH(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::LEFT_BACK_EDGE_3D: // eg
      starId = this->getVertexStarEG(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::RIGHT_BACK_EDGE_3D: // fh
      starId = this->getVertexStarFH(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_EDGE_3D: // ae
      starId = this->getVertexStarAE(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_EDGE_3D: // bf
      starId = this->getVertexStarBF(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
      starId = this->getVertexStarCG(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
      starId = this->getVertexStarDH(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
      starId = this->getVertexStarA(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
      starId = this->getVertexStarB(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
      starId = this->getVertexStarC(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
      starId = this->getVertexStarD(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
      starId = this->getVertexStarE(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
      starId = this->getVertexStarF(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
      starId = this->getVertexStarG(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
      starId = this->getVertexStarH(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::CENTER_2D:
      starId = this->getVertexStar2dABCD(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::TOP_EDGE_2D:
      starId = this->getVertexStar2dAB(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_EDGE_2D:
      starId = this->getVertexStar2dCD(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::LEFT_EDGE_2D:
      starId = this->getVertexStar2dAC(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::RIGHT_EDGE_2D:
      starId = this->getVertexStar2dBD(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::TOP_LEFT_CORNER_2D: // a
      starId = this->getVertexStar2dA(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::TOP_RIGHT_CORNER_2D: // b
      starId = this->getVertexStar2dB(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
      starId = this->getVertexStar2dC(p.data(), localStarId);
      break;
    case CRTPClass::VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
      starId = this->getVertexStar2dD(p.data(), localStarId);
      break;
    default: // 1D
      starId = -1;
      break;
  }

  return 0;
}

template<size_t card>
const vector<vector<SimplexId>> *
  ImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getVertexStars)() {

  if(this->vertexStarList_.empty()) {
    Timer t;
    this->vertexStarList_.resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; ++i) {
      this->vertexStarList_[i].resize(this->getVertexStarNumber(i));
      for(SimplexId j = 0; j < (SimplexId)this->vertexStarList_[i].size(); ++j)
        this->getVertexStar(i, j, this->vertexStarList_[i][j]);
    }

    this->printMsg("Built " + to_string(vertexNumber_) + " vertex stars.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->vertexStarList_;
}

template <size_t card, typename Derived>
int ImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  getVertexPoint)(const SimplexId &vertexId,
                  float &x,
                  float &y,
                  float &z) const {

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
int ImplicitTriangulationCRTP<card, Derived>::getEdgeVertexInternal(
  const SimplexId &edgeId,
  const int &localVertexId,
  SimplexId &vertexId) const {
  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= this->edgeNumber_)
    return -1;
  if(localVertexId < 0 or localVertexId >= 2)
    return -2;
#endif

  const auto &p = this->underlying().getEdgeCoords(edgeId);

  const auto helper3d = [&](const SimplexId a, const SimplexId b) -> SimplexId {
    if(this->isAccelerated_) {
      const auto tmp = p[0] + (p[1] << this->div_[0]) + (p[2] << this->div_[1]);
      return (localVertexId == 0) ? tmp + a : tmp + b;
    } else {
      const auto tmp = p[0] + (p[1] * this->vshift_[0]) + (p[2] * this->vshift_[1]);
      return (localVertexId == 0) ? tmp + a : tmp + b;
    }
  };

  const auto helper2d = [&](const SimplexId a, const SimplexId b) -> SimplexId {
    if(this->isAccelerated_) {
      const auto tmp = p[0] + (p[1] << this->div_[0]);
      return localVertexId == 0 ? tmp + a : tmp + b;
    } else {
      const auto tmp = p[0] + (p[1] * this->vshift_[0]);
      return localVertexId == 0 ? tmp + a : tmp + b;
    }
  };

  switch(this->underlying().getEdgePosition(edgeId)) {
  CASE_EDGE_POSITION_L_3D:
    vertexId = helper3d(0, 1);
    break;
  CASE_EDGE_POSITION_H_3D:
    vertexId = helper3d(0, this->vshift_[0]);
    break;
  CASE_EDGE_POSITION_P_3D:
    vertexId = helper3d(0, this->vshift_[1]);
    break;
  CASE_EDGE_POSITION_D1_3D:
    vertexId = helper3d(1, this->vshift_[0]);
    break;
  CASE_EDGE_POSITION_D2_3D:
    vertexId = helper3d(0, this->vshift_[0] + this->vshift_[1]);
    break;
  CASE_EDGE_POSITION_D3_3D:
    vertexId = helper3d(1, this->vshift_[1]);
    break;
    case CRTPClass::EdgePosition::D4_3D:
      vertexId = helper3d(1, this->vshift_[0] + this->vshift_[1]);
      break;

    CASE_EDGE_POSITION_L_2D:
      vertexId = helper2d(0, 1);
      break;
    CASE_EDGE_POSITION_H_2D:
      vertexId = helper2d(0, this->vshift_[0]);
      break;
    case CRTPClass::EdgePosition::D1_2D:
      vertexId = helper2d(1, this->vshift_[0]);
      break;

    case CRTPClass::EdgePosition::FIRST_EDGE_1D:
      vertexId = localVertexId == 0 ? 0 : 1;
      break;
    case CRTPClass::EdgePosition::LAST_EDGE_1D:
      vertexId = localVertexId == 0 ? this->edgeNumber_ - 1 : this->edgeNumber_;
      break;
    case CRTPClass::EdgePosition::CENTER_1D:
      vertexId = localVertexId == 0 ? edgeId : edgeId + 1;
      break;
  }

  return 0;
}

template<size_t card>
const vector<std::array<SimplexId, 2>> *
  ImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getEdges)() {

  if(this->edgeList_.empty()) {
    Timer t;

    this->edgeList_.resize(this->edgeNumber_);
    for(SimplexId i = 0; i < this->edgeNumber_; ++i) {
      SimplexId id0, id1;
      this->getEdgeVertexInternal(i, 0, id0);
      this->getEdgeVertexInternal(i, 1, id1);
      this->edgeList_[i] = {id0, id1};
    }

    this->printMsg(
      "Built " + to_string(this->edgeNumber_) + " edges.", 1, t.getElapsedTime(), 1);
  }

  return &this->edgeList_;
}

template <size_t card, typename Derived>
SimplexId ImplicitTriangulationCRTP<card, Derived>::getEdgeTriangleNumberInternal(
  const SimplexId &edgeId) const {
  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= this->edgeNumber_)
    return -1;
#endif

  switch(this->underlying().getEdgePosition(edgeId)) {
    case CRTPClass::EdgePosition::L_xnn_3D:
    case CRTPClass::EdgePosition::H_nyn_3D:
    case CRTPClass::EdgePosition::P_nnz_3D:
    case CRTPClass::EdgePosition::D4_3D:
      return 6;
    case CRTPClass::EdgePosition::L_x0n_3D:
    case CRTPClass::EdgePosition::L_xNn_3D:
    case CRTPClass::EdgePosition::L_xn0_3D:
    case CRTPClass::EdgePosition::L_xnN_3D:
    case CRTPClass::EdgePosition::H_ny0_3D:
    case CRTPClass::EdgePosition::H_nyN_3D:
    case CRTPClass::EdgePosition::H_0yn_3D:
    case CRTPClass::EdgePosition::H_Nyn_3D:
    case CRTPClass::EdgePosition::P_n0z_3D:
    case CRTPClass::EdgePosition::P_nNz_3D:
    case CRTPClass::EdgePosition::P_0nz_3D:
    case CRTPClass::EdgePosition::P_Nnz_3D:
    case CRTPClass::EdgePosition::D1_xyn_3D:
    case CRTPClass::EdgePosition::D2_nyz_3D:
    case CRTPClass::EdgePosition::D3_xnz_3D:
      return 4;
    case CRTPClass::EdgePosition::L_x00_3D:
    case CRTPClass::EdgePosition::L_xNN_3D:
    case CRTPClass::EdgePosition::H_0yN_3D:
    case CRTPClass::EdgePosition::H_Ny0_3D:
    case CRTPClass::EdgePosition::P_0Nz_3D:
    case CRTPClass::EdgePosition::P_N0z_3D:
    case CRTPClass::EdgePosition::D1_xy0_3D:
    case CRTPClass::EdgePosition::D1_xyN_3D:
    case CRTPClass::EdgePosition::D2_0yz_3D:
    case CRTPClass::EdgePosition::D2_Nyz_3D:
    case CRTPClass::EdgePosition::D3_x0z_3D:
    case CRTPClass::EdgePosition::D3_xNz_3D:
      return 3;
    case CRTPClass::EdgePosition::L_xN0_3D:
    case CRTPClass::EdgePosition::L_x0N_3D:
    case CRTPClass::EdgePosition::H_0y0_3D:
    case CRTPClass::EdgePosition::H_NyN_3D:
    case CRTPClass::EdgePosition::P_00z_3D:
    case CRTPClass::EdgePosition::P_NNz_3D:
    case CRTPClass::EdgePosition::L_xn_2D:
    case CRTPClass::EdgePosition::H_ny_2D:
    case CRTPClass::EdgePosition::D1_2D:
      return 2;
    case CRTPClass::EdgePosition::L_x0_2D:
    case CRTPClass::EdgePosition::L_xN_2D:
    case CRTPClass::EdgePosition::H_0y_2D:
    case CRTPClass::EdgePosition::H_Ny_2D:
      return 1;

    default: // 1D
      break;
  }

  return 0;
}

template <size_t card, typename Derived>
int ImplicitTriangulationCRTP<card, Derived>::getEdgeTriangleInternal(
  const SimplexId &edgeId,
  const int &localTriangleId,
  SimplexId &triangleId) const {
  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(localTriangleId < 0
     or localTriangleId >= getEdgeTriangleNumberInternal(edgeId))
    return -1;
#endif

  const auto &p = this->underlying().getEdgeCoords(edgeId);

  switch(this->underlying().getEdgePosition(edgeId)) {
    case CRTPClass::EdgePosition::L_xnn_3D:
      triangleId = this->getEdgeTriangleL_xnn(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::L_xn0_3D:
      triangleId = this->getEdgeTriangleL_xn0(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::L_xnN_3D:
      triangleId = this->getEdgeTriangleL_xnN(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::L_x0n_3D:
      triangleId = this->getEdgeTriangleL_x0n(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::L_x00_3D:
      triangleId = this->getEdgeTriangleL_x00(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::L_x0N_3D:
      triangleId = this->getEdgeTriangleL_x0N(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::L_xNn_3D:
      triangleId = this->getEdgeTriangleL_xNn(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::L_xN0_3D:
      triangleId = this->getEdgeTriangleL_xN0(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::L_xNN_3D:
      triangleId = this->getEdgeTriangleL_xNN(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::H_nyn_3D:
      triangleId = this->getEdgeTriangleH_nyn(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::H_ny0_3D:
      triangleId = this->getEdgeTriangleH_ny0(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::H_nyN_3D:
      triangleId = this->getEdgeTriangleH_nyN(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::H_0yn_3D:
      triangleId = this->getEdgeTriangleH_0yn(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::H_0y0_3D:
      triangleId = this->getEdgeTriangleH_0y0(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::H_0yN_3D:
      triangleId = this->getEdgeTriangleH_0yN(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::H_Nyn_3D:
      triangleId = this->getEdgeTriangleH_Nyn(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::H_Ny0_3D:
      triangleId = this->getEdgeTriangleH_Ny0(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::H_NyN_3D:
      triangleId = this->getEdgeTriangleH_NyN(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::P_nnz_3D:
      triangleId = this->getEdgeTriangleP_nnz(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::P_n0z_3D:
      triangleId = this->getEdgeTriangleP_n0z(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::P_nNz_3D:
      triangleId = this->getEdgeTriangleP_nNz(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::P_0nz_3D:
      triangleId = this->getEdgeTriangleP_0nz(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::P_00z_3D:
      triangleId = this->getEdgeTriangleP_00z(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::P_0Nz_3D:
      triangleId = this->getEdgeTriangleP_0Nz(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::P_Nnz_3D:
      triangleId = this->getEdgeTriangleP_Nnz(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::P_N0z_3D:
      triangleId = this->getEdgeTriangleP_N0z(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::P_NNz_3D:
      triangleId = this->getEdgeTriangleP_NNz(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::D1_xyn_3D:
      triangleId = this->getEdgeTriangleD1_xyn(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::D1_xy0_3D:
      triangleId = this->getEdgeTriangleD1_xy0(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::D1_xyN_3D:
      triangleId = this->getEdgeTriangleD1_xyN(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::D2_nyz_3D:
      triangleId = this->getEdgeTriangleD2_nyz(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::D2_0yz_3D:
      triangleId = this->getEdgeTriangleD2_0yz(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::D2_Nyz_3D:
      triangleId = this->getEdgeTriangleD2_Nyz(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::D3_xnz_3D:
      triangleId = this->getEdgeTriangleD3_xnz(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::D3_x0z_3D:
      triangleId = this->getEdgeTriangleD3_x0z(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::D3_xNz_3D:
      triangleId = this->getEdgeTriangleD3_xNz(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::D4_3D:
      triangleId = this->getEdgeTriangleD4_xyz(p.data(), localTriangleId);
      break;

    case CRTPClass::EdgePosition::L_xn_2D:
      triangleId = this->getEdgeTriangleL_xn(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::L_x0_2D:
      triangleId = this->getEdgeTriangleL_x0(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::L_xN_2D:
      triangleId = this->getEdgeTriangleL_xN(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::H_ny_2D:
      triangleId = this->getEdgeTriangleH_ny(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::H_0y_2D:
      triangleId = this->getEdgeTriangleH_0y(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::H_Ny_2D:
      triangleId = this->getEdgeTriangleH_Ny(p.data(), localTriangleId);
      break;
    case CRTPClass::EdgePosition::D1_2D:
      triangleId = this->getEdgeTriangleD1_xy(p.data(), localTriangleId);
      break;

    default: // 1D
      triangleId = -1;
      break;
  }

  return 0;
}

template<size_t card>
const vector<vector<SimplexId>> *
  ImplicitTriangulation<card>::getEdgeTrianglesInternal() {
  if(this->edgeTriangleList_.empty()) {
    Timer t;

    this->edgeTriangleList_.resize(this->edgeNumber_);
    for(SimplexId i = 0; i < this->edgeNumber_; ++i) {
      this->edgeTriangleList_[i].resize(this->getEdgeTriangleNumberInternal(i));
      for(SimplexId j = 0; j < (SimplexId)this->edgeTriangleList_[i].size(); ++j)
        this->getEdgeTriangleInternal(i, j, this->edgeTriangleList_[i][j]);
    }

    this->printMsg("Built " + to_string(this->edgeNumber_) + " edge triangles.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->edgeTriangleList_;
}

template<size_t card>
SimplexId ImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getEdgeLinkNumber)(
  const SimplexId &edgeId) const {
  return TTK_TRIANGULATION_INTERNAL(this->getEdgeStarNumber)(edgeId);
}

template <size_t card, typename Derived>
int ImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(getEdgeLink)(
  const SimplexId &edgeId, const int &localLinkId, SimplexId &linkId) const {

  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(localLinkId < 0 or localLinkId >= getEdgeLinkNumber(edgeId))
    return -1;
#endif

  const auto &p = this->underlying().getEdgeCoords(edgeId);

  switch(this->underlying().getEdgePosition(edgeId)) {
  CASE_EDGE_POSITION_L_3D:
    linkId = this->getEdgeLinkL(p.data(), localLinkId);
    break;
  CASE_EDGE_POSITION_H_3D:
    linkId = this->getEdgeLinkH(p.data(), localLinkId);
    break;
  CASE_EDGE_POSITION_P_3D:
    linkId = this->getEdgeLinkP(p.data(), localLinkId);
    break;
  CASE_EDGE_POSITION_D1_3D:
    linkId = this->getEdgeLinkD1(p.data(), localLinkId);
    break;
  CASE_EDGE_POSITION_D2_3D:
    linkId = this->getEdgeLinkD2(p.data(), localLinkId);
    break;
  CASE_EDGE_POSITION_D3_3D:
    linkId = this->getEdgeLinkD3(p.data(), localLinkId);
    break;
    case CRTPClass::EdgePosition::D4_3D:
      linkId = this->getEdgeLinkD4(p.data(), localLinkId);
      break;

    CASE_EDGE_POSITION_L_2D:
      linkId = this->getEdgeLink2dL(p.data(), localLinkId);
      break;
    CASE_EDGE_POSITION_H_2D:
      linkId = this->getEdgeLink2dH(p.data(), localLinkId);
      break;
    case CRTPClass::EdgePosition::D1_2D:
      linkId = this->getEdgeLink2dD1(p.data(), localLinkId);
      break;

    default: // 1D
      linkId = -1;
      break;
  }

  return 0;
}

template<size_t card>
const vector<vector<SimplexId>> *
  ImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getEdgeLinks)() {

  if(this->edgeLinkList_.empty()) {
    Timer t;

    this->edgeLinkList_.resize(this->edgeNumber_);
    for(SimplexId i = 0; i < this->edgeNumber_; ++i) {
      this->edgeLinkList_[i].resize(getEdgeLinkNumber(i));
      for(SimplexId j = 0; j < (SimplexId)this->edgeLinkList_[i].size(); ++j)
        this->getEdgeLink(i, j, this->edgeLinkList_[i][j]);
    }

    this->printMsg("Built " + to_string(this->edgeNumber_) + " edge links.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->edgeLinkList_;
}

template <size_t card, typename Derived>
SimplexId ImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  getEdgeStarNumber)(const SimplexId &edgeId) const {

  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(edgeId < 0 or edgeId >= this->edgeNumber_)
    return -1;
#endif

  switch(this->underlying().getEdgePosition(edgeId)) {
    case CRTPClass::EdgePosition::L_xnn_3D:
    case CRTPClass::EdgePosition::H_nyn_3D:
    case CRTPClass::EdgePosition::P_nnz_3D:
    case CRTPClass::EdgePosition::D4_3D:
      return 6;
    case CRTPClass::EdgePosition::D1_xyn_3D:
    case CRTPClass::EdgePosition::D2_nyz_3D:
    case CRTPClass::EdgePosition::D3_xnz_3D:
      return 4;
    case CRTPClass::EdgePosition::L_x0n_3D:
    case CRTPClass::EdgePosition::L_xNn_3D:
    case CRTPClass::EdgePosition::L_xn0_3D:
    case CRTPClass::EdgePosition::L_xnN_3D:
    case CRTPClass::EdgePosition::H_ny0_3D:
    case CRTPClass::EdgePosition::H_nyN_3D:
    case CRTPClass::EdgePosition::H_0yn_3D:
    case CRTPClass::EdgePosition::H_Nyn_3D:
    case CRTPClass::EdgePosition::P_n0z_3D:
    case CRTPClass::EdgePosition::P_nNz_3D:
    case CRTPClass::EdgePosition::P_0nz_3D:
    case CRTPClass::EdgePosition::P_Nnz_3D:
      return 3;
    case CRTPClass::EdgePosition::L_x00_3D:
    case CRTPClass::EdgePosition::L_xNN_3D:
    case CRTPClass::EdgePosition::H_0yN_3D:
    case CRTPClass::EdgePosition::H_Ny0_3D:
    case CRTPClass::EdgePosition::P_0Nz_3D:
    case CRTPClass::EdgePosition::P_N0z_3D:
    case CRTPClass::EdgePosition::D1_xy0_3D:
    case CRTPClass::EdgePosition::D1_xyN_3D:
    case CRTPClass::EdgePosition::D2_0yz_3D:
    case CRTPClass::EdgePosition::D2_Nyz_3D:
    case CRTPClass::EdgePosition::D3_x0z_3D:
    case CRTPClass::EdgePosition::D3_xNz_3D:
    case CRTPClass::EdgePosition::L_xn_2D:
    case CRTPClass::EdgePosition::H_ny_2D:
    case CRTPClass::EdgePosition::D1_2D:
      return 2;
    case CRTPClass::EdgePosition::L_xN0_3D:
    case CRTPClass::EdgePosition::L_x0N_3D:
    case CRTPClass::EdgePosition::H_0y0_3D:
    case CRTPClass::EdgePosition::H_NyN_3D:
    case CRTPClass::EdgePosition::P_00z_3D:
    case CRTPClass::EdgePosition::P_NNz_3D:
    case CRTPClass::EdgePosition::L_x0_2D:
    case CRTPClass::EdgePosition::L_xN_2D:
    case CRTPClass::EdgePosition::H_0y_2D:
    case CRTPClass::EdgePosition::H_Ny_2D:
      return 1;

    default: // 1D
      break;
  }

  return 0;
}

template <size_t card, typename Derived>
int ImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(getEdgeStar)(
  const SimplexId &edgeId, const int &localStarId, SimplexId &starId) const {

  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(localStarId < 0 or localStarId >= this->getEdgeStarNumber(edgeId))
    return -1;
#endif

  const auto &p = this->underlying().getEdgeCoords(edgeId);

  switch(this->underlying().getEdgePosition(edgeId)) {
  CASE_EDGE_POSITION_L_3D:
    starId = this->getEdgeStarL(p.data(), localStarId);
    break;
  CASE_EDGE_POSITION_H_3D:
    starId = this->getEdgeStarH(p.data(), localStarId);
    break;
  CASE_EDGE_POSITION_P_3D:
    starId = this->getEdgeStarP(p.data(), localStarId);
    break;
  CASE_EDGE_POSITION_D1_3D:
    starId = this->getEdgeStarD1(p.data(), localStarId);
    break;
  CASE_EDGE_POSITION_D2_3D:
    starId = this->getEdgeStarD2(p.data(), localStarId);
    break;
  CASE_EDGE_POSITION_D3_3D:
    starId = this->getEdgeStarD3(p.data(), localStarId);
    break;
    case CRTPClass::EdgePosition::D4_3D:
      starId
        = p[2] * this->tetshift_[1] + p[1] * this->tetshift_[0] + p[0] * 6 + localStarId;
      break;

    CASE_EDGE_POSITION_L_2D:
      starId = this->getEdgeStar2dL(p.data(), localStarId);
      break;
    CASE_EDGE_POSITION_H_2D:
      starId = this->getEdgeStar2dH(p.data(), localStarId);
      break;
    case CRTPClass::EdgePosition::D1_2D:
      starId = p[0] * 2 + p[1] * this->tshift_[0] + localStarId;
      break;

    default: // 1D
      starId = -1;
      break;
  }

  return 0;
}

template<size_t card>
const vector<vector<SimplexId>> *
  ImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getEdgeStars)() {

  if(this->edgeStarList_.empty()) {
    Timer t;

    this->edgeStarList_.resize(this->edgeNumber_);
    for(SimplexId i = 0; i < this->edgeNumber_; ++i) {
      this->edgeStarList_[i].resize(this->getEdgeStarNumber(i));
      for(SimplexId j = 0; j < (SimplexId)this->edgeStarList_[i].size(); ++j)
        this->getEdgeStar(i, j, this->edgeStarList_[i][j]);
    }

    this->printMsg("Built " + to_string(this->edgeNumber_) + " edge stars.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->edgeStarList_;
}

template <size_t card, typename Derived>
int ImplicitTriangulationCRTP<card, Derived>::getTriangleVertexInternal(
  const SimplexId &triangleId,
  const int &localVertexId,
  SimplexId &vertexId) const {
  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
  if(localVertexId < 0 or localVertexId >= 3)
    return -2;
#endif

  //    e--------f
  //   /|       /|
  //  / |      / |
  // a--g-----b--h
  // | /      | /
  // |/       |/
  // c--------d
  //
  // Classement des "Triangles" et dans cet ordre:
  // F: face (type abc/bcd)
  // C: cote (type abe/bef)
  // H: haut (type acg/aeg)
  // D1: diagonale1 (type bdg/beg)
  // D2: diagonale2 (type abg/bgh)
  // D3: diagonale3 (type bcg/bfg)

  const auto &p = this->underlying().getTriangleCoords(triangleId);
  vertexId = -1;

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
      switch(localVertexId) {
        break;
        case 0:
          vertexId = p[0] / 2 + p[1] * this->vshift_[0];
          break;
        case 1:
          vertexId = p[0] / 2 + p[1] * this->vshift_[0] + 1;
          break;
        case 2:
          vertexId = p[0] / 2 + p[1] * this->vshift_[0] + this->vshift_[0];
          break;
      }
      break;
    case CRTPClass::TrianglePosition::BOTTOM_2D:
      switch(localVertexId) {
        break;
        case 0:
          vertexId = p[0] / 2 + p[1] * this->vshift_[0] + 1;
          break;
        case 1:
          vertexId = p[0] / 2 + p[1] * this->vshift_[0] + this->vshift_[0] + 1;
          break;
        case 2:
          vertexId = p[0] / 2 + p[1] * this->vshift_[0] + this->vshift_[0];
          break;
      }
  }

  return 0;
}

template <size_t card, typename Derived>
int ImplicitTriangulationCRTP<card, Derived>::getTriangleEdgeInternal(
  const SimplexId &triangleId,
  const int &localEdgeId,
  SimplexId &edgeId) const {
  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
  if(localEdgeId < 0 or localEdgeId >= 3)
    return -2;
#endif

  const auto &p = this->underlying().getTriangleCoords(triangleId);
  const auto par = triangleId % 2;
  edgeId = -1;

  switch(this->underlying().getTrianglePosition(triangleId)) {
    case CRTPClass::TrianglePosition::F_3D:
      edgeId = (par == 1) ? this->getTriangleEdgeF_1(p.data(), localEdgeId)
                          : this->getTriangleEdgeF_0(p.data(), localEdgeId);
      break;
    case CRTPClass::TrianglePosition::H_3D:
      edgeId = (par == 1) ? this->getTriangleEdgeH_1(p.data(), localEdgeId)
                          : this->getTriangleEdgeH_0(p.data(), localEdgeId);
      break;
    case CRTPClass::TrianglePosition::C_3D:
      edgeId = (par == 1) ? this->getTriangleEdgeC_1(p.data(), localEdgeId)
                          : this->getTriangleEdgeC_0(p.data(), localEdgeId);
      break;
    case CRTPClass::TrianglePosition::D1_3D:
      edgeId = (par == 1) ? this->getTriangleEdgeD1_1(p.data(), localEdgeId)
                          : this->getTriangleEdgeD1_0(p.data(), localEdgeId);
      break;
    case CRTPClass::TrianglePosition::D2_3D:
      edgeId = (par == 1) ? this->getTriangleEdgeD2_1(p.data(), localEdgeId)
                          : this->getTriangleEdgeD2_0(p.data(), localEdgeId);
      break;
    case CRTPClass::TrianglePosition::D3_3D:
      edgeId = (par == 1) ? this->getTriangleEdgeD3_1(p.data(), localEdgeId)
                          : this->getTriangleEdgeD3_0(p.data(), localEdgeId);
      break;
    case CRTPClass::TrianglePosition::TOP_2D:
      switch(localEdgeId) {
        break;
        case 0:
          edgeId = p[0] / 2 + p[1] * this->eshift_[0];
          break;
        case 1:
          edgeId = this->esetshift_[0] + p[0] / 2 + p[1] * this->eshift_[2];
          break;
        case 2:
          edgeId = this->esetshift_[1] + p[0] / 2 + p[1] * this->eshift_[4];
          break;
      }
      break;
    case CRTPClass::TrianglePosition::BOTTOM_2D:
      switch(localEdgeId) {
        break;
        case 0:
          edgeId = p[0] / 2 + (p[1] + 1) * this->eshift_[0];
          break;
        case 1:
          edgeId = this->esetshift_[0] + (p[0] + 1) / 2 + p[1] * this->eshift_[2];
          break;
        case 2:
          edgeId = this->esetshift_[1] + p[0] / 2 + p[1] * this->eshift_[4];
          break;
      }
  }

  return 0;
}

template<size_t card>
int ImplicitTriangulation<card>::getTriangleEdgesInternal(
  vector<vector<SimplexId>> &edges) const {
  edges.resize(triangleNumber_);
  for(SimplexId i = 0; i < triangleNumber_; ++i) {
    edges[i].resize(3);
    for(int j = 0; j < 3; ++j)
      this->getTriangleEdgeInternal(i, j, edges[i][j]);
  }
  return 0;
}

template<size_t card>
const vector<vector<SimplexId>> *
  ImplicitTriangulation<card>::getTriangleEdgesInternal() {
  if(this->triangleEdgeVector_.empty()) {
    Timer t;

    getTriangleEdgesInternal(this->triangleEdgeVector_);

    this->printMsg("Built " + to_string(triangleNumber_) + " triangle edges.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->triangleEdgeVector_;
}

template<size_t card>
const vector<std::array<SimplexId, 3>> *
  ImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getTriangles)() {

  if(this->triangleList_.empty()) {
    Timer t;

    this->triangleList_.resize(triangleNumber_);
    for(SimplexId i = 0; i < triangleNumber_; ++i) {
      for(int j = 0; j < 3; ++j)
        this->getTriangleVertexInternal(i, j, this->triangleList_[i][j]);
    }

    this->printMsg("Built " + to_string(triangleNumber_) + " triangles.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->triangleList_;
}

template <size_t card, typename Derived>
int ImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  getTriangleLink)(const SimplexId &triangleId,
                   const int &localLinkId,
                   SimplexId &linkId) const {

  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(localLinkId < 0 or localLinkId >= getTriangleLinkNumber(triangleId))
    return -1;
#endif

  const auto &p = this->underlying().getTriangleCoords(triangleId);

  switch(this->underlying().getTrianglePosition(triangleId)) {
    case CRTPClass::TrianglePosition::F_3D:
      linkId = this->getTriangleLinkF(p.data(), localLinkId);
      break;
    case CRTPClass::TrianglePosition::H_3D:
      linkId = this->getTriangleLinkH(p.data(), localLinkId);
      break;
    case CRTPClass::TrianglePosition::C_3D:
      linkId = this->getTriangleLinkC(p.data(), localLinkId);
      break;
    case CRTPClass::TrianglePosition::D1_3D:
      linkId = this->getTriangleLinkD1(p.data(), localLinkId);
      break;
    case CRTPClass::TrianglePosition::D2_3D:
      linkId = this->getTriangleLinkD2(p.data(), localLinkId);
      break;
    case CRTPClass::TrianglePosition::D3_3D:
      linkId = this->getTriangleLinkD3(p.data(), localLinkId);
      break;
    default: // 2D
      linkId = -1;
      break;
  }

  return 0;
}

template <size_t card>
SimplexId ImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(
  getTriangleLinkNumber)(const SimplexId &triangleId) const {
  return TTK_TRIANGULATION_INTERNAL(this->getTriangleStarNumber)(triangleId);
}

template<size_t card>
const vector<vector<SimplexId>> *
  ImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getTriangleLinks)() {
  if(this->triangleLinkList_.empty()) {
    Timer t;

    this->triangleLinkList_.resize(triangleNumber_);
    for(SimplexId i = 0; i < triangleNumber_; ++i) {
      this->triangleLinkList_[i].resize(getTriangleLinkNumber(i));
      for(SimplexId j = 0; j < (SimplexId)this->triangleLinkList_[i].size(); ++j)
        this->getTriangleLink(i, j, this->triangleLinkList_[i][j]);
    }

    this->printMsg("Built " + to_string(triangleNumber_) + " triangle links.", 1,
             t.getElapsedTime(), 1);
  }
  return &this->triangleLinkList_;
}

template <size_t card, typename Derived>
SimplexId ImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  getTriangleStarNumber)(const SimplexId &triangleId) const {

  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
#endif

  const auto &p = this->underlying().getTriangleCoords(triangleId);

  switch(this->underlying().getTrianglePosition(triangleId)) {
    case CRTPClass::TrianglePosition::F_3D:
      return (p[2] > 0 and p[2] < this->nbvoxels_[2]) ? 2 : 1;
    case CRTPClass::TrianglePosition::H_3D:
      return (p[1] > 0 and p[1] < this->nbvoxels_[1]) ? 2 : 1;
    case CRTPClass::TrianglePosition::C_3D:
      return (p[0] < 2 or p[0] >= (this->dimensions_[0] * 2 - 2)) ? 1 : 2;

    case CRTPClass::TrianglePosition::D1_3D:
    case CRTPClass::TrianglePosition::D2_3D:
    case CRTPClass::TrianglePosition::D3_3D:
      return 2;
    default: // 2D
      break;
  }
  return 0;
}

template <size_t card, typename Derived>
int ImplicitTriangulationCRTP<card, Derived>::TTK_TRIANGULATION_INTERNAL(
  getTriangleStar)(const SimplexId &triangleId,
                   const int &localStarId,
                   SimplexId &starId) const {

  using CRTPClass = ttk::ImplicitTriangulationCRTP<card, Derived>;
#ifndef TTK_ENABLE_KAMIKAZE
  if(localStarId < 0 or localStarId >= this->getTriangleStarNumber(triangleId))
    return -1;
#endif

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
    default: // 2D
      starId = -1;
      break;
  }

  return 0;
}

template<size_t card>
const vector<vector<SimplexId>> *
  ImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getTriangleStars)() {

  if(this->triangleStarList_.empty()) {
    Timer t;

    this->triangleStarList_.resize(triangleNumber_);
    for(SimplexId i = 0; i < triangleNumber_; ++i) {
      this->triangleStarList_[i].resize(this->getTriangleStarNumber(i));
      for(SimplexId j = 0; j < (SimplexId)this->triangleStarList_[i].size(); ++j)
        this->getTriangleStar(i, j, this->triangleStarList_[i][j]);
    }

    this->printMsg("Built " + to_string(triangleNumber_) + " triangle stars.", 1,
             t.getElapsedTime(), 1);
  }
  return &this->triangleStarList_;
}

template <size_t card, typename Derived>
SimplexId ImplicitTriangulationCRTP<card, Derived>::getTriangleNeighborNumber(
  const SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(triangleId < 0 or triangleId >= triangleNumber_)
    return -1;
#endif

  if(this->dimensionality_ == 2) {
    const auto &p = this->underlying().getTriangleCoords(triangleId);
    const SimplexId id = triangleId % 2;

    if(id) {
      if(p[0] / 2 == this->nbvoxels_[this->Di_] - 1 and p[1] == this->nbvoxels_[this->Dj_] - 1)
        return 1;
      else if(p[0] / 2 == this->nbvoxels_[this->Di_] - 1 or p[1] == this->nbvoxels_[this->Dj_] - 1)
        return 2;
      else
        return 3;
    } else {
      if(p[0] == 0 and p[1] == 0)
        return 1;
      else if(p[0] == 0 or p[1] == 0)
        return 2;
      else
        return 3;
    }
  }

  return 0;
}

template <size_t card, typename Derived>
int ImplicitTriangulationCRTP<card, Derived>::getTriangleNeighbor(
  const SimplexId &triangleId,
  const int &localNeighborId,
  SimplexId &neighborId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localNeighborId < 0
     or localNeighborId >= getTriangleNeighborNumber(triangleId))
    return -1;
#endif

  neighborId = -1;

  if(this->dimensionality_ == 2) {
    const auto &p = this->underlying().getTriangleCoords(triangleId);
    const SimplexId id = triangleId % 2;

    if(id) {
      if(p[0] / 2 == this->nbvoxels_[this->Di_] - 1 and p[1] == this->nbvoxels_[this->Dj_] - 1)
        neighborId = triangleId - 1;
      else if(p[0] / 2 == this->nbvoxels_[this->Di_] - 1) {
        switch(localNeighborId) {
          case 0:
            neighborId = triangleId - 1;
            break;
          case 1:
            neighborId = triangleId + this->tshift_[0] - 1;
            break;
        }
      } else if(p[1] == this->nbvoxels_[this->Dj_] - 1) {
        switch(localNeighborId) {
          case 0:
            neighborId = triangleId - 1;
            break;
          case 1:
            neighborId = triangleId + 1;
            break;
        }
      } else {
        switch(localNeighborId) {
          case 0:
            neighborId = triangleId - 1;
            break;
          case 1:
            neighborId = triangleId + 1;
            break;
          case 2:
            neighborId = triangleId + this->tshift_[0] - 1;
            break;
        }
      }
    } else {
      if(p[0] == 0 and p[1] == 0)
        neighborId = triangleId + 1;
      else if(p[0] == 0) {
        switch(localNeighborId) {
          case 0:
            neighborId = triangleId + 1;
            break;
          case 1:
            neighborId = triangleId - this->tshift_[0] + 1;
            break;
        }
      } else if(p[1] == 0) {
        switch(localNeighborId) {
          case 0:
            neighborId = triangleId + 1;
            break;
          case 1:
            neighborId = triangleId - 1;
            break;
        }
      } else {
        switch(localNeighborId) {
          case 0:
            neighborId = triangleId + 1;
            break;
          case 1:
            neighborId = triangleId - 1;
            break;
          case 2:
            neighborId = triangleId - this->tshift_[0] + 1;
            break;
        }
      }
    }
  }

  return 0;
}

template<size_t card>
int ImplicitTriangulation<card>::getTriangleNeighbors(
  vector<vector<SimplexId>> &neighbors) {
  neighbors.resize(triangleNumber_);
  for(SimplexId i = 0; i < triangleNumber_; ++i) {
    neighbors[i].resize(getTriangleNeighborNumber(i));
    for(SimplexId j = 0; j < (SimplexId)neighbors[i].size(); ++j)
      getTriangleNeighbor(i, j, neighbors[i][j]);
  }
  return 0;
}

template <size_t card, typename Derived>
int ImplicitTriangulationCRTP<card, Derived>::getTetrahedronVertex(
  const SimplexId &tetId, const int &localVertexId, SimplexId &vertexId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(tetId < 0 or tetId >= tetrahedronNumber_)
    return -1;
  if(localVertexId < 0 or localVertexId >= 4)
    return -2;
#endif

  vertexId = -1;

  if(this->dimensionality_ == 3) {
    const SimplexId id = tetId % 6;
    const auto &c = this->underlying().getTetrahedronCoords(tetId);
    const auto p{c.data()};

    switch(id) {
      case 0:
        vertexId = this->getTetrahedronVertexABCG(p, localVertexId);
        break;
      case 1:
        vertexId = this->getTetrahedronVertexBCDG(p, localVertexId);
        break;
      case 2:
        vertexId = this->getTetrahedronVertexABEG(p, localVertexId);
        break;
      case 3:
        vertexId = this->getTetrahedronVertexBEFG(p, localVertexId);
        break;
      case 4:
        vertexId = this->getTetrahedronVertexBFGH(p, localVertexId);
        break;
      case 5:
        vertexId = this->getTetrahedronVertexBDGH(p, localVertexId);
        break;
    }
  }
  return 0;
}

template <size_t card, typename Derived>
int ImplicitTriangulationCRTP<card, Derived>::getTetrahedronEdge(
  const SimplexId &tetId, const int &localEdgeId, SimplexId &edgeId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(tetId < 0 or tetId >= tetrahedronNumber_)
    return -1;
  if(localEdgeId < 0 or localEdgeId >= 6)
    return -2;
#endif

  edgeId = -1;

  if(this->dimensionality_ == 3) {
    const SimplexId id = tetId % 6;
    const auto &c = this->underlying().getTetrahedronCoords(tetId);
    const auto p{c.data()};

    switch(id) {
      case 0:
        edgeId = this->getTetrahedronEdgeABCG(p, localEdgeId);
        break;
      case 1:
        edgeId = this->getTetrahedronEdgeBCDG(p, localEdgeId);
        break;
      case 2:
        edgeId = this->getTetrahedronEdgeABEG(p, localEdgeId);
        break;
      case 3:
        edgeId = this->getTetrahedronEdgeBEFG(p, localEdgeId);
        break;
      case 4:
        edgeId = this->getTetrahedronEdgeBFGH(p, localEdgeId);
        break;
      case 5:
        edgeId = this->getTetrahedronEdgeBDGH(p, localEdgeId);
        break;
    }
  }

  return 0;
}

template<size_t card>
int ImplicitTriangulation<card>::getTetrahedronEdges(
  vector<vector<SimplexId>> &edges) const {
  edges.resize(tetrahedronNumber_);
  for(SimplexId i = 0; i < tetrahedronNumber_; ++i) {
    edges[i].resize(6);
    for(int j = 0; j < 6; ++j)
      getTetrahedronEdge(i, j, edges[i][j]);
  }

  return 0;
}

template <size_t card, typename Derived>
int ImplicitTriangulationCRTP<card, Derived>::getTetrahedronTriangle(
  const SimplexId &tetId,
  const int &localTriangleId,
  SimplexId &triangleId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(tetId < 0 or tetId >= tetrahedronNumber_)
    return -1;
  if(localTriangleId < 0 or localTriangleId >= 4)
    return -2;
#endif

  triangleId = -1;

  if(this->dimensionality_ == 3) {
    const SimplexId id = tetId % 6;
    const auto &c = this->underlying().getTetrahedronCoords(tetId);
    const auto p{c.data()};

    switch(id) {
      case 0:
        triangleId = this->getTetrahedronTriangleABCG(p, localTriangleId);
        break;
      case 1:
        triangleId = this->getTetrahedronTriangleBCDG(p, localTriangleId);
        break;
      case 2:
        triangleId = this->getTetrahedronTriangleABEG(p, localTriangleId);
        break;
      case 3:
        triangleId = this->getTetrahedronTriangleBEFG(p, localTriangleId);
        break;
      case 4:
        triangleId = this->getTetrahedronTriangleBFGH(p, localTriangleId);
        break;
      case 5:
        triangleId = this->getTetrahedronTriangleBDGH(p, localTriangleId);
        break;
    }
  }

  return 0;
}

template<size_t card>
int ImplicitTriangulation<card>::getTetrahedronTriangles(
  vector<vector<SimplexId>> &triangles) const {
  triangles.resize(tetrahedronNumber_);
  for(SimplexId i = 0; i < tetrahedronNumber_; ++i) {
    triangles[i].resize(4);
    for(int j = 0; j < 4; ++j)
      getTetrahedronTriangle(i, j, triangles[i][j]);
  }

  return 0;
}

template <size_t card, typename Derived>
SimplexId ImplicitTriangulationCRTP<card, Derived>::getTetrahedronNeighborNumber(
  const SimplexId &tetId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(tetId < 0 or tetId >= tetrahedronNumber_)
    return -1;
#endif

  if(this->dimensionality_ == 3) {
    const SimplexId id = tetId % 6;
    const auto &c = this->underlying().getTetrahedronCoords(tetId);
    const auto p{c.data()};

    switch(id) {
      case 0: // ABCG
        if(p[0] == 0 and p[2] == 0)
          return 2;
        else if(p[0] == 0 or p[2] == 0)
          return 3;
        else
          return 4;
        break;
      case 1: // BCDG
        if(p[1] == this->nbvoxels_[1] - 1 and p[2] == 0)
          return 2;
        else if(p[1] == this->nbvoxels_[1] - 1 or p[2] == 0)
          return 3;
        else
          return 4;
        break;
      case 2: // ABEG
        if(p[0] == 0 and p[1] == 0)
          return 2;
        else if(p[0] == 0 or p[1] == 0)
          return 3;
        else
          return 4;
        break;
      case 3: // BEFG
        if(p[1] == 0 and p[2] == this->nbvoxels_[2] - 1)
          return 2;
        else if(p[1] == 0 or p[2] == this->nbvoxels_[2] - 1)
          return 3;
        else
          return 4;
        break;
      case 4: // BFGH
        if(p[0] == this->nbvoxels_[0] - 1 and p[2] == this->nbvoxels_[2] - 1)
          return 2;
        else if(p[0] == this->nbvoxels_[0] - 1 or p[2] == this->nbvoxels_[2] - 1)
          return 3;
        else
          return 4;
        break;
      case 5: // BDGH
        if(p[0] == this->nbvoxels_[0] - 1 and p[1] == this->nbvoxels_[1] - 1)
          return 2;
        else if(p[0] == this->nbvoxels_[0] - 1 or p[1] == this->nbvoxels_[1] - 1)
          return 3;
        else
          return 4;
        break;
    }
  }

  return 0;
}

template <size_t card, typename Derived>
int ImplicitTriangulationCRTP<card, Derived>::getTetrahedronNeighbor(
  const SimplexId &tetId,
  const int &localNeighborId,
  SimplexId &neighborId) const {
#ifndef TTK_ENABLE_KAMIKAZE
  if(localNeighborId < 0
     or localNeighborId >= this->getTetrahedronNeighborNumber(tetId))
    return -1;
#endif

  neighborId = -1;

  if(this->dimensionality_ == 3) {
    const SimplexId id = tetId % 6;
    const auto &c = this->underlying().getTetrahedronCoords(tetId);
    const auto p{c.data()};

    switch(id) {
      case 0:
        neighborId = this->getTetrahedronNeighborABCG(tetId, p, localNeighborId);
        break;
      case 1:
        neighborId = this->getTetrahedronNeighborBCDG(tetId, p, localNeighborId);
        break;
      case 2:
        neighborId = this->getTetrahedronNeighborABEG(tetId, p, localNeighborId);
        break;
      case 3:
        neighborId = this->getTetrahedronNeighborBEFG(tetId, p, localNeighborId);
        break;
      case 4:
        neighborId = this->getTetrahedronNeighborBFGH(tetId, p, localNeighborId);
        break;
      case 5:
        neighborId = this->getTetrahedronNeighborBDGH(tetId, p, localNeighborId);
        break;
    }
  }

  return 0;
}

template<size_t card>
int ImplicitTriangulation<card>::getTetrahedronNeighbors(
  vector<vector<SimplexId>> &neighbors) {
  neighbors.resize(tetrahedronNumber_);
  for(SimplexId i = 0; i < tetrahedronNumber_; ++i) {
    neighbors[i].resize(getTetrahedronNeighborNumber(i));
    for(SimplexId j = 0; j < (SimplexId)neighbors[i].size(); ++j)
      getTetrahedronNeighbor(i, j, neighbors[i][j]);
  }

  return 0;
}

template <size_t card>
SimplexId ImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(
  getCellVertexNumber)(const SimplexId & /*cellId*/) const {
  return this->dimensionality_ + 1;
}

template<size_t card>
int ImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getCellVertex)(
  const SimplexId &cellId,
  const int &localVertexId,
  SimplexId &vertexId) const {

  if(this->dimensionality_ == 3)
    getTetrahedronVertex(cellId, localVertexId, vertexId);
  else if(this->dimensionality_ == 2)
    this->getTriangleVertexInternal(cellId, localVertexId, vertexId);
  else if(this->dimensionality_ == 1)
    this->getEdgeVertexInternal(cellId, localVertexId, vertexId);

  return 0;
}

template <size_t card>
SimplexId ImplicitTriangulation<card>::getCellEdgeNumberInternal(
  const SimplexId & /*cellId*/) const {
  if(this->dimensionality_ == 3)
    return 6;
  else if(this->dimensionality_ == 2)
    return 3;

  return 0;
}

template<size_t card>
int ImplicitTriangulation<card>::getCellEdgeInternal(const SimplexId &cellId,
                                               const int &localEdgeId,
                                               SimplexId &edgeId) const {
  if(this->dimensionality_ == 3)
    getTetrahedronEdge(cellId, localEdgeId, edgeId);
  else if(this->dimensionality_ == 2)
    this->getTriangleEdgeInternal(cellId, localEdgeId, edgeId);
  else if(this->dimensionality_ == 1)
    getCellNeighbor(cellId, localEdgeId, edgeId);

  return 0;
}

template <size_t card>
const vector<vector<SimplexId>> *ImplicitTriangulation<card>::getCellEdgesInternal() {
  if(this->cellEdgeVector_.empty()) {
    Timer t;

    if(this->dimensionality_ == 3)
      getTetrahedronEdges(this->cellEdgeVector_);
    else if(this->dimensionality_ == 2)
      getTriangleEdgesInternal(this->cellEdgeVector_);

    this->printMsg("Built " + to_string(cellNumber_) + " cell edges.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->cellEdgeVector_;
}

template<size_t card>
int ImplicitTriangulation<card>::getCellTriangleInternal(
  const SimplexId &cellId,
  const int &localTriangleId,
  SimplexId &triangleId) const {
  if(this->dimensionality_ == 3)
    getTetrahedronTriangle(cellId, localTriangleId, triangleId);

  return 0;
}

template<size_t card>
const vector<vector<SimplexId>> *
  ImplicitTriangulation<card>::getCellTrianglesInternal() {
  if(this->cellTriangleVector_.empty()) {
    Timer t;

    if(this->dimensionality_ == 3)
      getTetrahedronTriangles(this->cellTriangleVector_);

    this->printMsg("Built " + to_string(cellNumber_) + " cell triangles.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->cellTriangleVector_;
}

template <size_t card>
SimplexId ImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(
  getCellNeighborNumber)(const SimplexId &cellId) const {
  if(this->dimensionality_ == 3)
    return getTetrahedronNeighborNumber(cellId);
  else if(this->dimensionality_ == 2)
    return getTriangleNeighborNumber(cellId);
  else if(this->dimensionality_ == 1) {
    this->printErr("getCellNeighborNumber() not implemented in 1D! (TODO)");
    return -1;
  }

  return 0;
}

template<size_t card>
int ImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getCellNeighbor)(
  const SimplexId &cellId,
  const int &localNeighborId,
  SimplexId &neighborId) const {
  if(this->dimensionality_ == 3)
    getTetrahedronNeighbor(cellId, localNeighborId, neighborId);
  else if(this->dimensionality_ == 2)
    getTriangleNeighbor(cellId, localNeighborId, neighborId);
  else if(this->dimensionality_ == 1) {
    this->printErr("getCellNeighbor() not implemented in 1D! (TODO)");
    return -1;
  }

  return 0;
}

template<size_t card>
const vector<vector<SimplexId>> *
  ImplicitTriangulation<card>::TTK_TRIANGULATION_INTERNAL(getCellNeighbors)() {
  if(this->cellNeighborList_.empty()) {
    Timer t;

    if(this->dimensionality_ == 3)
      getTetrahedronNeighbors(this->cellNeighborList_);
    else if(this->dimensionality_ == 2)
      getTriangleNeighbors(this->cellNeighborList_);
    else if(this->dimensionality_ == 1) {
      this->printErr("getCellNeighbors() not implemented in 1D! (TODO)");
      return nullptr;
    }

    this->printMsg("Built " + to_string(cellNumber_) + " cell neighbors.", 1,
             t.getElapsedTime(), 1);
  }

  return &this->cellNeighborList_;
}

template<size_t card>
int ImplicitTriangulation<card>::preconditionVertexNeighborsInternal() {
  // V(abcdefgh)=V(g)+V(d)::{g,h}+V(h)::{g}+V(b)::{c,d,g,h}
  this->vertexNeighborABCDEFGH_ = {
    -this->vshift_[0] - this->vshift_[1], // a
    1 - this->vshift_[0] - this->vshift_[1], // b
    -this->vshift_[1], // c
    1 - this->vshift_[1], // d
    -this->vshift_[0], // e
    1 - this->vshift_[0], // f
    1, // h
    -1 + this->vshift_[1], // V(d)::{g}
    this->vshift_[1], // V(d)::{h}
    -1, // V(h)::{g}
    -1 + this->vshift_[0], // V(b)::{c}
    this->vshift_[0], // V(b)::{d}
    -1 + this->vshift_[0] + this->vshift_[1], // V(b)::{g}
    this->vshift_[0] + this->vshift_[1] // V(b)::{h}
  };

  // V(abdc)=V(b)+V(d)::{b}+V(c)::{b}+V(a)::{b}
  this->vertexNeighborABCD_ = {
    -1, // a
    -1 + this->vshift_[0], // c
    this->vshift_[0], // d
    -1 + this->vshift_[1], // e
    this->vshift_[1], // f
    -1 + this->vshift_[0] + this->vshift_[1], // g
    this->vshift_[0] + this->vshift_[1], // h
    -this->vshift_[0], // V(d)::{b}
    1 - this->vshift_[0], // V(c)::{b}
    1, // V(a)::{b}
  };
  // V(efhg)=V(g)+V(h)::{g}+V(f)::{g,h}
  this->vertexNeighborEFGH_ = {
    -this->vshift_[0] - this->vshift_[1], // a
    1 - this->vshift_[0] - this->vshift_[1], // b
    -this->vshift_[1], // c
    1 - this->vshift_[1], // d
    -this->vshift_[0], // e
    1 - this->vshift_[0], // f
    1, // h
    -1, // V(h)::{g}
    -1 + this->vshift_[0], // V(f)::{g}
    this->vshift_[0], // V(f)::{h}
  };
  // V(aefb)=V(b)+V(a)::{b}+V(e)::{b}+V(f)::{b}
  this->vertexNeighborAEFB_ = {
    -1, // a
    -1 + this->vshift_[0], // c
    this->vshift_[0], // d
    -1 + this->vshift_[1], // e
    this->vshift_[1], // f
    -1 + this->vshift_[0] + this->vshift_[1], // g
    this->vshift_[0] + this->vshift_[1], // h
    1, // V(a)::{b}
    1 - this->vshift_[1], // V(e)::{b}
    -this->vshift_[1], // V(f)::{b}
  };
  // V(ghdc)=V(g)+V(h)::{g}+V(d)::{g,h}
  this->vertexNeighborGHDC_ = {
    -this->vshift_[0] - this->vshift_[1], // a
    1 - this->vshift_[0] - this->vshift_[1], // b
    -this->vshift_[1], // c
    1 - this->vshift_[1], // d
    -this->vshift_[0], // e
    1 - this->vshift_[0], // f
    1, // h
    -1, // V(h)::{g}
    -1 + this->vshift_[1], // V(d)::{g}
    this->vshift_[1] // V(d)::{h}
  };
  // V(aegc)=V(g)+V(a)::{c,g}+V(c)::{g}
  this->vertexNeighborAEGC_ = {
    -this->vshift_[0] - this->vshift_[1], // a
    1 - this->vshift_[0] - this->vshift_[1], // b
    -this->vshift_[1], // c
    1 - this->vshift_[1], // d
    -this->vshift_[0], // e
    1 - this->vshift_[0], // f
    1, // h
    this->vshift_[0], // V(a)::{c}
    this->vshift_[0] + this->vshift_[1], // V(a)::{g}
    this->vshift_[1], // V(c)::{g}
  };
  // V(bfhd)=V(b)+V(f)::{b}+V(h)::{b}+V(d)::{b}
  this->vertexNeighborBFHD_ = {
    -1, // a
    -1 + this->vshift_[0], // c
    this->vshift_[0], // d
    -1 + this->vshift_[1], // e
    this->vshift_[1], // f
    -1 + this->vshift_[0] + this->vshift_[1], // g
    this->vshift_[0] + this->vshift_[1], // h
    -this->vshift_[1], // V(f)::{b}
    -this->vshift_[0] - this->vshift_[1], // V(h)::{b}
    -this->vshift_[0], // V(d)::{b}
  };

  // V(ab)=V(b)+V(a)::{b}
  this->vertexNeighborAB_ = {
    -1, // a
    -1 + this->vshift_[0], // c
    this->vshift_[0], // d
    -1 + this->vshift_[1], // e
    this->vshift_[1], // f
    -1 + this->vshift_[0] + this->vshift_[1], // g
    this->vshift_[0] + this->vshift_[1], // h
    1, // V(a)::{b}
  };
  // V(bd)=V(b)+V(d)::{b}
  this->vertexNeighborBD_ = {
    -1, // a
    -1 + this->vshift_[0], // c
    this->vshift_[0], // d
    -1 + this->vshift_[1], // e
    this->vshift_[1], // f
    -1 + this->vshift_[0] + this->vshift_[1], // g
    this->vshift_[0] + this->vshift_[1], // h
    -this->vshift_[0], // V(d)::{b}
  };
  // V(gh)=V(g)+V(h)::{g}
  this->vertexNeighborGH_ = {
    -this->vshift_[0] - this->vshift_[1], // a
    1 - this->vshift_[0] - this->vshift_[1], // b
    -this->vshift_[1], // c
    1 - this->vshift_[1], // d
    -this->vshift_[0], // e
    1 - this->vshift_[0], // f
    1, // h
    -1, // V(h)::{g}
  };
  // V(eg)=V(g)+V(e)::{g}
  this->vertexNeighborEG_ = {
    -this->vshift_[0] - this->vshift_[1], // a
    1 - this->vshift_[0] - this->vshift_[1], // b
    -this->vshift_[1], // c
    1 - this->vshift_[1], // d
    -this->vshift_[0], // e
    1 - this->vshift_[0], // f
    1, // h
    this->vshift_[0], // V(e)::{g}
  };
  // V(cg)=V(g)+V(c)::{g}
  this->vertexNeighborCG_ = {
    -this->vshift_[0] - this->vshift_[1], // a
    1 - this->vshift_[0] - this->vshift_[1], // b
    -this->vshift_[1], // c
    1 - this->vshift_[1], // d
    -this->vshift_[0], // e
    1 - this->vshift_[0], // f
    1, // h
    this->vshift_[1], // V(c)::{g}
  };
  // V(bf)=V(b)+V(f)::{b}
  this->vertexNeighborBF_ = {
    -1, // a
    -1 + this->vshift_[0], // c
    this->vshift_[0], // d
    -1 + this->vshift_[1], // e
    this->vshift_[1], // f
    -1 + this->vshift_[0] + this->vshift_[1], // g
    this->vshift_[0] + this->vshift_[1], // h
    -this->vshift_[1], // V(f)::{b}
  };

  // V(b)={a,c,d,e,f,g,h}
  this->vertexNeighborB_ = {
    -1, // a
    -1 + this->vshift_[0], // c
    this->vshift_[0], // d
    -1 + this->vshift_[1], // e
    this->vshift_[1], // f
    -1 + this->vshift_[0] + this->vshift_[1], // g
    this->vshift_[0] + this->vshift_[1], // h
  };
  // V(g)={a,b,c,d,e,f,h}
  this->vertexNeighborG_ = {
    -this->vshift_[0] - this->vshift_[1], // a
    1 - this->vshift_[0] - this->vshift_[1], // b
    -this->vshift_[1], // c
    1 - this->vshift_[1], // d
    -this->vshift_[0], // e
    1 - this->vshift_[0], // f
    1, // h
  };

  // V(ef)=V(f)+V(e)::{b,f}
  this->vertexNeighborEF_ = {
    -this->vshift_[1], // b
    -1, // e
    -1 + this->vshift_[0], // g
    this->vshift_[0], // h
    1 - this->vshift_[1], // V(e)::{b}
    1, // V(e)::{f}
  };
  // V(cd)=V(d)+V(c)::{b,d}
  this->vertexNeighborCD_ = {
    -this->vshift_[0], // b
    -1, // c
    -1 + this->vshift_[1], // g
    this->vshift_[1], // h
    1 - this->vshift_[0], // V(c)::{b}
    1, // V(c)::{d}
  };
  // V(ac)=V(c)+V(a)::{c,g}
  this->vertexNeighborAC_ = {
    -this->vshift_[0], // a
    1 - this->vshift_[0], // b
    1, // d
    this->vshift_[1], // g
    this->vshift_[0], // V(a)::{c}
    this->vshift_[0] + this->vshift_[1], // V(a)::{c}
  };
  // V(ae)=V(a)+V(e)::{a,b}
  this->vertexNeighborAE_ = {
    1, // b
    this->vshift_[0], // c
    this->vshift_[1], // e
    this->vshift_[0] + this->vshift_[1], // g
    -this->vshift_[1], // V(e)::{a}
    1 - this->vshift_[1], // V(e)::{b}
  };
  // V(fh)=V(f)+V(h)::{b,f}
  this->vertexNeighborFH_ = {
    -this->vshift_[1], // b
    -1, // e
    -1 + this->vshift_[0], // g
    this->vshift_[0], // h
    -this->vshift_[0] - this->vshift_[1], // V(h)::{b}
    -this->vshift_[0], // V(h)::{f}
  };
  // V(dh)=V(d)+V(h)::{b,d}
  this->vertexNeighborDH_ = {
    -this->vshift_[0], // b
    -1, // c
    -1 + this->vshift_[1], // g
    this->vshift_[1], // h
    -this->vshift_[0] - this->vshift_[1], // V(h)::{b}
    -this->vshift_[1], // V(h)::{d}
  };

  // V(a)={b,c,e,g}
  this->vertexNeighborA_ = {
    1, // b
    this->vshift_[0], // c
    this->vshift_[1], // e
    this->vshift_[0] + this->vshift_[1], // g
  };
  // V(c)={a,b,d,g}
  this->vertexNeighborC_ = {
    -this->vshift_[0], // a
    1 - this->vshift_[0], // b
    1, // d
    +this->vshift_[1], // g
  };
  // V(d)={b,c,g,h}
  this->vertexNeighborD_ = {
    -this->vshift_[0], // b
    -1, // c
    -1 + this->vshift_[1], // g
    this->vshift_[1], // h
  };
  // V(e)={a,b,f,g}
  this->vertexNeighborE_ = {
    -this->vshift_[1], // a
    1 - this->vshift_[1], // b
    1, // f
    +this->vshift_[0], // g
  };
  // V(f)={b,e,g,h}
  this->vertexNeighborF_ = {
    -this->vshift_[1], // b
    -1, // e
    -1 + this->vshift_[0], // g
    this->vshift_[0], // h
  };
  // V(h)={b,d,f,g}
  this->vertexNeighborH_ = {
    -this->vshift_[0] - this->vshift_[1], // b
    -this->vshift_[1], // d
    -this->vshift_[0], // f
    -1, // g
  };

  // V(abcd)=V(d)::{b,c}+V(c)::{b,d}+V(a)::{b}+V(b)::{c}
  this->vertexNeighbor2dABCD_ = {
    -1, -this->vshift_[0], -this->vshift_[0] + 1, 1, this->vshift_[0], this->vshift_[0] - 1,
  };
  // V(ab)=V(b)::{a,c,d}+V(a)::{b}
  this->vertexNeighbor2dAB_ = {
    -1, // V(b)::a
    this->vshift_[0] - 1, // V(b)::c
    this->vshift_[0], // V(b)::d
    +1, // V(a)::b
  };
  // V(cd)=V(c)::{a,b,d}+V(d)::{c}
  this->vertexNeighbor2dCD_ = {
    -1, // V(d)::c
    -this->vshift_[0], // V(c)::a
    -this->vshift_[0] + 1, // V(c)::b
    1, // V(c)::d
  };
  // V(ac)=V(c)::{a,b,d}+V(a)::{c}
  this->vertexNeighbor2dAC_ = {
    -this->vshift_[0], // V(c)::{a}
    -this->vshift_[0] + 1, // V(c)::{b}
    1, // V(c)::{d}
    this->vshift_[0], // V(a)::{c}
  };
  // V(bd)=V(b)::{c,d}+V(d)::{b,c}
  this->vertexNeighbor2dBD_ = {
    this->vshift_[0] - 1, // V(b)::{c}
    this->vshift_[0], // V(b)::{d}
    -this->vshift_[0], // V(d)::{b}
    -1, // V(d)::{c}
  };
  // V(b)={a,c,d}
  this->vertexNeighbor2dB_ = {
    -1, // a
    this->vshift_[0], // d
    this->vshift_[0] - 1, // c
  };
  // V(c)={a,b,d}
  this->vertexNeighbor2dC_ = {
    1, // d
    -this->vshift_[0], // a
    -this->vshift_[0] + 1, // b
  };
  this->vertexNeighbor2dA_ = {};
  // V(a)={b,c}
  this->vertexNeighbor2dA_ = {
    1, // b
    this->vshift_[0] // c
  };
  // V(d)={c,b}
  this->vertexNeighbor2dD_ = {
    -1, // c
    -this->vshift_[0], // b

  };

  return 0;
}

#ifdef TTK_ENABLE_MPI
template <size_t card>
int ttk::ImplicitTriangulation<card>::preconditionDistributedCells() {
  if(this->hasPreconditionedDistributedCells_) {
    return 0;
  }
  if(!ttk::hasInitializedMPI()) {
    return -1;
  }
  if(this->cellGhost_ == nullptr) {
    if(ttk::isRunningWithMPI()) {
      this->printErr("Missing cell ghost array!");
    }
    return -3;
  }

  Timer tm{};

  // number of local cells (with ghost cells...)
  const auto nLocCells{this->getNumberOfCells()};

  // there are 6 tetrahedra per cubic cell (and 2 triangles per square)
  const int nTetraPerCube{this->dimensionality_ == 3 ? 6 : 2};
  std::vector<unsigned char> fillCells(nLocCells / nTetraPerCube);

  this->ghostCellsPerOwner_.resize(ttk::MPIsize_);

  this->neighborCellBBoxes_.resize(ttk::MPIsize_);
  auto &localBBox{this->neighborCellBBoxes_[ttk::MPIrank_]};
  // "good" starting values?
  localBBox = {
    this->localGridOffset_[0] + this->dimensions_[0], this->localGridOffset_[0],
    this->localGridOffset_[1] + this->dimensions_[1], this->localGridOffset_[1],
    this->localGridOffset_[2] + this->dimensions_[2], this->localGridOffset_[2],
  };

  for(SimplexId lcid = 0; lcid < nLocCells; ++lcid) {
    // only keep non-ghost cells
    if(this->cellGhost_[lcid / nTetraPerCube] == 1) {
      continue;
    }
    // local vertex coordinates
    std::array<SimplexId, 3> p{};
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
      localBBox[0] = p[0];
    }
    if(p[0] > localBBox[1]) {
      localBBox[1] = p[0];
    }
    if(p[1] < localBBox[2]) {
      localBBox[2] = p[1];
    }
    if(p[1] > localBBox[3]) {
      localBBox[3] = p[1];
    }
    if(p[2] < localBBox[4]) {
      localBBox[4] = p[2];
    }
    if(p[2] > localBBox[5]) {
      localBBox[5] = p[2];
    }
  }
  localBBox[1]++;
  localBBox[3]++;
  localBBox[5]++;

  for(size_t i = 0; i < this->neighborRanks_.size(); ++i) {
    const auto neigh{this->neighborRanks_[i]};
    MPI_Sendrecv(this->neighborCellBBoxes_[ttk::MPIrank_].data(), 6,
                 ttk::getMPIType(SimplexId{}), neigh, ttk::MPIrank_,
                 this->neighborCellBBoxes_[neigh].data(), 6,
                 ttk::getMPIType(SimplexId{}), neigh, neigh, ttk::MPIcomm_,
                 MPI_STATUS_IGNORE);
  }

  this->hasPreconditionedDistributedCells_ = true;

  return 0;
}


template <size_t card>
void ttk::ImplicitTriangulation<card>::createMetaGrid(const double *const bounds) {

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

  const std::array<int, 3> dimensions = {
    static_cast<int>(
      std::round((globalBounds[1] - globalBounds[0]) / this->spacing_[0]))
      + 1,
    static_cast<int>(
      std::round((globalBounds[3] - globalBounds[2]) / this->spacing_[1]))
      + 1,
    static_cast<int>(
      std::round((globalBounds[5] - globalBounds[4]) / this->spacing_[2]))
      + 1,
  };

  this->localGridOffset_ = {
    static_cast<SimplexId>(
      std::round((this->origin_[0] - globalBounds[0]) / this->spacing_[0])),
    static_cast<SimplexId>(
      std::round((this->origin_[1] - globalBounds[2]) / this->spacing_[1])),
    static_cast<SimplexId>(
      std::round((this->origin_[2] - globalBounds[4]) / this->spacing_[2])),
  };

  this->metaGrid_ = std::make_shared<ImplicitNoPreconditions<card> >();
  this->metaGrid_->setInputGrid(globalBounds[0], globalBounds[1],
                                globalBounds[2], this->spacing_[0],
                                this->spacing_[1], this->spacing_[2],
                                dimensions[0], dimensions[1], dimensions[2]);
  this->metaGrid_->preconditionBoundaryVertices();
  this->metaGrid_->preconditionBoundaryEdges();
  this->metaGrid_->preconditionBoundaryTriangles();
}

template <size_t card>
std::array<SimplexId, 3>
  ttk::ImplicitTriangulation<card>::getVertGlobalCoords(const SimplexId lvid) const {

  // local vertex coordinates
  std::array<SimplexId, 3> p{};
  if(this->dimensionality_ == 3) {
    this->vertexToPosition(lvid, p.data());
  } else if(this->dimensionality_ == 2) {
    this->vertexToPosition2d(lvid, p.data());
  }

  // global vertex coordinates
  p[0] += this->localGridOffset_[0];
  p[1] += this->localGridOffset_[1];
  p[2] += this->localGridOffset_[2];

  return p;
}

template <size_t card>
std::array<SimplexId, 3>
  ttk::ImplicitTriangulation<card>::getVertLocalCoords(const SimplexId gvid) const {

  // global vertex coordinates
  std::array<SimplexId, 3> p{};
  if(this->dimensionality_ == 3) {
    this->metaGrid_->vertexToPosition(gvid, p.data());
  } else if(this->dimensionality_ == 2) {
    this->metaGrid_->vertexToPosition2d(gvid, p.data());
  }

  // local vertex coordinates
  p[0] -= this->localGridOffset_[0];
  p[1] -= this->localGridOffset_[1];
  p[2] -= this->localGridOffset_[2];

  return p;
}

template<size_t card>
bool ImplicitTriangulation<card>::isVertexOnGlobalBoundaryInternal(
  const SimplexId lvid) const {

  if(!ttk::isRunningWithMPI()) {
    return this->isVertexOnBoundary(lvid);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(lvid > this->TTK_TRIANGULATION_INTERNAL(getNumberOfVertices)() - 1
     || lvid < 0) {
    return false;
  }
  if(this->metaGrid_ == nullptr) {
    return false;
  }
#endif // TTK_ENABLE_KAMIKAZE

  const auto gvid{this->getVertexGlobalIdInternal(lvid)};
  if(gvid == -1) {
    return false;
  }
  return this->metaGrid_->isVertexOnBoundary(gvid);
}

template<size_t card>
bool ImplicitTriangulation<card>::isEdgeOnGlobalBoundaryInternal(
  const SimplexId leid) const {

  if(!ttk::isRunningWithMPI()) {
    return this->isEdgeOnBoundary(leid);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(leid > this->getNumberOfEdgesInternal() - 1 || leid < 0) {
    return false;
  }
  if(this->metaGrid_ == nullptr) {
    return false;
  }
#endif // TTK_ENABLE_KAMIKAZE

  const auto geid{this->getEdgeGlobalIdInternal(leid)};
  if(geid == -1) {
    return false;
  }
  return this->metaGrid_->isEdgeOnBoundary(geid);
}

template<size_t card>
bool ImplicitTriangulation<card>::isTriangleOnGlobalBoundaryInternal(
  const SimplexId ltid) const {

  if(!ttk::isRunningWithMPI()) {
    return this->isTriangleOnBoundary(ltid);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(ltid > this->getNumberOfTrianglesInternal() - 1 || ltid < 0) {
    return false;
  }
  if(this->metaGrid_ == nullptr) {
    return false;
  }
#endif // TTK_ENABLE_KAMIKAZE

  const auto gtid{this->getTriangleGlobalIdInternal(ltid)};
  if(gtid == -1) {
    return false;
  }
  return this->metaGrid_->isTriangleOnBoundary(gtid);
}

template <size_t card>
int ttk::ImplicitTriangulation<card>::getCellRankInternal(
  const SimplexId lcid) const {

  const int nTetraPerCube{this->dimensionality_ == 3 ? 6 : 2};
  const auto locCubeId{lcid / nTetraPerCube};

  if(this->cellGhost_[locCubeId] == 0) {
    return ttk::MPIrank_;
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(this->neighborRanks_.empty() && ttk::MPIsize_ > 1) {
    this->printErr("Empty neighborsRanks_!");
    return -1;
  }
#endif // TTK_ENABLE_KAMIKAZE

  const auto nVertsCell{this->getCellVertexNumber(lcid)};
  std::vector<bool> inRank(nVertsCell);
  for(const auto neigh : this->neighborRanks_) {
    std::fill(inRank.begin(), inRank.end(), false);
    const auto &bbox{this->neighborCellBBoxes_[neigh]};
    for(SimplexId i = 0; i < nVertsCell; ++i) {
      SimplexId v{};
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
  }

  return -1;
}

#endif // TTK_ENABLE_MPI

// explicit instantiations
template class ttk::ImplicitTriangulation<0>;
template class ttk::ImplicitTriangulation<1>;
template class ttk::ImplicitTriangulation<2>;
template class ttk::ImplicitTriangulation<3>;
template class ttk::ImplicitTriangulationCRTP<0, ttk::ImplicitWithPreconditions<0>>;
template class ttk::ImplicitTriangulationCRTP<0, ttk::ImplicitNoPreconditions<0>>;
template class ttk::ImplicitTriangulationCRTP<1, ttk::ImplicitWithPreconditions<1>>;
template class ttk::ImplicitTriangulationCRTP<1, ttk::ImplicitNoPreconditions<1>>;
template class ttk::ImplicitTriangulationCRTP<2, ttk::ImplicitWithPreconditions<2>>;
template class ttk::ImplicitTriangulationCRTP<2, ttk::ImplicitNoPreconditions<2>>;
template class ttk::ImplicitTriangulationCRTP<3, ttk::ImplicitWithPreconditions<3>>;
template class ttk::ImplicitTriangulationCRTP<3, ttk::ImplicitNoPreconditions<3>>;
/// @endcond
