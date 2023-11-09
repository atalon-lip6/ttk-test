#pragma once

#include <PeriodicImplicitTriangulation.h>

namespace ttk {
  /**
   * @brief Periodic implicit Triangulation class with preconditioning
   */
  template<size_t card>
  class PeriodicWithPreconditions final
    : public PeriodicImplicitTriangulationCRTP<card, PeriodicWithPreconditions<card>> {
  public:
    PeriodicWithPreconditions() {
      this->setDebugMsgPrefix("PeriodicTriangulationWithPreconditions");
    }

    using EdgePosition = typename PeriodicWithPreconditions<card>::EdgePosition; 
    using TrianglePosition = typename PeriodicWithPreconditions<card>::TrianglePosition; 
    int preconditionVerticesInternal() override;
    int preconditionEdgesInternal() override;
    int preconditionTrianglesInternal() override;
    int preconditionTetrahedronsInternal() override;

    inline std::array<SimplexId, 3> const &
      getVertexCoords(const SimplexId v) const {
      return this->vertexCoords_[v];
    }
    inline EdgePosition getEdgePosition(const SimplexId e) const {
      return this->edgePositions_[e];
    }
    inline std::array<SimplexId, 3> const &
      getEdgeCoords(const SimplexId e) const {
      return this->edgeCoords_[e];
    }
    inline TrianglePosition getTrianglePosition(const SimplexId t) const {
      return this->trianglePositions_[t];
    }
    inline std::array<SimplexId, 3> const &
      getTriangleCoords(const SimplexId t) const {
      return this->triangleCoords_[t];
    }
    inline std::array<SimplexId, 3> const &
      getTetrahedronCoords(const SimplexId t) const {
      return this->tetrahedronCoords_[t];
    }
    inline SimplexId getEdgeVertexAccelerated(const SimplexId e) const {
      return this->edgeVertexAccelerated_[e];
    }

    inline int clear() {
      vertexCoords_ = std::vector<std::array<SimplexId, 3>>{};
      edgePositions_ = std::vector<EdgePosition>{};
      edgeCoords_ = std::vector<std::array<SimplexId, 3>>{};
      trianglePositions_ = std::vector<TrianglePosition>{};
      triangleCoords_ = std::vector<std::array<SimplexId, 3>>{};
      tetrahedronCoords_ = std::vector<std::array<SimplexId, 3>>{};
      edgeVertexAccelerated_ = std::vector<SimplexId>{};
      this->hasPreconditionedVerticesAndCells_ = false;
      return AbstractTriangulation::clear();
    }

  private:
    // for  every vertex, its coordinates on the grid
    std::vector<std::array<SimplexId, 3>> vertexCoords_{};
    // for every edge, its position on the grid
    std::vector<EdgePosition> edgePositions_{};
    // for every edge, its coordinates on the grid
    std::vector<std::array<SimplexId, 3>> edgeCoords_{};
    // for every triangle, its position on the grid
    std::vector<TrianglePosition> trianglePositions_{};
    // for every triangle, its coordinates on the grid
    std::vector<std::array<SimplexId, 3>> triangleCoords_{};
    // for every tetrahedron, its coordinates on the grid
    std::vector<std::array<SimplexId, 3>> tetrahedronCoords_{};

    // cache some edge vertex computation wrt acceleration
    std::vector<SimplexId> edgeVertexAccelerated_{};
  };

  /**
   * @brief Periodic implicit Triangulation class without preconditioning
   */
  template <size_t card>
  class PeriodicNoPreconditions final
    : public PeriodicImplicitTriangulationCRTP<0, PeriodicNoPreconditions<card>> {
  public:
    PeriodicNoPreconditions() {
      this->setDebugMsgPrefix("PeriodicTriangulationNoPreconditions");
    }

    using EdgePosition = typename PeriodicNoPreconditions<card>::EdgePosition; 
    using TrianglePosition = typename PeriodicNoPreconditions<card>::TrianglePosition; 
    inline int preconditionVerticesInternal() override {
      return 0;
    }
    inline int preconditionEdgesInternal() override {
      return 0;
    }
    inline int preconditionTrianglesInternal() override {
      return 0;
    }
    inline int preconditionTetrahedronsInternal() override {
      return 0;
    }

    inline std::array<SimplexId, 3> getVertexCoords(const SimplexId v) const {
      std::array<SimplexId, 3> p{};
      if(this->dimensionality_ == 2) {
        this->vertexToPosition2d(v, p.data());
      } else if(this->dimensionality_ == 3) {
        this->vertexToPosition(v, p.data());
      }
      return p;
    }

    EdgePosition getEdgePosition(const SimplexId e) const;
    std::array<SimplexId, 3> getEdgeCoords(const SimplexId e) const;
    TrianglePosition getTrianglePosition(const SimplexId t) const;
    std::array<SimplexId, 3> getTriangleCoords(const SimplexId t) const;

    inline std::array<SimplexId, 3>
      getTetrahedronCoords(const SimplexId t) const {
      std::array<SimplexId, 3> p{};
      this->tetrahedronToPosition(t, p.data());
      return p;
    }

    inline SimplexId getEdgeVertexAccelerated(const SimplexId e) const {
      const auto p{this->getEdgeCoords(e)};
      if(this->isAccelerated_) {
        return (p[1] << this->div_[0]) + (p[2] << this->div_[1]);
      } else {
        return p[1] * this->vshift_[0] + p[2] * this->vshift_[1];
      }
    }
  };

// -- was in cpp

template<size_t card>
int ttk::PeriodicWithPreconditions<card>::preconditionVerticesInternal() {
  vertexCoords_.resize(this->vertexNumber_);

  if(this->dimensionality_ == 1) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->vertexNumber_; ++i) {
      vertexCoords_[i][0] = i;
    }
  } else if(this->dimensionality_ == 2) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->vertexNumber_; ++i) {
      auto &p = vertexCoords_[i];
      this->vertexToPosition2d(i, p.data());
    }
  } else if(this->dimensionality_ == 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->vertexNumber_; ++i) {
      auto &p = vertexCoords_[i];
      this->vertexToPosition(i, p.data());
    }
  }

  return 0;
}

template<size_t card>
int ttk::PeriodicWithPreconditions<card>::preconditionEdgesInternal() {
  edgePositions_.resize(this->edgeNumber_);
  edgeCoords_.resize(this->edgeNumber_);

  if(this->dimensionality_ == 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->edgeNumber_; ++i) {
      auto &p = edgeCoords_[i];
      if(i < this->esetshift_[0]) {
        this->edgeToPosition(i, 0, p.data());
        edgePositions_[i] = EdgePosition::L_3D;
      } else if(i < this->esetshift_[1]) {
        this->edgeToPosition(i, 1, p.data());
        edgePositions_[i] = EdgePosition::H_3D;
      } else if(i < this->esetshift_[2]) {
        this->edgeToPosition(i, 2, p.data());
        edgePositions_[i] = EdgePosition::P_3D;
      } else if(i < this->esetshift_[3]) {
        this->edgeToPosition(i, 3, p.data());
        edgePositions_[i] = EdgePosition::D1_3D;
      } else if(i < this->esetshift_[4]) {
        this->edgeToPosition(i, 4, p.data());
        edgePositions_[i] = EdgePosition::D2_3D;
      } else if(i < this->esetshift_[5]) {
        this->edgeToPosition(i, 5, p.data());
        edgePositions_[i] = EdgePosition::D3_3D;
      } else if(i < this->esetshift_[6]) {
        this->edgeToPosition(i, 6, p.data());
        edgePositions_[i] = EdgePosition::D4_3D;
      }
    }

  } else if(this->dimensionality_ == 2) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->edgeNumber_; ++i) {
      auto &p = edgeCoords_[i];
      if(i < this->esetshift_[0]) {
        this->edgeToPosition2d(i, 0, p.data());
        edgePositions_[i] = EdgePosition::L_2D;
      } else if(i < this->esetshift_[1]) {
        this->edgeToPosition2d(i, 1, p.data());
        edgePositions_[i] = EdgePosition::H_2D;
      } else if(i < this->esetshift_[2]) {
        this->edgeToPosition2d(i, 2, p.data());
        edgePositions_[i] = EdgePosition::D1_2D;
      }
    }
  }

  else if(this->dimensionality_ == 1) {
    edgePositions_[0] = EdgePosition::FIRST_EDGE_1D;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 1; i < this->edgeNumber_ - 1; ++i) {
      edgePositions_[i] = EdgePosition::CENTER_1D;
    }
    edgePositions_[this->edgeNumber_ - 1] = EdgePosition::LAST_EDGE_1D;
  }

  edgeVertexAccelerated_.resize(this->edgeNumber_);

  if(this->isAccelerated_) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->edgeNumber_; ++i) {
      const auto &p = edgeCoords_[i];
      edgeVertexAccelerated_[i] = (p[1] << this->div_[0]) + (p[2] << this->div_[1]);
    }
  } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->edgeNumber_; ++i) {
      const auto &p = edgeCoords_[i];
      edgeVertexAccelerated_[i] = p[1] * this->vshift_[0] + p[2] * this->vshift_[1];
    }
  }

  return 0;
}

template<size_t card>
int ttk::PeriodicWithPreconditions<card>::preconditionTrianglesInternal() {
  if(this->dimensionality_ != 3 && this->dimensionality_ != 2) {
    return 1;
  }

  trianglePositions_.resize(this->triangleNumber_);
  triangleCoords_.resize(this->triangleNumber_);

  if(this->dimensionality_ == 3) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->triangleNumber_; ++i) {
      if(i < this->tsetshift_[0]) {
        this->triangleToPosition(i, 0, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::F_3D;
      } else if(i < this->tsetshift_[1]) {
        this->triangleToPosition(i, 1, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::H_3D;
      } else if(i < this->tsetshift_[2]) {
        this->triangleToPosition(i, 2, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::C_3D;
      } else if(i < this->tsetshift_[3]) {
        this->triangleToPosition(i, 3, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::D1_3D;
      } else if(i < this->tsetshift_[4]) {
        this->triangleToPosition(i, 4, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::D2_3D;
      } else if(i < this->tsetshift_[5]) {
        this->triangleToPosition(i, 5, triangleCoords_[i].data());
        trianglePositions_[i] = TrianglePosition::D3_3D;
      }
    }

  } else if(this->dimensionality_ == 2) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
    for(SimplexId i = 0; i < this->triangleNumber_; ++i) {
      this->triangleToPosition2d(i, triangleCoords_[i].data());
      if(i % 2 == 0) {
        trianglePositions_[i] = TrianglePosition::TOP_2D;
      } else {
        trianglePositions_[i] = TrianglePosition::BOTTOM_2D;
      }
    }
  }

  return 0;
}


template<size_t card>
int ttk::PeriodicWithPreconditions<card>::preconditionTetrahedronsInternal() {

  if(this->dimensionality_ != 3) {
    return 1;
  }
  tetrahedronCoords_.resize(this->tetrahedronNumber_);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i = 0; i < this->tetrahedronNumber_; ++i) {
    this->tetrahedronToPosition(i, tetrahedronCoords_[i].data());
  }
  return 0;
}

template<size_t card>
typename ttk::PeriodicNoPreconditions<card>::EdgePosition
  ttk::PeriodicNoPreconditions<card>::getEdgePosition(const SimplexId e) const {

  std::array<SimplexId, 3> p{};
  if(this->dimensionality_ == 3) {
    if(e < this->esetshift_[0]) {
      this->edgeToPosition(e, 0, p.data());
      return EdgePosition::L_3D;
    } else if(e < this->esetshift_[1]) {
      this->edgeToPosition(e, 1, p.data());
      return EdgePosition::H_3D;
    } else if(e < this->esetshift_[2]) {
      this->edgeToPosition(e, 2, p.data());
      return EdgePosition::P_3D;
    } else if(e < this->esetshift_[3]) {
      this->edgeToPosition(e, 3, p.data());
      return EdgePosition::D1_3D;
    } else if(e < this->esetshift_[4]) {
      this->edgeToPosition(e, 4, p.data());
      return EdgePosition::D2_3D;
    } else if(e < this->esetshift_[5]) {
      this->edgeToPosition(e, 5, p.data());
      return EdgePosition::D3_3D;
    } else if(e < this->esetshift_[6]) {
      this->edgeToPosition(e, 6, p.data());
      return EdgePosition::D4_3D;
    }
  } else if(this->dimensionality_ == 2) {
    if(e < this->esetshift_[0]) {
      this->edgeToPosition2d(e, 0, p.data());
      return EdgePosition::L_2D;
    } else if(e < this->esetshift_[1]) {
      this->edgeToPosition2d(e, 1, p.data());
      return EdgePosition::H_2D;
    } else if(e < this->esetshift_[2]) {
      this->edgeToPosition2d(e, 2, p.data());
      return EdgePosition::D1_2D;
    }
  } else if(this->dimensionality_ == 1) {
    if(e == 0) {
      return EdgePosition::FIRST_EDGE_1D;
    } else if(e == this->edgeNumber_ - 1)
      return EdgePosition::CENTER_1D;
  } else {
    return EdgePosition::LAST_EDGE_1D;
  }

  return EdgePosition::CENTER_1D;
}

template<size_t card>
std::array<ttk::SimplexId, 3>
  ttk::PeriodicNoPreconditions<card>::getEdgeCoords(const SimplexId e) const {
  std::array<SimplexId, 3> p{};
  if(this->dimensionality_ == 3) {
    if(e < this->esetshift_[0]) {
      this->edgeToPosition(e, 0, p.data());
    } else if(e < this->esetshift_[1]) {
      this->edgeToPosition(e, 1, p.data());
    } else if(e < this->esetshift_[2]) {
      this->edgeToPosition(e, 2, p.data());
    } else if(e < this->esetshift_[3]) {
      this->edgeToPosition(e, 3, p.data());
    } else if(e < this->esetshift_[4]) {
      this->edgeToPosition(e, 4, p.data());
    } else if(e < this->esetshift_[5]) {
      this->edgeToPosition(e, 5, p.data());
    } else if(e < this->esetshift_[6]) {
      this->edgeToPosition(e, 6, p.data());
    }

  } else if(this->dimensionality_ == 2) {
    if(e < this->esetshift_[0]) {
      this->edgeToPosition2d(e, 0, p.data());
    } else if(e < this->esetshift_[1]) {
      this->edgeToPosition2d(e, 1, p.data());
    } else if(e < this->esetshift_[2]) {
      this->edgeToPosition2d(e, 2, p.data());
    }
  }
  return p;
}

template<size_t card>
typename ttk::PeriodicNoPreconditions<card>::TrianglePosition
  ttk::PeriodicNoPreconditions<card>::getTrianglePosition(const SimplexId t) const {
  if(this->dimensionality_ == 2) {
    if(t % 2 == 0) {
      return TrianglePosition::TOP_2D;
    } else {
      return TrianglePosition::BOTTOM_2D;
    }
  } else if(this->dimensionality_ == 3) {
    if(t < this->tsetshift_[0]) {
      return TrianglePosition::F_3D;
    } else if(t < this->tsetshift_[1]) {
      return TrianglePosition::H_3D;
    } else if(t < this->tsetshift_[2]) {
      return TrianglePosition::C_3D;
    } else if(t < this->tsetshift_[3]) {
      return TrianglePosition::D1_3D;
    } else if(t < this->tsetshift_[4]) {
      return TrianglePosition::D2_3D;
    } else if(t < this->tsetshift_[5]) {
      return TrianglePosition::D3_3D;
    }
  }
  return TrianglePosition::C_3D;
}

template<size_t card>
std::array<ttk::SimplexId, 3>
  ttk::PeriodicNoPreconditions<card>::getTriangleCoords(const SimplexId t) const {
  std::array<SimplexId, 3> p{};
  if(this->dimensionality_ == 2) {
    this->triangleToPosition2d(t, p.data());
  } else if(this->dimensionality_ == 3) {
    if(t < this->tsetshift_[0]) {
      this->triangleToPosition(t, 0, p.data());
    } else if(t < this->tsetshift_[1]) {
      this->triangleToPosition(t, 1, p.data());
    } else if(t < this->tsetshift_[2]) {
      this->triangleToPosition(t, 2, p.data());
    } else if(t < this->tsetshift_[3]) {
      this->triangleToPosition(t, 3, p.data());
    } else if(t < this->tsetshift_[4]) {
      this->triangleToPosition(t, 4, p.data());
    } else if(t < this->tsetshift_[5]) {
      this->triangleToPosition(t, 5, p.data());
    }
  }
  return p;
}

} // namespace ttk
