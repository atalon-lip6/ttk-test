/// \ingroup base
/// \class ttk::PeriodicImplicitTriangulation
/// \author Talha Bin Masood <talha.bin.masood@liu.se>
/// \date September 2019
///
/// \brief TTK triangulation class for grids with periodic boundary conditions
/// implemented in all directions.
///
/// \sa ttk::Triangulation
/// \sa ttk::Triangulation::setPeriodicBoundaryConditions
///

#pragma once

// base code includes
#include <RegularGridTriangulation.h>

#include <array>
#include <numeric>

namespace ttk {

  template<size_t card>
  class PeriodicImplicitTriangulation : public RegularGridTriangulation<card> {

  public:
    PeriodicImplicitTriangulation();
    ~PeriodicImplicitTriangulation() override;

    PeriodicImplicitTriangulation(const PeriodicImplicitTriangulation &)
      = default;
    PeriodicImplicitTriangulation(PeriodicImplicitTriangulation &&) = default;
    PeriodicImplicitTriangulation &
      operator=(const PeriodicImplicitTriangulation &)
      = default;
    PeriodicImplicitTriangulation &operator=(PeriodicImplicitTriangulation &&)
      = default;

    int getCellEdgeInternal(const SimplexId &cellId,
                            const int &id,
                            SimplexId &edgeId) const override;

    SimplexId getCellEdgeNumberInternal(const SimplexId &cellId) const override;

    const std::vector<std::vector<SimplexId>> *getCellEdgesInternal() override;

    int TTK_TRIANGULATION_INTERNAL(getCellNeighbor)(
      const SimplexId &cellId,
      const int &localNeighborId,
      SimplexId &neighborId) const override;

    SimplexId TTK_TRIANGULATION_INTERNAL(getCellNeighborNumber)(
      const SimplexId &cellId) const override;

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getCellNeighbors)() override;

    int getCellTriangleInternal(const SimplexId &cellId,
                                const int &id,
                                SimplexId &triangleId) const override;

    SimplexId getCellTriangleNumberInternal(
      const SimplexId &ttkNotUsed(cellId)) const override {
      // NOTE: the output is always 4 here. let's keep the function in there
      // in case of further generalization to CW-complexes
      return 4;
    }

    const std::vector<std::vector<SimplexId>> *
      getCellTrianglesInternal() override;

    int TTK_TRIANGULATION_INTERNAL(getCellVertex)(
      const SimplexId &cellId,
      const int &localVertexId,
      SimplexId &vertexId) const override;

    SimplexId TTK_TRIANGULATION_INTERNAL(getCellVertexNumber)(
      const SimplexId &cellId) const override;

    int TTK_TRIANGULATION_INTERNAL(getDimensionality)() const override {
      return dimensionality_;
    }

    SimplexId TTK_TRIANGULATION_INTERNAL(getEdgeLinkNumber)(
      const SimplexId &edgeId) const override;

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeLinks)() override;

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeStars)() override;

    const std::vector<std::vector<SimplexId>> *
      getEdgeTrianglesInternal() override;

    const std::vector<std::array<SimplexId, 2>> *
      TTK_TRIANGULATION_INTERNAL(getEdges)() override;

    SimplexId TTK_TRIANGULATION_INTERNAL(getNumberOfCells)() const override {
      return cellNumber_;
    }

    SimplexId getNumberOfEdgesInternal() const override {
      return edgeNumber_;
    }

    SimplexId getNumberOfTrianglesInternal() const override {
      return triangleNumber_;
    }

    SimplexId TTK_TRIANGULATION_INTERNAL(getNumberOfVertices)() const override {
      return vertexNumber_;
    }

    virtual int getTetrahedronEdge(const SimplexId &tetId,
                                   const int &id,
                                   SimplexId &edgeId) const = 0;

    int getTetrahedronEdges(std::vector<std::vector<SimplexId>> &edges) const;

    virtual int getTetrahedronTriangle(const SimplexId &tetId,
                                       const int &id,
                                       SimplexId &triangleId) const = 0;

    int getTetrahedronTriangles(
      std::vector<std::vector<SimplexId>> &triangles) const;

    virtual int getTetrahedronNeighbor(const SimplexId &tetId,
                                       const int &localNeighborId,
                                       SimplexId &neighborId) const = 0;

    SimplexId getTetrahedronNeighborNumber(const SimplexId &tetId) const;

    int getTetrahedronNeighbors(std::vector<std::vector<SimplexId>> &neighbors);

    virtual int getTetrahedronVertex(const SimplexId &tetId,
                                     const int &localVertexId,
                                     SimplexId &vertexId) const = 0;

    SimplexId getTriangleEdgeNumberInternal(
      const SimplexId &ttkNotUsed(triangleId)) const override {
      // NOTE: the output is always 3 here. let's keep the function in there
      // in case of further generalization to CW-complexes
      return 3;
    }

    const std::vector<std::vector<SimplexId>> *
      getTriangleEdgesInternal() override;

    int getTriangleEdgesInternal(
      std::vector<std::vector<SimplexId>> &edges) const;

    SimplexId TTK_TRIANGULATION_INTERNAL(getTriangleLinkNumber)(
      const SimplexId &triangleId) const override;

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleLinks)() override;

    virtual int getTriangleNeighbor(const SimplexId &triangleId,
                                    const int &localNeighborId,
                                    SimplexId &neighborId) const = 0;

    SimplexId getTriangleNeighborNumber(const SimplexId &triangleId) const;

    int getTriangleNeighbors(std::vector<std::vector<SimplexId>> &neighbors);

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleStars)() override;

    const std::vector<std::array<SimplexId, 3>> *
      TTK_TRIANGULATION_INTERNAL(getTriangles)() override;

    SimplexId
      getVertexEdgeNumberInternal(const SimplexId &vertexId) const override;

    const std::vector<std::vector<SimplexId>> *
      getVertexEdgesInternal() override;

    SimplexId TTK_TRIANGULATION_INTERNAL(getVertexLinkNumber)(
      const SimplexId &vertexId) const override;

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexLinks)() override;

    SimplexId TTK_TRIANGULATION_INTERNAL(getVertexNeighborNumber)(
      const SimplexId &vertexId) const override;

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexNeighbors)() override;

    SimplexId TTK_TRIANGULATION_INTERNAL(getVertexStarNumber)(
      const SimplexId &vertexId) const override;

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexStars)() override;

    SimplexId
      getVertexTriangleNumberInternal(const SimplexId &vertexId) const override;

    const std::vector<std::vector<SimplexId>> *
      getVertexTrianglesInternal() override;

    bool TTK_TRIANGULATION_INTERNAL(isEdgeOnBoundary)(
      const SimplexId &edgeId) const override;

    inline bool isEmpty() const override {
      return !vertexNumber_;
    }

    bool TTK_TRIANGULATION_INTERNAL(isTriangleOnBoundary)(
      const SimplexId &triangleId) const override;

    bool TTK_TRIANGULATION_INTERNAL(isVertexOnBoundary)(
      const SimplexId &vertexId) const override;

    int setInputGrid(const float &xOrigin,
                     const float &yOrigin,
                     const float &zOrigin,
                     const float &xSpacing,
                     const float &ySpacing,
                     const float &zSpacing,
                     const SimplexId &xDim,
                     const SimplexId &yDim,
                     const SimplexId &zDim) override;

    inline const std::array<ttk::SimplexId, 3> &
      getGridDimensions() const override {
      return this->dimensions_;
    }

    virtual int preconditionVerticesInternal() = 0;
    int preconditionEdgesInternal() override = 0;
    int preconditionTrianglesInternal() override = 0;
    virtual int preconditionTetrahedronsInternal() = 0;

    inline int preconditionCellsInternal() {
      if(dimensionality_ == 3) {
        return this->preconditionTetrahedronsInternal();
      } else if(dimensionality_ == 2 && !this->hasPreconditionedTriangles_) {
        this->hasPreconditionedTriangles_ = true;
        return this->preconditionTrianglesInternal();
      }
      return 0;
    }

    inline int preconditionVerticesAndCells() {
      if(!this->hasPreconditionedVerticesAndCells_) {
        this->preconditionVerticesInternal();
        this->preconditionCellsInternal();
        this->hasPreconditionedVerticesAndCells_ = true;
      }
      return 0;
    }

    inline int getCellVTKIDInternal(const int &ttkId,
                                    int &vtkId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(ttkId < 0) {
        return -1;
      }
#endif // TTK_ENABLE_KAMIKAZE
      const SimplexId nSimplexPerCell{this->getDimensionality() == 3 ? 6 : 2};
      vtkId = ttkId / nSimplexPerCell;
      return 0;
    }

  protected:
    int dimensionality_; //
    float origin_[3]; //
    float spacing_[3]; //
    SimplexId nbvoxels_[3]; // nombre de voxels par axe
    SimplexId wrap_[3];

    // Vertex helper //
    SimplexId vshift_[2]; // VertexShift

    // Edge helper //
    SimplexId esetdims_[7]; // EdgeSetDimensions
    SimplexId esetshift_[7]; // EdgeSetShift
    SimplexId eshift_[14]; // EdgeShift

    // Triangle helper //
    SimplexId tsetdims_[6]; // TriangleSetDimensions
    SimplexId tsetshift_[6]; // TriangleSetShift
    SimplexId tshift_[12]; // TriangleShift

    // Tetrahedron helper //
    SimplexId tetshift_[2]; // TetrahedronShift

    SimplexId cellNumber_; // number of cells
    SimplexId vertexNumber_; // number of vertices
    SimplexId edgeNumber_; // number of edges
    SimplexId triangleNumber_; // number of triangles
    SimplexId tetrahedronNumber_; // number of tetrahedra

    // 2d helpers
    SimplexId Di_{};
    SimplexId Dj_{};

    // acceleration variables
    bool isAccelerated_;
    SimplexId mod_[2];
    SimplexId div_[2];

    enum class EdgePosition : char {
      //    e--------f
      //   /|       /|
      //  / |      / |
      // a--------b  |
      // |  g-----|--h
      // | /      | /
      // |/       |/
      // c--------d

      // length (ab)
      L_3D,
      // height (ac)
      H_3D,
      // depth (ae)
      P_3D,
      // diagonal1 (bc)
      D1_3D,
      // diagonal2 (ag)
      D2_3D,
      // diagonal3 (be)
      D3_3D,
      // diagonal4 (bg)
      D4_3D,

      // length (ab)
      L_2D,
      // height (ac)
      H_2D,
      // diagonal1 (bc)
      D1_2D,

      FIRST_EDGE_1D,
      LAST_EDGE_1D,
      CENTER_1D,
    };

    enum class TrianglePosition : char {
      //    e--------f
      //   /|       /|
      //  / |      / |
      // a--------b  |
      // |  g-----|--h
      // | /      | /
      // |/       |/
      // c--------d

      F_3D, // face (abc, bcd)
      C_3D, // side (abe, bef)
      H_3D, // top (acg, aeg)
      D1_3D, // diagonal1 (bdg, beg)
      D2_3D, // diagonal2 (abg, bgh)
      D3_3D, // diagonal3 (bcg, bfg)

      TOP_2D, // abc
      BOTTOM_2D, // bcd
    };

    bool hasPreconditionedVerticesAndCells_{false};

    // acceleration functions
    int checkAcceleration();
    bool isPowerOfTwo(unsigned long long int v, unsigned long long int &r);

    //\cond
    // 2D //
    void vertexToPosition2d(const SimplexId vertex,
                            SimplexId p[2]) const override;
    void
      edgeToPosition2d(const SimplexId edge, const int k, SimplexId p[2]) const;
    void triangleToPosition2d(const SimplexId triangle,
                              SimplexId p[2]) const override;

    SimplexId getVertexNeighbor2d(const SimplexId p[2],
                                  const SimplexId v,
                                  const int id) const;
    SimplexId getVertexEdge2d(const SimplexId p[2], const int id) const;
    SimplexId getVertexStar2d(const SimplexId p[2], const int id) const;
    SimplexId getVertexLink2d(const SimplexId p[2], const int id) const;

    SimplexId getEdgeTriangle2dL(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangle2dH(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangle2dD1(const SimplexId p[3], const int id) const;

    SimplexId getEdgeLink2dL(const SimplexId p[2], const int id) const;
    SimplexId getEdgeLink2dH(const SimplexId p[2], const int id) const;
    SimplexId getEdgeLink2dD1(const SimplexId p[2], const int id) const;

    SimplexId getEdgeStar2dL(const SimplexId p[2], const int id) const;
    SimplexId getEdgeStar2dH(const SimplexId p[2], const int id) const;

    // 3D //
    void vertexToPosition(const SimplexId vertex,
                          SimplexId p[3]) const override;
    void
      edgeToPosition(const SimplexId edge, const int k, SimplexId p[3]) const;
    void triangleToPosition(const SimplexId triangle,
                            const int k,
                            SimplexId p[3]) const override;
    void tetrahedronToPosition(const SimplexId tetrahedron,
                               SimplexId p[3]) const override;

    SimplexId getVertexNeighbor3d(const SimplexId p[3],
                                  const SimplexId v,
                                  const int id) const;
    SimplexId getVertexEdge3d(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangle3d(const SimplexId p[3], const int id) const;
    SimplexId getVertexLink3d(const SimplexId p[3], const int id) const;
    SimplexId getVertexStar3d(const SimplexId p[3], const int id) const;

    SimplexId getEdgeTriangle3dL(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangle3dH(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangle3dP(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangle3dD1(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangle3dD2(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangle3dD3(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangle3dD4(const SimplexId p[3], const int id) const;

    SimplexId getEdgeLinkL(const SimplexId p[3], const int id) const;
    SimplexId getEdgeLinkH(const SimplexId p[3], const int id) const;
    SimplexId getEdgeLinkP(const SimplexId p[3], const int id) const;
    SimplexId getEdgeLinkD1(const SimplexId p[3], const int id) const;
    SimplexId getEdgeLinkD2(const SimplexId p[3], const int id) const;
    SimplexId getEdgeLinkD3(const SimplexId p[3], const int id) const;
    SimplexId getEdgeLinkD4(const SimplexId p[3], const int id) const;

    SimplexId getEdgeStarL(const SimplexId p[3], const int id) const;
    SimplexId getEdgeStarH(const SimplexId p[3], const int id) const;
    SimplexId getEdgeStarP(const SimplexId p[3], const int id) const;
    SimplexId getEdgeStarD1(const SimplexId p[3], const int id) const;
    SimplexId getEdgeStarD2(const SimplexId p[3], const int id) const;
    SimplexId getEdgeStarD3(const SimplexId p[3], const int id) const;

    SimplexId getTriangleVertexF(const SimplexId p[3], const int id) const;
    SimplexId getTriangleVertexH(const SimplexId p[3], const int id) const;
    SimplexId getTriangleVertexC(const SimplexId p[3], const int id) const;
    SimplexId getTriangleVertexD1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleVertexD2(const SimplexId p[3], const int id) const;
    SimplexId getTriangleVertexD3(const SimplexId p[3], const int id) const;

    SimplexId getTriangleEdgeF_0(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeF_1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeH_0(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeH_1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeC_0(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeC_1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeD1_0(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeD1_1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeD2_0(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeD2_1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeD3_0(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeD3_1(const SimplexId p[3], const int id) const;

    SimplexId getTriangleLinkF(const SimplexId p[3], const int id) const;
    SimplexId getTriangleLinkH(const SimplexId p[3], const int id) const;
    SimplexId getTriangleLinkC(const SimplexId p[3], const int id) const;
    SimplexId getTriangleLinkD1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleLinkD2(const SimplexId p[3], const int id) const;
    SimplexId getTriangleLinkD3(const SimplexId p[3], const int id) const;

    SimplexId getTriangleStarF(const SimplexId p[3], const int id) const;
    SimplexId getTriangleStarH(const SimplexId p[3], const int id) const;
    SimplexId getTriangleStarC(const SimplexId p[3], const int id) const;
    SimplexId getTriangleStarD1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleStarD2(const SimplexId p[3], const int id) const;
    SimplexId getTriangleStarD3(const SimplexId p[3], const int id) const;

    SimplexId getTetrahedronVertexABCG(const SimplexId p[3],
                                       const int id) const;
    SimplexId getTetrahedronVertexBCDG(const SimplexId p[3],
                                       const int id) const;
    SimplexId getTetrahedronVertexABEG(const SimplexId p[3],
                                       const int id) const;
    SimplexId getTetrahedronVertexBEFG(const SimplexId p[3],
                                       const int id) const;
    SimplexId getTetrahedronVertexBFGH(const SimplexId p[3],
                                       const int id) const;
    SimplexId getTetrahedronVertexBDGH(const SimplexId p[3],
                                       const int id) const;

    SimplexId getTetrahedronEdgeABCG(const SimplexId p[3], const int id) const;
    SimplexId getTetrahedronEdgeBCDG(const SimplexId p[3], const int id) const;
    SimplexId getTetrahedronEdgeABEG(const SimplexId p[3], const int id) const;
    SimplexId getTetrahedronEdgeBEFG(const SimplexId p[3], const int id) const;
    SimplexId getTetrahedronEdgeBFGH(const SimplexId p[3], const int id) const;
    SimplexId getTetrahedronEdgeBDGH(const SimplexId p[3], const int id) const;

    SimplexId getTetrahedronTriangleABCG(const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronTriangleBCDG(const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronTriangleABEG(const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronTriangleBEFG(const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronTriangleBFGH(const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronTriangleBDGH(const SimplexId p[3],
                                         const int id) const;

    SimplexId getTetrahedronNeighborABCG(const SimplexId t,
                                         const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronNeighborBCDG(const SimplexId t,
                                         const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronNeighborABEG(const SimplexId t,
                                         const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronNeighborBEFG(const SimplexId t,
                                         const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronNeighborBFGH(const SimplexId t,
                                         const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronNeighborBDGH(const SimplexId t,
                                         const SimplexId p[3],
                                         const int id) const;
    //\endcond

#ifdef TTK_ENABLE_MPI

  protected:
    int preconditionDistributedCells() override;
    // std::shared_ptr<PeriodicImplicitTriangulation> metaGrid_;
    std::array<unsigned char, 6> isBoundaryPeriodic{};

  public:
    void createMetaGrid(const double *const bounds) override;
    void setIsBoundaryPeriodic(std::array<unsigned char, 6> boundary);
    int getCellRankInternal(const SimplexId lcid) const override;

  protected:
    std::array<SimplexId, 3>
      getVertGlobalCoords(const SimplexId lvid) const override;
    std::array<SimplexId, 3>
      getVertLocalCoords(const SimplexId gvid) const override;

#endif // TTK_ENABLE_MPI
  };

  template <size_t card, typename Derived>
  class PeriodicImplicitTriangulationCRTP
    : public PeriodicImplicitTriangulation<card> {
    inline Derived &underlying() {
      return static_cast<Derived &>(*this);
    }
    inline Derived const &underlying() const {
      return static_cast<Derived const &>(*this);
    }

  public:
    int TTK_TRIANGULATION_INTERNAL(getVertexNeighbor)(
      const SimplexId &vertexId,
      const int &localNeighborId,
      SimplexId &neighborId) const override;

    int getVertexEdgeInternal(const SimplexId &vertexId,
                              const int &id,
                              SimplexId &edgeId) const override;

    int getVertexTriangleInternal(const SimplexId &vertexId,
                                  const int &id,
                                  SimplexId &triangleId) const override;

    int TTK_TRIANGULATION_INTERNAL(getVertexLink)(
      const SimplexId &vertexId,
      const int &localLinkId,
      SimplexId &linkId) const override;

    int TTK_TRIANGULATION_INTERNAL(getVertexStar)(
      const SimplexId &vertexId,
      const int &localStarId,
      SimplexId &starId) const override;

    int TTK_TRIANGULATION_INTERNAL(getVertexPoint)(const SimplexId &vertexId,
                                                   float &x,
                                                   float &y,
                                                   float &z) const override;

    int getEdgeVertexInternal(const SimplexId &edgeId,
                              const int &localVertexId,
                              SimplexId &vertexId) const override;

    SimplexId
      getEdgeTriangleNumberInternal(const SimplexId &edgeId) const override;

    int getEdgeTriangleInternal(const SimplexId &edgeId,
                                const int &id,
                                SimplexId &triangleId) const override;

    int
      TTK_TRIANGULATION_INTERNAL(getEdgeLink)(const SimplexId &edgeId,
                                              const int &localLinkId,
                                              SimplexId &linkId) const override;

    SimplexId TTK_TRIANGULATION_INTERNAL(getEdgeStarNumber)(
      const SimplexId &edgeId) const override;

    int
      TTK_TRIANGULATION_INTERNAL(getEdgeStar)(const SimplexId &edgeId,
                                              const int &localStarId,
                                              SimplexId &starId) const override;

    int getTriangleVertexInternal(const SimplexId &triangleId,
                                  const int &localVertexId,
                                  SimplexId &vertexId) const override;

    int getTriangleEdgeInternal(const SimplexId &triangleId,
                                const int &id,
                                SimplexId &edgeId) const override;

    int TTK_TRIANGULATION_INTERNAL(getTriangleLink)(
      const SimplexId &triangleId,
      const int &localLinkId,
      SimplexId &linkId) const override;

    int TTK_TRIANGULATION_INTERNAL(getTriangleStar)(
      const SimplexId &triangleId,
      const int &localStarId,
      SimplexId &starId) const override;

    SimplexId TTK_TRIANGULATION_INTERNAL(getTriangleStarNumber)(
      const SimplexId &triangleId) const override;

    int getTriangleNeighbor(const SimplexId &triangleId,
                            const int &localNeighborId,
                            SimplexId &neighborId) const override;

    int getTetrahedronVertex(const SimplexId &tetId,
                             const int &localVertexId,
                             SimplexId &vertexId) const override;

    int getTetrahedronEdge(const SimplexId &tetId,
                           const int &id,
                           SimplexId &edgeId) const override;

    int getTetrahedronTriangle(const SimplexId &tetId,
                               const int &id,
                               SimplexId &triangleId) const override;

    int getTetrahedronNeighbor(const SimplexId &tetId,
                               const int &localNeighborId,
                               SimplexId &neighborId) const override;

    /**
     * Compute the barycenter of the points of the given edge identifier.
     */
    virtual int getEdgeIncenter(SimplexId edgeId, float incenter[3]) const {
      SimplexId v0{}, v1{};
      getEdgeVertexInternal(edgeId, 0, v0);
      getEdgeVertexInternal(edgeId, 1, v1);

      std::array<float, 3> p0{}, p1{};
      this->getVertexPointInternal(v0, p0[0], p0[1], p0[2]);
      this->getVertexPointInternal(v1, p1[0], p1[1], p1[2]);

      const auto &ind0 = this->underlying().getVertexCoords(v0);
      const auto &ind1 = this->underlying().getVertexCoords(v1);

      for(int i = 0; i < this->dimensionality_; ++i) {
        if(ind1[i] == this->nbvoxels_[i]) {
          p0[i] += (ind0[i] == 0) * this->dimensions_[i] * this->spacing_[i];
        } else if(ind0[i] == this->nbvoxels_[i]) {
          p1[i] += (ind1[i] == 0) * this->dimensions_[i] * this->spacing_[i];
        }
      }

      for(int i = 0; i < 3; ++i) {
        incenter[i] = 0.5f * (p0[i] + p1[i]);
      }

      return 0;
    }

    /**
     * Compute the incenter of the points of the given triangle
     * identifier.
     */
    virtual int getTriangleIncenter(SimplexId triangleId,
                                    float incenter[3]) const {

      SimplexId v0{}, v1{}, v2{};
      getTriangleVertexInternal(triangleId, 0, v0);
      getTriangleVertexInternal(triangleId, 1, v1);
      getTriangleVertexInternal(triangleId, 2, v2);

      std::array<float, 3> p0{}, p1{}, p2{};
      this->getVertexPointInternal(v0, p0[0], p0[1], p0[2]);
      this->getVertexPointInternal(v1, p1[0], p1[1], p1[2]);
      this->getVertexPointInternal(v2, p2[0], p2[1], p2[2]);

      const auto &ind0 = this->underlying().getVertexCoords(v0);
      const auto &ind1 = this->underlying().getVertexCoords(v1);
      const auto &ind2 = this->underlying().getVertexCoords(v2);

      for(int i = 0; i < this->dimensionality_; ++i) {
        if(ind0[i] == this->nbvoxels_[i]) {
          p1[i] += (ind1[i] == 0) * this->dimensions_[i] * this->spacing_[i];
          p2[i] += (ind2[i] == 0) * this->dimensions_[i] * this->spacing_[i];
        } else if(ind1[i] == this->nbvoxels_[i]) {
          p0[i] += (ind0[i] == 0) * this->dimensions_[i] * this->spacing_[i];
          p2[i] += (ind2[i] == 0) * this->dimensions_[i] * this->spacing_[i];
        } else if(ind2[i] == this->nbvoxels_[i]) {
          p0[i] += (ind0[i] == 0) * this->dimensions_[i] * this->spacing_[i];
          p1[i] += (ind1[i] == 0) * this->dimensions_[i] * this->spacing_[i];
        }
      }

      std::array<float, 3> d{Geometry::distance(p1.data(), p2.data()),
                             Geometry::distance(p2.data(), p0.data()),
                             Geometry::distance(p0.data(), p1.data())};
      const float sum = d[0] + d[1] + d[2];
      for(int i = 0; i < 3; ++i) {
        incenter[i] = (d[0] * p0[i] + d[1] * p1[i] + d[2] * p2[i]) / sum;
      }

      return 0;
    }

    /**
     * Compute the barycenter of the incenters of the triangles of the
     * given tetra identifier.
     */
    virtual int getTetraIncenter(SimplexId tetraId, float incenter[3]) const {

      SimplexId v0{}, v1{}, v2{}, v3{};
      this->getCellVertexInternal(tetraId, 0, v0);
      this->getCellVertexInternal(tetraId, 1, v1);
      this->getCellVertexInternal(tetraId, 2, v2);
      this->getCellVertexInternal(tetraId, 3, v3);

      std::array<float, 3> p0{}, p1{}, p2{}, p3{};
      this->getVertexPointInternal(v0, p0[0], p0[1], p0[2]);
      this->getVertexPointInternal(v1, p1[0], p1[1], p1[2]);
      this->getVertexPointInternal(v2, p2[0], p2[1], p2[2]);
      this->getVertexPointInternal(v3, p3[0], p3[1], p3[2]);

      const auto &ind0 = this->underlying().getVertexCoords(v0);
      const auto &ind1 = this->underlying().getVertexCoords(v1);
      const auto &ind2 = this->underlying().getVertexCoords(v2);
      const auto &ind3 = this->underlying().getVertexCoords(v3);

      for(int i = 0; i < this->dimensionality_; ++i) {
        if(ind0[i] == this->nbvoxels_[i]) {
          p1[i] += (ind1[i] == 0) * this->dimensions_[i] * this->spacing_[i];
          p2[i] += (ind2[i] == 0) * this->dimensions_[i] * this->spacing_[i];
          p3[i] += (ind3[i] == 0) * this->dimensions_[i] * this->spacing_[i];
        } else if(ind1[i] == this->nbvoxels_[i]) {
          p0[i] += (ind0[i] == 0) * this->dimensions_[i] * this->spacing_[i];
          p2[i] += (ind2[i] == 0) * this->dimensions_[i] * this->spacing_[i];
          p3[i] += (ind3[i] == 0) * this->dimensions_[i] * this->spacing_[i];
        } else if(ind2[i] == this->nbvoxels_[i]) {
          p0[i] += (ind0[i] == 0) * this->dimensions_[i] * this->spacing_[i];
          p1[i] += (ind1[i] == 0) * this->dimensions_[i] * this->spacing_[i];
          p3[i] += (ind3[i] == 0) * this->dimensions_[i] * this->spacing_[i];
        } else if(ind3[i] == this->nbvoxels_[i]) {
          p0[i] += (ind0[i] == 0) * this->dimensions_[i] * this->spacing_[i];
          p1[i] += (ind1[i] == 0) * this->dimensions_[i] * this->spacing_[i];
          p2[i] += (ind2[i] == 0) * this->dimensions_[i] * this->spacing_[i];
        }
      }

      for(int i = 0; i < 3; ++i) {
        incenter[i] = 0.25f * (p0[i] + p1[i] + p2[i] + p3[i]);
      }
      return 0;
    }
  };

} // namespace ttk

/// @cond
template <size_t card>
inline void
  ttk::PeriodicImplicitTriangulation<card>::vertexToPosition2d(const SimplexId vertex,
                                                         SimplexId p[2]) const {
  if(isAccelerated_) {
    p[0] = vertex & mod_[0];
    p[1] = vertex >> div_[0];
  } else {
    p[0] = vertex % vshift_[0];
    p[1] = vertex / vshift_[0];
  }
}

template <size_t card>
inline void ttk::PeriodicImplicitTriangulation<card>::edgeToPosition2d(
  const SimplexId edge, const int k, SimplexId p[2]) const {
  const int e = (k) ? edge - esetshift_[k - 1] : edge;
  p[0] = e % eshift_[2 * k];
  p[1] = e / eshift_[2 * k];
}

template <size_t card>
inline void ttk::PeriodicImplicitTriangulation<card>::triangleToPosition2d(
  const SimplexId triangle, SimplexId p[2]) const {
  p[0] = triangle % tshift_[0];
  p[1] = triangle / tshift_[0];
}

template <size_t card>
inline ttk::SimplexId ttk::PeriodicImplicitTriangulation<card>::getVertexNeighbor2d(
  const SimplexId p[2], const SimplexId v, const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapYBottom = 0;
  if(p[0] == 0)
    wrapXLeft = wrap_[0];
  if(p[0] == nbvoxels_[Di_])
    wrapXRight = -wrap_[0];
  if(p[1] == 0)
    wrapYTop = wrap_[1];
  if(p[1] == nbvoxels_[Dj_])
    wrapYBottom = -wrap_[1];
  switch(id) {
    case 0:
      return v - 1 + wrapXLeft;
    case 1:
      return v - vshift_[0] + wrapYTop;
    case 2:
      return v - vshift_[0] + 1 + wrapXRight + wrapYTop;
    case 3:
      return v + 1 + wrapXRight;
    case 4:
      return v + vshift_[0] + wrapYBottom;
    case 5:
      return v + vshift_[0] - 1 + wrapXLeft + wrapYBottom;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getVertexEdge2d(const SimplexId p[2],
                                                      const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapYTop = 0;
  if(p[0] == 0)
    wrapXLeft = wrap_[0];
  if(p[1] == 0)
    wrapYTop = wrap_[1];
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2] + wrapYTop;
    case 1:
      return p[0] + p[1] * eshift_[0] - 1 + wrapXLeft;
    case 2:
      return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4] + wrapYTop;
    case 3:
      return p[0] + p[1] * eshift_[0];
    case 4:
      return esetshift_[0] + p[0] + p[1] * eshift_[2];
    case 5:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] - 1 + wrapXLeft;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getVertexStar2d(const SimplexId p[2],
                                                      const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapYTop = 0;
  if(p[0] == 0)
    wrapXLeft = 2 * wrap_[0];
  if(p[1] == 0)
    wrapYTop = 2 * wrap_[1];
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + wrapXLeft;
    case 1:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1 + wrapXLeft;
    case 2:
      return p[0] * 2 + p[1] * tshift_[0];
    case 3:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + wrapYTop;
    case 4:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1 + wrapYTop;
    case 5:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + 1 + wrapXLeft
             + wrapYTop;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getVertexLink2d(const SimplexId p[2],
                                                      const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapYBottom = 0;
  if(p[0] == 0)
    wrapXLeft = wrap_[0];
  if(p[0] == nbvoxels_[Di_])
    wrapXRight = -wrap_[0];
  if(p[1] == 0)
    wrapYTop = wrap_[1];
  if(p[1] == nbvoxels_[Dj_])
    wrapYBottom = -wrap_[1];
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] - 1 + wrapXLeft;
    case 1:
      return p[0] + (p[1] + 1) * eshift_[0] - 1 + wrapXLeft + wrapYBottom;
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4];
    case 3:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2] + 1 + wrapXRight
             + wrapYTop;
    case 4:
      return p[0] + (p[1] - 1) * eshift_[0] + wrapYTop;
    case 5:
      return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4] - 1 + wrapXLeft
             + wrapYTop;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeTriangle2dL(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapYTop = 0;
  if(p[1] == 0)
    wrapYTop = 2 * wrap_[1];
  switch(id) {
    case 0:
      return p[Di_] * 2 + p[Dj_] * tshift_[0];
    case 1:
      return p[Di_] * 2 + (p[Dj_] - 1) * tshift_[0] + 1 + wrapYTop;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeTriangle2dH(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapXLeft = 0;
  if(p[0] == 0)
    wrapXLeft = 2 * wrap_[0];
  switch(id) {
    case 0:
      return p[Di_] * 2 + p[Dj_] * tshift_[0];
    case 1:
      return (p[Di_] - 1) * 2 + p[Dj_] * tshift_[0] + 1 + wrapXLeft;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeTriangle2dD1(const SimplexId p[3],
                                                          const int id) const {
  switch(id) {
    case 0:
      return p[Di_] * 2 + p[Dj_] * tshift_[0];
    case 1:
      return p[Di_] * 2 + p[Dj_] * tshift_[0] + 1;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeLink2dL(const SimplexId p[2],
                                                     const int id) const {
  if(p[1] > 0 and p[1] < nbvoxels_[Dj_]) {
    switch(id) {
      case 0:
        return p[0] + (p[1] + 1) * vshift_[0];
      case 1:
        return ((p[0] < nbvoxels_[Di_])
                  ? (p[0] + (p[1] - 1) * vshift_[0] + 1)
                  : (p[0] + (p[1] - 1) * vshift_[0] + 1 - wrap_[0]));
    }
    return -1;
  } else if(p[1] == 0) {
    switch(id) {
      case 0:
        return p[0] + (p[1] + 1) * vshift_[0];
      case 1:
        return ((p[0] < nbvoxels_[Di_])
                  ? (p[0] + (p[1] - 1) * vshift_[0] + 1)
                  : (p[0] + (p[1] - 1) * vshift_[0] + 1 - wrap_[0]))
               + wrap_[1];
    }
    return -1;
  } else {
    switch(id) {
      case 0:
        return p[0] + (p[1] + 1) * vshift_[0] - wrap_[1];
      case 1:
        return ((p[0] < nbvoxels_[Di_])
                  ? (p[0] + (p[1] - 1) * vshift_[0] + 1)
                  : (p[0] + (p[1] - 1) * vshift_[0] + 1 - wrap_[0]));
    }
    return -1;
  }
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeLink2dH(const SimplexId p[2],
                                                     const int id) const {
  if(p[0] > 0 and p[0] < nbvoxels_[Di_]) {
    switch(id) {
      case 0:
        return p[0] + p[1] * vshift_[0] + 1;
      case 1:
        return ((p[1] < nbvoxels_[Dj_])
                  ? (p[0] + (p[1] + 1) * vshift_[0] - 1)
                  : (p[0] + (p[1] + 1) * vshift_[0] - 1 - wrap_[1]));
    }
    return -1;
  } else if(p[0] == 0) {
    switch(id) {
      case 0:
        return p[0] + p[1] * vshift_[0] + 1 + wrap_[0];
      case 1:
        return ((p[1] < nbvoxels_[Dj_])
                  ? (p[0] + (p[1] + 1) * vshift_[0] - 1)
                  : (p[0] + (p[1] + 1) * vshift_[0] - 1 - wrap_[1]));
    }
    return -1;
  } else {
    switch(id) {
      case 0:
        return p[0] + p[1] * vshift_[0] + 1;
      case 1:
        return ((p[1] < nbvoxels_[Dj_])
                  ? (p[0] + (p[1] + 1) * vshift_[0] - 1)
                  : (p[0] + (p[1] + 1) * vshift_[0] - 1 - wrap_[1]))
               - wrap_[0];
    }
    return -1;
  }
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeLink2dD1(const SimplexId p[2],
                                                      const int id) const {
  const SimplexId wrapX = (p[0] < nbvoxels_[Di_]) ? 0 : wrap_[0];
  const SimplexId wrapY = (p[1] < nbvoxels_[Dj_]) ? 0 : wrap_[1];
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0];
    case 1:
      return p[0] + (p[1] + 1) * vshift_[0] + 1 - wrapX - wrapY;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeStar2dL(const SimplexId p[2],
                                                     const int id) const {
  if(p[1] == 0) {
    switch(id) {
      case 0:
        return p[0] * 2 + p[1] * tshift_[0];
      case 1:
        return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1 + 2 * wrap_[1];
    }
    return -1;
  } else {
    switch(id) {
      case 0:
        return p[0] * 2 + p[1] * tshift_[0];
      case 1:
        return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
    }
    return -1;
  }
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeStar2dH(const SimplexId p[2],
                                                     const int id) const {
  if(p[0] == 0) {
    switch(id) {
      case 0:
        return p[0] * 2 + p[1] * tshift_[0];
      case 1:
        return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1 + 2 * wrap_[0];
    }
    return -1;
  } else {
    switch(id) {
      case 0:
        return p[0] * 2 + p[1] * tshift_[0];
      case 1:
        return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
    }
    return -1;
  }
}

template <size_t card>
inline void
  ttk::PeriodicImplicitTriangulation<card>::vertexToPosition(const SimplexId vertex,
                                                       SimplexId p[3]) const {
  if(isAccelerated_) {
    p[0] = vertex & mod_[0];
    p[1] = (vertex & mod_[1]) >> div_[0];
    p[2] = vertex >> div_[1];
  } else {
    p[0] = vertex % vshift_[0];
    p[1] = (vertex % vshift_[1]) / vshift_[0];
    p[2] = vertex / vshift_[1];
  }
}

template <size_t card>
inline void ttk::PeriodicImplicitTriangulation<card>::edgeToPosition(
  const SimplexId edge, const int k, SimplexId p[3]) const {
  const int e = (k) ? edge - esetshift_[k - 1] : edge;
  p[0] = e % eshift_[2 * k];
  p[1] = (e % eshift_[2 * k + 1]) / eshift_[2 * k];
  p[2] = e / eshift_[2 * k + 1];
}

template <size_t card>
inline void ttk::PeriodicImplicitTriangulation<card>::triangleToPosition(
  const SimplexId triangle, const int k, SimplexId p[3]) const {
  const SimplexId t = (k) ? triangle - tsetshift_[k - 1] : triangle;
  p[0] = t % tshift_[2 * k];
  p[1] = (t % tshift_[2 * k + 1]) / tshift_[2 * k];
  p[2] = t / tshift_[2 * k + 1];
}

template <size_t card>
inline void ttk::PeriodicImplicitTriangulation<card>::tetrahedronToPosition(
  const SimplexId tetrahedron, SimplexId p[3]) const {
  p[0] = (tetrahedron % tetshift_[0]) / 6;
  p[1] = (tetrahedron % tetshift_[1]) / tetshift_[0];
  p[2] = tetrahedron / tetshift_[1];
}

template <size_t card>
inline ttk::SimplexId ttk::PeriodicImplicitTriangulation<card>::getVertexNeighbor3d(
  const SimplexId p[3], const SimplexId v, const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZBack = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == 0)
    wrapXLeft = wrap_[0];
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == 0)
    wrapYTop = wrap_[1];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == 0)
    wrapZBack = wrap_[2];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return v - vshift_[0] - vshift_[1] + wrapYTop + wrapZBack;
    case 1:
      return v + 1 - vshift_[0] - vshift_[1] + wrapXRight + wrapYTop
             + wrapZBack;
    case 2:
      return v - vshift_[1] + wrapZBack;
    case 3:
      return v + 1 - vshift_[1] + wrapXRight + wrapZBack;
    case 4:
      return v - vshift_[0] + wrapYTop;
    case 5:
      return v + 1 - vshift_[0] + wrapXRight + wrapYTop;
    case 6:
      return v + 1 + wrapXRight;
    case 7:
      return v - 1 + vshift_[1] + wrapXLeft + wrapZFront;
    case 8:
      return v + vshift_[1] + wrapZFront;
    case 9:
      return v - 1 + wrapXLeft;
    case 10:
      return v - 1 + vshift_[0] + wrapXLeft + wrapYBottom;
    case 11:
      return v + vshift_[0] + wrapYBottom;
    case 12:
      return v - 1 + vshift_[0] + vshift_[1] + wrapXLeft + wrapYBottom
             + wrapZFront;
    case 13:
      return v + vshift_[0] + vshift_[1] + wrapYBottom + wrapZFront;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getVertexEdge3d(const SimplexId p[3],
                                                      const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapZBack = 0;
  if(p[0] == 0)
    wrapXLeft = wrap_[0];
  if(p[1] == 0)
    wrapYTop = wrap_[1];
  if(p[2] == 0)
    wrapZBack = wrap_[2];
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9] + wrapYTop + wrapZBack;
    case 1:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + (p[2] - 1) * eshift_[13] + wrapYTop + wrapZBack;
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] + (p[2] - 1) * eshift_[5]
             + wrapZBack;
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11] + wrapZBack;
    case 4:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2] + p[2] * eshift_[3]
             + wrapYTop;
    case 5:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6] + p[2] * eshift_[7]
             + wrapYTop;
    case 6:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1];
    case 7:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11] - 1
             + wrapXLeft;
    case 8:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] + p[2] * eshift_[5];
    case 9:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1 + wrapXLeft;
    case 10:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7] - 1
             + wrapXLeft;
    case 11:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] + p[2] * eshift_[3];
    case 12:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13] - 1
             + wrapXLeft;
    case 13:
      return esetshift_[3] + p[0] + p[1] * eshift_[8] + p[2] * eshift_[9];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getVertexLink3d(const SimplexId p[3],
                                                      const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZBack = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == 0)
    wrapXLeft = 2 * wrap_[0];
  if(p[0] == nbvoxels_[0])
    wrapXRight = -2 * wrap_[0];
  if(p[1] == 0)
    wrapYTop = 2 * wrap_[1];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -2 * wrap_[1];
  if(p[2] == 0)
    wrapZBack = 2 * wrap_[2];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -2 * wrap_[2];
  switch(id) {
    case 0:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11];
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1;
    case 2:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + wrapXLeft;
    case 3:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + 1 + wrapXLeft;
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + wrapXLeft + wrapYBottom;
    case 5:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1 + wrapXLeft + wrapYBottom;
    case 6:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
             + wrapXLeft + wrapZFront;
    case 7:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1] + 1
             + wrapXLeft + wrapZFront;
    case 8:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + wrapYTop;
    case 9:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7] + wrapYTop;
    case 10:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11] + wrapXLeft + wrapYTop;
    case 11:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + 1 + wrapXLeft + wrapYTop;
    case 12:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + wrapZBack;
    case 13:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1 + wrapZBack;
    case 14:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1 + wrapXLeft + wrapZBack;
    case 15:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1 + wrapXLeft + wrapZBack;
    case 16:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + wrapXRight + wrapYTop + wrapZBack;
    case 17:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1 + wrapXRight + wrapYTop + wrapZBack;
    case 18:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + wrapYTop + wrapZBack;
    case 19:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1 + wrapYTop + wrapZBack;
    case 20:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1]
             + wrapYTop + wrapZBack;
    case 21:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1] + 1
             + wrapYTop + wrapZBack;
    case 22:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1 + wrapXLeft + wrapYTop + wrapZBack;
    case 23:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7] + wrapXLeft + wrapYTop + wrapZBack;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getVertexTriangle3d(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapZBack = 0;
  if(p[0] == 0)
    wrapXLeft = 2 * wrap_[0];
  if(p[1] == 0)
    wrapYTop = 2 * wrap_[1];
  if(p[2] == 0)
    wrapZBack = 2 * wrap_[2];
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1
             + wrapXLeft;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11] + wrapXLeft;
    case 2:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + p[2] * tshift_[7] + wrapXLeft;
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + p[2] * tshift_[9] + 1 + wrapXLeft;
    case 4:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 5:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11] + 1 + wrapXLeft;
    case 7:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1 + wrapXLeft;
    case 8:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + p[2] * tshift_[7] + 1 + wrapXLeft;
    case 9:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + p[2] * tshift_[9] + wrapXLeft;
    case 10:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + wrapXLeft;
    case 11:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + wrapXLeft;
    case 12:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7] + wrapYTop + wrapZBack;
    case 13:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1 + wrapYTop + wrapZBack;
    case 14:
      return p[0] * 2 + tsetshift_[1] + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + wrapYTop + wrapZBack;
    case 15:
      return p[0] * 2 + tsetshift_[1] + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1 + wrapYTop + wrapZBack;
    case 16:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + wrapYTop + wrapZBack;
    case 17:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1 + wrapYTop + wrapZBack;
    case 18:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + wrapYTop + wrapZBack;
    case 19:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1 + wrapYTop + wrapZBack;
    case 20:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + wrapYTop;
    case 21:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1
             + wrapYTop;
    case 22:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + wrapZBack;
    case 23:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1 + wrapZBack;
    case 24:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1
             + wrapXLeft + wrapYTop;
    case 25:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7] + wrapXLeft + wrapYTop;
    case 26:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5] + wrapYTop;
    case 27:
      return p[0] * 2 + tsetshift_[2] + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1 + wrapZBack;
    case 28:
      return p[0] * 2 + tsetshift_[1] + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1 + wrapZBack;
    case 29:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 30:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1 + wrapXLeft + wrapYTop + wrapZBack;
    case 31:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1 + wrapXLeft + wrapZBack;
    case 32:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 33:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
    case 34:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1 + wrapXLeft + wrapZBack;
    case 35:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11] + wrapYTop;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getVertexStar3d(const SimplexId p[3],
                                                      const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapZBack = 0;
  if(p[0] == 0)
    wrapXLeft = 6 * wrap_[0];
  if(p[1] == 0)
    wrapYTop = 6 * wrap_[1];
  if(p[2] == 0)
    wrapZBack = 6 * wrap_[2];
  switch(id) {
    case 0:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + wrapXLeft;
    case 1:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 1
             + wrapXLeft;
    case 2:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 2
             + wrapXLeft;
    case 3:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 3
             + wrapXLeft;
    case 4:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 4
             + wrapXLeft;
    case 5:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 5
             + wrapXLeft;
    case 6:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + wrapYTop + wrapZBack;
    case 7:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 1 + wrapYTop + wrapZBack;
    case 8:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 2 + wrapYTop + wrapZBack;
    case 9:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3 + wrapYTop + wrapZBack;
    case 10:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4 + wrapYTop + wrapZBack;
    case 11:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 5 + wrapYTop + wrapZBack;
    case 12:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1];
    case 13:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 2;
    case 14:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + wrapYTop;
    case 15:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1] + 1
             + wrapYTop;
    case 16:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1 + wrapXLeft + wrapYTop;
    case 17:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5 + wrapXLeft + wrapYTop;
    case 18:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1] + 2
             + wrapZBack;
    case 19:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1] + 3
             + wrapZBack;
    case 20:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3 + wrapXLeft + wrapZBack;
    case 21:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4 + wrapXLeft + wrapZBack;
    case 22:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 4 + wrapXLeft + wrapYTop + wrapZBack;
    case 23:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 5 + wrapXLeft + wrapYTop + wrapZBack;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeTriangle3dL(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapYTop = 0;
  SimplexId wrapZBack = 0;
  if(p[1] == 0)
    wrapYTop = 2 * wrap_[1];
  if(p[2] == 0)
    wrapZBack = 2 * wrap_[2];
  switch(id) {
    case 0:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1
             + wrapYTop;
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 2:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1 + wrapZBack;
    case 3:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1 + wrapYTop + wrapZBack;
    case 4:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 5:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeTriangle3dH(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapZBack = 0;
  if(p[0] == 0)
    wrapXLeft = 2 * wrap_[0];
  if(p[2] == 0)
    wrapZBack = 2 * wrap_[2];
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1
             + wrapXLeft;
    case 1:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + p[2] * tshift_[7] + wrapXLeft;
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 3:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1 + wrapZBack;
    case 4:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1 + wrapZBack;
    case 5:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeTriangle3dP(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapYTop = 0;
  if(p[0] == 0)
    wrapXLeft = 2 * wrap_[0];
  if(p[1] == 0)
    wrapYTop = 2 * wrap_[1];
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1 + wrapXLeft;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11] + 1 + wrapXLeft;
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 3:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 4:
      return tsetshift_[4] + p[0] * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11] + wrapYTop;
    case 5:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5] + wrapYTop;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeTriangle3dD1(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapZBack = 0;
  if(p[2] == 0)
    wrapZBack = 2 * wrap_[2];
  switch(id) {
    case 0:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 1:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 2:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11];
    case 3:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1 + wrapZBack;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeTriangle3dD2(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapXLeft = 0;
  if(p[0] == 0)
    wrapXLeft = 2 * wrap_[0];
  switch(id) {
    case 0:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 1:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + p[2] * tshift_[9] + 1 + wrapXLeft;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeTriangle3dD3(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapYTop = 0;
  if(p[1] == 0)
    wrapYTop = 2 * wrap_[1];
  switch(id) {
    case 0:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3]
             + 1;
    case 2:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1;
    case 3:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7] + wrapYTop;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeTriangle3dD4(const SimplexId p[3],
                                                          const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7];
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1;
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
    case 3:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9]
             + 1;
    case 4:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11];
    case 5:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11]
             + 1;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeLinkL(const SimplexId p[3],
                                                   const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZBack = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == 0)
    wrapYTop = wrap_[1];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == 0)
    wrapZBack = wrap_[2];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] + (p[1] + 1) * eshift_[4] + p[2] * eshift_[5]
             + wrapYBottom;
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] + (p[2] + 1) * eshift_[3]
             + wrapZFront;
    case 2:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + p[2] * eshift_[13] + wrapYTop;
    case 3:
      return esetshift_[1] + p[0] + 1 + (p[1] - 1) * eshift_[4]
             + (p[2] - 1) * eshift_[5] + wrapXRight + wrapYTop + wrapZBack;
    case 4:
      return esetshift_[0] + p[0] + 1 + (p[1] - 1) * eshift_[2]
             + (p[2] - 1) * eshift_[3] + wrapXRight + wrapYTop + wrapZBack;
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12]
             + (p[2] - 1) * eshift_[13] + wrapZBack;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeLinkH(const SimplexId p[3],
                                                   const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZBack = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == 0)
    wrapXLeft = wrap_[0];
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == 0)
    wrapZBack = wrap_[2];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + (p[2] - 1) * eshift_[1] + wrapZBack;
    case 1:
      return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1] - 1
             + wrapXLeft + wrapYBottom + wrapZFront;
    case 2:
      return esetshift_[1] + p[0] - 1 + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5] + wrapXLeft + wrapYBottom;
    case 3:
      return esetshift_[1] + p[0] + 1 + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5] + wrapXRight + wrapZBack;
    case 4:
      return esetshift_[5] + (p[0] - 1) + p[1] * eshift_[12]
             + (p[2] - 1) * eshift_[13] + wrapXLeft + wrapZBack;
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeLinkP(const SimplexId p[3],
                                                   const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == 0)
    wrapXLeft = wrap_[0];
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == 0)
    wrapYTop = wrap_[1];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + (p[1] - 1) * eshift_[0] + p[2] * eshift_[1] + wrapYTop;
    case 1:
      return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1] - 1
             + wrapXLeft + wrapYBottom + wrapZFront;
    case 2:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
    case 3:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + p[2] * eshift_[13] - 1 + wrapXLeft + wrapYTop;
    case 4:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] + (p[2] + 1) * eshift_[3]
             - 1 + wrapXLeft + wrapZFront;
    case 5:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2] + p[2] * eshift_[3]
             + 1 + wrapXRight + wrapYTop;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeLinkD1(const SimplexId p[3],
                                                    const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZBack = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == 0)
    wrapZBack = wrap_[2];
  switch(id) {
    case 0:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11] + wrapZBack;
    case 1:
      return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
             + p[2] * eshift_[11] + wrapYBottom;
    case 2:
      return esetshift_[3] + p[0] + p[1] * eshift_[8] + p[2] * eshift_[9];
    case 3:
      return esetshift_[3] + p[0] + p[1] * eshift_[8] + (p[2] - 1) * eshift_[9]
             + 1 + wrapXRight + wrapZBack;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeLinkD2(const SimplexId p[3],
                                                    const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == 0)
    wrapXLeft = wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11];
    case 1:
      return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
             + p[2] * eshift_[11] - 1 + wrapXLeft + wrapYBottom;
    case 2:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7];
    case 3:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + (p[2] + 1) * eshift_[7]
             - 1 + wrapXLeft + wrapZFront;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeLinkD3(const SimplexId p[3],
                                                    const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == 0)
    wrapYTop = wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6] + p[2] * eshift_[7]
             + wrapYTop;
    case 1:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + (p[2] + 1) * eshift_[7]
             + wrapZFront;
    case 2:
      return esetshift_[3] + p[0] + p[1] * eshift_[8] + p[2] * eshift_[9];
    case 3:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8] + p[2] * eshift_[9]
             + 1 + wrapXRight + wrapYTop;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeLinkD4(const SimplexId p[3],
                                                    const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + (p[1] + 1) * eshift_[0] + p[2] * eshift_[1] + wrapYBottom;
    case 1:
      return p[0] + p[1] * eshift_[0] + (p[2] + 1) * eshift_[1] + wrapZFront;
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] + p[2] * eshift_[5];
    case 3:
      return esetshift_[1] + p[0] + 1 + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5] + wrapXRight + wrapYBottom;
    case 4:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] + p[2] * eshift_[3];
    case 5:
      return esetshift_[0] + p[0] + 1 + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3] + wrapXRight + wrapZFront;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeStarL(const SimplexId p[3],
                                                   const int id) const {
  SimplexId wrapYTop = 0;
  SimplexId wrapZBack = 0;
  if(p[1] == 0)
    wrapYTop = 6 * wrap_[1];
  if(p[2] == 0)
    wrapZBack = 6 * wrap_[2];
  switch(id) {
    case 0:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6;
    case 1:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2;
    case 2:
      return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6 + 1
             + wrapYTop;
    case 3:
      return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 3
             + wrapZBack;
    case 4:
      return (p[2] - 1) * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
             + 4 + wrapYTop + wrapZBack;
    case 5:
      return (p[2] - 1) * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
             + 5 + wrapYTop + wrapZBack;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeStarH(const SimplexId p[3],
                                                   const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapZBack = 0;
  if(p[0] == 0)
    wrapXLeft = 6 * wrap_[0];
  if(p[2] == 0)
    wrapZBack = 6 * wrap_[2];
  switch(id) {
    case 0:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6;
    case 1:
      return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2
             + wrapZBack;
    case 2:
      return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 3
             + wrapZBack;
    case 3:
      return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
             + 4 + wrapXLeft + wrapZBack;
    case 4:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6 + 1
             + wrapXLeft;
    case 5:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6 + 5
             + wrapXLeft;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeStarP(const SimplexId p[3],
                                                   const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapYTop = 0;
  if(p[0] == 0)
    wrapXLeft = 6 * wrap_[0];
  if(p[1] == 0)
    wrapYTop = 6 * wrap_[1];
  switch(id) {
    case 0:
      return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + (p[0] - 1) * 6
             + 5 + wrapXLeft + wrapYTop;
    case 1:
      return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
             + wrapYTop;
    case 2:
      return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6 + 1
             + wrapYTop;
    case 3:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2;
    case 4:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6 + 3
             + wrapXLeft;
    case 5:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6 + 4
             + wrapXLeft;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeStarD1(const SimplexId p[3],
                                                    const int id) const {
  SimplexId wrapZBack = 0;
  if(p[2] == 0)
    wrapZBack = 6 * wrap_[2];
  switch(id) {
    case 0:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6;
    case 1:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 1;
    case 2:
      return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 3
             + wrapZBack;
    case 3:
      return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 4
             + wrapZBack;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeStarD2(const SimplexId p[3],
                                                    const int id) const {
  SimplexId wrapXLeft = 0;
  if(p[0] == 0)
    wrapXLeft = 6 * wrap_[0];
  switch(id) {
    case 0:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6;
    case 1:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2;
    case 2:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6 + 5
             + wrapXLeft;
    case 3:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6 + 4
             + wrapXLeft;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getEdgeStarD3(const SimplexId p[3],
                                                    const int id) const {
  SimplexId wrapYTop = 0;
  if(p[1] == 0)
    wrapYTop = 6 * wrap_[1];
  switch(id) {
    case 0:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2;
    case 1:
      return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 3;
    case 2:
      return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6 + 1
             + wrapYTop;
    case 3:
      return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6 + 5
             + wrapYTop;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleVertexF(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
               + wrapYBottom;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0] + 1
               + wrapXRight + wrapYBottom;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
               + wrapYBottom;
    }
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleVertexH(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1] + 1
               + wrapXRight + wrapZFront;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
               + wrapZFront;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
               + wrapZFront;
    }
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleVertexC(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
               + wrapZFront;
      case 2:
        return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
               + vshift_[0] + wrapYBottom + wrapZFront;
    }
  } else {
    switch(id) {
      case 0:
        return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
               + wrapYBottom;
      case 2:
        return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
               + vshift_[0] + wrapYBottom + wrapZFront;
    }
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleVertexD1(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
               + wrapZFront;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
               + vshift_[0] + wrapYBottom + wrapZFront;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0] + 1
               + wrapXRight + wrapYBottom;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
               + vshift_[0] + wrapYBottom + wrapZFront;
    }
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleVertexD2(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
               + vshift_[1] + 1 + wrapXRight + wrapYBottom + wrapZFront;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
               + vshift_[1] + wrapYBottom + wrapZFront;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
               + vshift_[1] + wrapYBottom + wrapZFront;
    }
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleVertexD3(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1] + 1
               + wrapXRight + wrapZFront;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
               + vshift_[1] + wrapYBottom + wrapZFront;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
               + wrapYBottom;
      case 2:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
               + vshift_[1] + wrapYBottom + wrapZFront;
    }
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleEdgeF_0(const SimplexId p[3],
                                                         const int id) const {
  switch(id) {
    case 0:
      return p[0] / 2 + p[1] * eshift_[0] + p[2] * eshift_[1];
    case 1:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2] + p[2] * eshift_[3];
    case 2:
      return esetshift_[2] + p[0] / 2 + p[1] * eshift_[6] + p[2] * eshift_[7];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleEdgeF_1(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapXRight = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  SimplexId wrapYBottom = 0;
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  switch(id) {
    case 0:
      return p[0] / 2 + (p[1] + 1) * eshift_[0] + p[2] * eshift_[1]
             + wrapYBottom;
    case 1:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2] + p[2] * eshift_[3]
             + 1 + wrapXRight;
    case 2:
      return esetshift_[2] + p[0] / 2 + p[1] * eshift_[6] + p[2] * eshift_[7];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleEdgeH_0(const SimplexId p[3],
                                                         const int id) const {
  switch(id) {
    case 0:
      return p[0] / 2 + p[1] * eshift_[0] + p[2] * eshift_[1];
    case 1:
      return esetshift_[1] + p[0] / 2 + p[1] * eshift_[4] + p[2] * eshift_[5];
    case 2:
      return esetshift_[4] + p[0] / 2 + p[1] * eshift_[10] + p[2] * eshift_[11];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleEdgeH_1(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapXRight = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  SimplexId wrapZFront = 0;
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] / 2 + p[1] * eshift_[0] + (p[2] + 1) * eshift_[1]
             + wrapZFront;
    case 1:
      return esetshift_[1] + p[0] / 2 + p[1] * eshift_[4] + p[2] * eshift_[5]
             + 1 + wrapXRight;
    case 2:
      return esetshift_[4] + p[0] / 2 + p[1] * eshift_[10] + p[2] * eshift_[11];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleEdgeC_0(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapYBottom = 0;
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2] + p[2] * eshift_[3];
    case 1:
      return esetshift_[1] + p[0] / 2 + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5] + wrapYBottom;
    case 2:
      return esetshift_[3] + p[0] / 2 + p[1] * eshift_[8] + p[2] * eshift_[9];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleEdgeC_1(const SimplexId p[3],
                                                         const int id) const {
  SimplexId wrapZFront = 0;
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3] + wrapZFront;
    case 1:
      return esetshift_[1] + p[0] / 2 + p[1] * eshift_[4] + p[2] * eshift_[5];
    case 2:
      return esetshift_[3] + p[0] / 2 + p[1] * eshift_[8] + p[2] * eshift_[9];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleEdgeD1_0(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapXRight = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  SimplexId wrapYBottom = 0;
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2] + p[2] * eshift_[3]
             + 1 + wrapXRight;
    case 1:
      return esetshift_[4] + p[0] / 2 + (p[1] + 1) * eshift_[10]
             + p[2] * eshift_[11] + wrapYBottom;
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleEdgeD1_1(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapZFront = 0;
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3] + wrapZFront;
    case 1:
      return esetshift_[4] + p[0] / 2 + p[1] * eshift_[10] + p[2] * eshift_[11];
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleEdgeD2_0(const SimplexId p[3],
                                                          const int id) const {
  switch(id) {
    case 0:
      return p[0] / 2 + p[1] * eshift_[0] + p[2] * eshift_[1];
    case 1:
      return esetshift_[3] + p[0] / 2 + p[1] * eshift_[8] + p[2] * eshift_[9];
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleEdgeD2_1(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] / 2 + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1]
             + wrapYBottom + wrapZFront;
    case 1:
      return esetshift_[3] + p[0] / 2 + p[1] * eshift_[8] + p[2] * eshift_[9]
             + 1 + wrapXRight;
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleEdgeD3_0(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapYBottom = 0;
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] / 2 + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5] + wrapYBottom;
    case 1:
      return esetshift_[2] + p[0] / 2 + p[1] * eshift_[6] + p[2] * eshift_[7];
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleEdgeD3_1(const SimplexId p[3],
                                                          const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] / 2 + p[1] * eshift_[4] + p[2] * eshift_[5]
             + 1 + wrapXRight;
    case 1:
      return esetshift_[2] + p[0] / 2 + p[1] * eshift_[6]
             + (p[2] + 1) * eshift_[7] + wrapZFront;
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleLinkF(const SimplexId p[3],
                                                       const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZBack = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == 0)
    wrapZBack = wrap_[2];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1]
             + wrapYBottom + wrapZFront;
    case 1:
      return p[0] / 2 + p[1] * vshift_[0] + (p[2] - 1) * vshift_[1] + 1
             + wrapXRight + wrapZBack;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleLinkH(const SimplexId p[3],
                                                       const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYTop = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == 0)
    wrapYTop = wrap_[1];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1]
             + wrapYBottom + wrapZFront;
    case 1:
      return p[0] / 2 + (p[1] - 1) * vshift_[0] + p[2] * vshift_[1] + 1
             + wrapXRight + wrapYTop;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleLinkC(const SimplexId p[3],
                                                       const int id) const {
  SimplexId wrapXLeft = 0;
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == 0)
    wrapXLeft = wrap_[0];
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1 + wrapXRight;
    case 1:
      return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1] - 1
             + wrapXLeft + wrapYBottom + wrapZFront;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleLinkD1(const SimplexId p[3],
                                                        const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + (p[2] + 1) * vshift_[1] + 1
               + wrapXRight + wrapZFront;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + p[2] * vshift_[1]
               + wrapYBottom;
      case 1:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1] + 1
               + wrapXRight + wrapYBottom + wrapZFront;
    }
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleLinkD2(const SimplexId p[3],
                                                        const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight + wrapYBottom;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + (p[2] + 1) * vshift_[1] + 1
               + wrapXRight + wrapZFront;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + p[2] * vshift_[1]
               + wrapYBottom;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + (p[2] + 1) * vshift_[1]
               + wrapZFront;
    }
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleLinkD3(const SimplexId p[3],
                                                        const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] / 2 == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + (p[2] + 1) * vshift_[1]
               + wrapZFront;
      case 1:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1] + 1
               + wrapXRight + wrapYBottom + wrapZFront;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + p[2] * vshift_[1] + 1
               + wrapXRight + wrapYBottom;
    }
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleStarF(const SimplexId p[3],
                                                       const int id) const {
  SimplexId wrapZBack = 0;
  if(p[2] == 0)
    wrapZBack = 6 * wrap_[2];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 1;
      case 1:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
               + 4 + wrapZBack;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1];
      case 1:
        return p[0] * 3 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1] + 3
               + wrapZBack;
    }
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleStarH(const SimplexId p[3],
                                                       const int id) const {
  SimplexId wrapYTop = 0;
  if(p[1] == 0)
    wrapYTop = 6 * wrap_[1];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 3;
      case 1:
        return (p[0] - 1) * 3 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
               + 5 + wrapYTop;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 2;
      case 1:
        return p[0] * 3 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1] + 1
               + wrapYTop;
    }
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleStarC(const SimplexId p[3],
                                                       const int id) const {
  SimplexId wrapXLeft = 0;
  if(p[0] < 2)
    wrapXLeft = 6 * wrap_[0];
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 2;
      case 1:
        return (p[0] / 2 - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
               + 4 + wrapXLeft;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1];
      case 1:
        return (p[0] / 2 - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
               + 5 + wrapXLeft;
    }
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleStarD1(const SimplexId p[3],
                                                        const int id) const {
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 2;
      case 1:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 3;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 1;
      case 1:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 5;
    }
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleStarD2(const SimplexId p[3],
                                                        const int id) const {
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 5;
      case 1:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 4;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1];
      case 1:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 2;
    }
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTriangleStarD3(const SimplexId p[3],
                                                        const int id) const {
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 3;
      case 1:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 4;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1];
      case 1:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 1;
    }
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronVertexABCG(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1];
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1 + wrapXRight;
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + wrapYBottom;
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + wrapYBottom + wrapZFront;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronVertexBCDG(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1 + wrapXRight;
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + wrapYBottom;
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0] + 1
             + wrapXRight + wrapYBottom;
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + wrapYBottom + wrapZFront;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronVertexABEG(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1];
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1 + wrapXRight;
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
             + wrapZFront;
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + wrapYBottom + wrapZFront;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronVertexBEFG(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1 + wrapXRight;
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
             + wrapZFront;
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1] + 1
             + wrapXRight + wrapZFront;
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + wrapYBottom + wrapZFront;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronVertexBFGH(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1 + wrapXRight;
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1] + 1
             + wrapXRight + wrapZFront;
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + wrapYBottom + wrapZFront;
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + 1 + wrapXRight + wrapYBottom + wrapZFront;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronVertexBDGH(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1 + wrapXRight;
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0] + 1
             + wrapXRight + wrapYBottom;
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + wrapYBottom + wrapZFront;
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + 1 + wrapXRight + wrapYBottom + wrapZFront;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronEdgeABCG(
    const SimplexId p[3], const int id) const {
  SimplexId wrapYBottom = 0;
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1];
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] + p[2] * eshift_[3];
    case 2:
      return esetshift_[1] + p[0] + (p[1] + 1) * eshift_[4] + p[2] * eshift_[5]
             + wrapYBottom;
    case 3:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7];
    case 4:
      return esetshift_[3] + p[0] + p[1] * eshift_[8] + p[2] * eshift_[9];
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronEdgeBCDG(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  switch(id) {
    case 0:
      return p[0] + (p[1] + 1) * eshift_[0] + p[2] * eshift_[1] + wrapYBottom;
    case 1:
      return esetshift_[0] + (p[0] + 1) + p[1] * eshift_[2] + p[2] * eshift_[3]
             + wrapXRight;
    case 2:
      return esetshift_[1] + p[0] + (p[1] + 1) * eshift_[4] + p[2] * eshift_[5]
             + wrapYBottom;
    case 3:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7];
    case 4:
      return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
             + p[2] * eshift_[11] + wrapYBottom;
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronEdgeABEG(
    const SimplexId p[3], const int id) const {
  SimplexId wrapZFront = 0;
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1];
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] + (p[2] + 1) * eshift_[3]
             + wrapZFront;
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] + p[2] * eshift_[5];
    case 3:
      return esetshift_[3] + p[0] + p[1] * eshift_[8] + p[2] * eshift_[9];
    case 4:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11];
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronEdgeBEFG(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + (p[2] + 1) * eshift_[1] + wrapZFront;
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] + (p[2] + 1) * eshift_[3]
             + wrapZFront;
    case 2:
      return esetshift_[1] + (p[0] + 1) + p[1] * eshift_[4] + p[2] * eshift_[5]
             + wrapXRight;
    case 3:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + (p[2] + 1) * eshift_[7]
             + wrapZFront;
    case 4:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11];
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronEdgeBFGH(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1]
             + wrapYBottom + wrapZFront;
    case 1:
      return esetshift_[0] + (p[0] + 1) + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3] + wrapXRight + wrapZFront;
    case 2:
      return esetshift_[1] + (p[0] + 1) + p[1] * eshift_[4] + p[2] * eshift_[5]
             + wrapXRight;
    case 3:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + (p[2] + 1) * eshift_[7]
             + wrapZFront;
    case 4:
      return esetshift_[3] + (p[0] + 1) + p[1] * eshift_[8] + p[2] * eshift_[9]
             + wrapXRight;
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronEdgeBDGH(
    const SimplexId p[3], const int id) const {
  SimplexId wrapXRight = 0;
  SimplexId wrapYBottom = 0;
  SimplexId wrapZFront = 0;
  if(p[0] == nbvoxels_[0])
    wrapXRight = -wrap_[0];
  if(p[1] == nbvoxels_[1])
    wrapYBottom = -wrap_[1];
  if(p[2] == nbvoxels_[2])
    wrapZFront = -wrap_[2];
  switch(id) {
    case 0:
      return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1]
             + wrapYBottom + wrapZFront;
    case 1:
      return esetshift_[0] + (p[0] + 1) + p[1] * eshift_[2] + p[2] * eshift_[3]
             + wrapXRight;
    case 2:
      return esetshift_[1] + (p[0] + 1) + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5] + wrapXRight + wrapYBottom;
    case 3:
      return esetshift_[3] + (p[0] + 1) + p[1] * eshift_[8] + p[2] * eshift_[9]
             + wrapXRight;
    case 4:
      return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
             + p[2] * eshift_[11] + wrapYBottom;
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronTriangleABCG(
    const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 1:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
    case 3:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11];
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronTriangleBCDG(
    const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11];
    case 2:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7];
    case 3:
      return (p[1] < nbvoxels_[1])
               ? tsetshift_[0] + p[0] * 2 + (p[1] + 1) * tshift_[2]
                   + p[2] * tshift_[3]
               : tsetshift_[0] + p[0] * 2 + (p[1] + 1) * tshift_[2]
                   + p[2] * tshift_[3] - wrap_[1] * 2;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronTriangleABEG(
    const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 1:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 3:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronTriangleBEFG(
    const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1;
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3]
             + 1;
    case 2:
      return (p[2] < nbvoxels_[2])
               ? p[0] * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
               : p[0] * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
                   - wrap_[2] * 2;
    case 3:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11]
             + 1;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronTriangleBFGH(
    const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9]
             + 1;
    case 1:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11]
             + 1;
    case 2:
      return (p[0] < nbvoxels_[0])
               ? tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
                   + p[2] * tshift_[5] + 3
               : tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
                   + p[2] * tshift_[5] + 3 - wrap_[0] * 2;
    case 3:
      return (p[2] < nbvoxels_[2])
               ? p[0] * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1] + 1
               : p[0] * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1] + 1
                   - wrap_[2] * 2;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronTriangleBDGH(
    const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return (p[0] < nbvoxels_[0])
               ? tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
                   + p[2] * tshift_[5] + 2
               : tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
                   + p[2] * tshift_[5] + 2 - wrap_[0] * 2;
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7];
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9]
             + 1;
    case 3:
      return (p[1] < nbvoxels_[1])
               ? tsetshift_[0] + p[0] * 2 + (p[1] + 1) * tshift_[2]
                   + p[2] * tshift_[3] + 1
               : tsetshift_[0] + p[0] * 2 + (p[1] + 1) * tshift_[2]
                   + p[2] * tshift_[3] + 1 - wrap_[1] * 2;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronNeighborABCG(
    const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t + 1;
    case 1:
      return t + 2;
    case 2:
      return p[0] > 0 ? t - 1 : t - 1 + wrap_[0] * 6;
    case 3:
      return p[2] > 0 ? t - tetshift_[1] + 3
                      : t - tetshift_[1] + 3 + wrap_[2] * 6;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronNeighborBCDG(
    const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 1;
    case 1:
      return t + 4;
    case 2:
      return p[2] > 0 ? t - tetshift_[1] + 3
                      : t - tetshift_[1] + 3 + wrap_[2] * 6;
    case 3:
      return p[1] < nbvoxels_[1] ? t + tetshift_[0] + 1
                                 : t + tetshift_[0] + 1 - wrap_[1] * 6;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronNeighborABEG(
    const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 2;
    case 1:
      return t + 1;
    case 2:
      return p[0] > 0 ? t - 4 : t - 4 + wrap_[0] * 6;
    case 3:
      return p[1] > 0 ? t - tetshift_[0] - 1
                      : t - tetshift_[0] - 1 + wrap_[1] * 6;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronNeighborBEFG(
    const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 1;
    case 1:
      return t + 1;
    case 2:
      return p[1] > 0 ? t - tetshift_[0] + 2
                      : t - tetshift_[0] + 2 + wrap_[1] * 6;
    case 3:
      return p[2] < nbvoxels_[2] ? t + tetshift_[1] - 3
                                 : t + tetshift_[1] - 3 - wrap_[2] * 6;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronNeighborBFGH(
    const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 1;
    case 1:
      return t + 1;
    case 2:
      return p[0] < nbvoxels_[0] ? t + 4 : t + 4 - wrap_[0] * 6;
    case 3:
      return p[2] < nbvoxels_[2] ? t + tetshift_[1] - 3
                                 : t + tetshift_[1] - 3 - wrap_[2] * 6;
  }
  return -1;
}

template <size_t card>
inline ttk::SimplexId
  ttk::PeriodicImplicitTriangulation<card>::getTetrahedronNeighborBDGH(
    const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 1;
    case 1:
      return t - 4;
    case 2:
      return p[0] < nbvoxels_[0] ? t + 1 : t + 1 - wrap_[0] * 6;
    case 3:
      return p[1] < nbvoxels_[1] ? t + tetshift_[0] - 2
                                 : t + tetshift_[0] - 2 - wrap_[1] * 6;
  }
  return -1;
}

#include <PeriodicPreconditions.h>


// -- was in cpptemplate<size_t card>
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
      localBBox[0] = max(p[0], static_cast<ttk::SimplexId>(0));
    }
    if(p[0] > localBBox[1]) {
      localBBox[1] = min(p[0], dims[0]);
    }
    if(p[1] < localBBox[2]) {
      localBBox[2] = max(p[1], static_cast<ttk::SimplexId>(0));
    }
    if(p[1] > localBBox[3]) {
      localBBox[3] = min(p[1], dims[1]);
    }
    if(p[2] < localBBox[4]) {
      localBBox[4] = max(p[2], static_cast<ttk::SimplexId>(0));
    }
    if(p[2] > localBBox[5]) {
      localBBox[5] = min(p[2], dims[2]);
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
       && pGlobal[i] == this->metaGrid_->this->dimensions_[i] - 1) {
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
/*
template class ttk::PeriodicImplicitTriangulationCRTP<0,
  ttk::PeriodicWithPreconditions<0>>;
template class ttk::PeriodicImplicitTriangulationCRTP<0,
  ttk::PeriodicNoPreconditions<0>>;*/
/// @endcond
