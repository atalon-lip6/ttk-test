/// \ingroup base
/// \class ttk::ExplicitTriangulation
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date January 2016.
///
/// \brief ExplicitTriangulation is a class that provides time efficient
/// traversal methods on triangulations of piecewise linear manifolds.
/// \sa Triangulation

#pragma once

// base code includes
#include <AbstractTriangulation.h>
#include <CellArray.h>
#include <FlatJaggedArray.h>

#include <OneSkeleton.h>
#include <ThreeSkeleton.h>
#include <TwoSkeleton.h>
#include <ZeroSkeleton.h>

#include <cstring>
#include <numeric>

#include <memory>

namespace ttk {

  template <size_t card>
  class ExplicitTriangulation final : public AbstractTriangulation {

  public:
    ExplicitTriangulation();

    ~ExplicitTriangulation() override;

    ExplicitTriangulation(const ExplicitTriangulation &) = default;
    ExplicitTriangulation(ExplicitTriangulation &&) = default;
    ExplicitTriangulation &operator=(const ExplicitTriangulation &) = default;
    ExplicitTriangulation &operator=(ExplicitTriangulation &&) = default;

    int clear();

    size_t footprint(size_t size = 0) const;


    inline int getCellEdge(const SimplexId &cellId,
                                   const int &localEdgeId,
                                   SimplexId &edgeId) const override {
      #ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      edgeId = -1;

      if(!hasPreconditionedCellEdges())
        return -1;
#endif

      if constexpr(card == 1)
        return getCellNeighbor(cellId, localEdgeId, edgeId);

      else if constexpr(card == 2)
        return getTriangleEdgeInternal(cellId, localEdgeId, edgeId);

      return getCellEdgeInternal(cellId, localEdgeId, edgeId);
    }

    inline int getCellEdgeInternal(const SimplexId &cellId,
                                   const int &localEdgeId,
                                   SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)tetraEdgeList_.size()))
        return -1;
      if((localEdgeId < 0)
         || (localEdgeId >= (SimplexId)tetraEdgeList_[cellId].size()))
        return -2;
#endif
      edgeId = tetraEdgeList_[cellId][localEdgeId];
      return 0;
    }


    inline SimplexId
      getCellEdgeNumberInternal(const SimplexId &cellId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)tetraEdgeList_.size()))
        return -1;
#endif
      return tetraEdgeList_[cellId].size();
    }

    template <std::size_t N>
    inline void
      convertToVector(const std::vector<std::array<SimplexId, N>> &table,
                      std::vector<std::vector<SimplexId>> &vec) {
      for(size_t i = 0; i < table.size(); ++i) {
        vec[i] = {table[i].begin(), table[i].end()};
      }
    }


    inline SimplexId getCellEdgeNumber(const SimplexId &cellId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedCellEdges())
        return -1;
#endif
      if constexpr(card == 1)
        return getCellNeighborNumber(cellId);

      else if constexpr(card == 2)
        return getTriangleEdgeNumber(cellId);

      return getCellEdgeNumberInternal(cellId);
    }

    inline const std::vector<std::vector<SimplexId>> *getCellEdges() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedCellEdges())
        return nullptr;
#endif
      if constexpr(card == 1)
        return getCellNeighbors();

      else if constexpr(card == 2)
        return getTriangleEdgesInternal();

      return getCellEdgesInternal();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getCellEdgesInternal() override {

      convertToVector(tetraEdgeList_, cellEdgeVector_);
      return &cellEdgeVector_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getCellNeighbor)(
      const SimplexId &cellId,
      const int &localNeighborId,
      SimplexId &neighborId) const override {

      neighborId = cellNeighborData_[cellId][localNeighborId];
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getCellNeighborNumber)(
      const SimplexId &cellId) const override {
      return cellNeighborData_.size(cellId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getCellNeighbors)() override {
      cellNeighborData_.copyTo(cellNeighborList_);
      return &cellNeighborList_;
    }

    inline int getCellTriangle(const SimplexId &cellId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      triangleId = -1;

      if constexpr(card == 1)
        return -1;

      if(!hasPreconditionedCellTriangles())
        return -2;
#endif
      if constexpr(card == 2)
        return getCellNeighbor(cellId, localTriangleId, triangleId);

      return getCellTriangleInternal(cellId, localTriangleId, triangleId);
    }


    inline int getCellTriangleInternal(const SimplexId &cellId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)tetraTriangleList_.size()))
        return -1;
      if((localTriangleId < 0)
         || (localTriangleId >= (SimplexId)tetraTriangleList_[cellId].size()))
        return -2;
#endif
      triangleId = tetraTriangleList_[cellId][localTriangleId];

      return 0;
    }

    inline SimplexId
      getCellTriangleNumber(const SimplexId &cellId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if constexpr(card == 1)
        return -1;

      if(!hasPreconditionedCellTriangles())
        return -2;
#endif
      if constexpr(card == 2)
        return getCellNeighborNumber(cellId);

      return getCellTriangleNumberInternal(cellId);
    }


    inline SimplexId
      getCellTriangleNumberInternal(const SimplexId &cellId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((cellId < 0) || (cellId >= (SimplexId)tetraTriangleList_.size()))
        return -1;
#endif

      return tetraTriangleList_[cellId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getCellTriangles() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if constexpr(card == 1)
        return nullptr;

      if(!hasPreconditionedCellTriangles())
        return nullptr;
#endif
      if constexpr(card == 2)
        return getCellNeighbors();

      return getCellTrianglesInternal();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getCellTrianglesInternal() override {

      convertToVector(tetraTriangleList_, cellTriangleVector_);
      return &cellTriangleVector_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getCellVertex)(
      const SimplexId &cellId,
      const int &localVertexId,
      SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((!cellArray_) || (!cellNumber_))
        return -1;
#endif
      vertexId = cellArray_->getCellVertex(cellId, localVertexId);
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getCellVertexNumber)(
      const SimplexId &cellId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((!cellArray_) || (!cellNumber_))
        return -1;
#endif
      return cellArray_->getCellVertexNumber(cellId);
    }

    int TTK_TRIANGULATION_INTERNAL(getDimensionality)() const override {
      return maxCellDim_;
    }

    inline const std::vector<std::array<SimplexId, 2>> *
      TTK_TRIANGULATION_INTERNAL(getEdges)() override {
      return &edgeList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getEdgeLink)(
      const SimplexId &edgeId,
      const int &localLinkId,
      SimplexId &linkId) const override {

      linkId = edgeLinkData_[edgeId][localLinkId];
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getEdgeLinkNumber)(
      const SimplexId &edgeId) const override {
      return edgeLinkData_.size(edgeId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeLinks)() override {
      edgeLinkData_.copyTo(edgeLinkList_);
      return &edgeLinkList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getEdgeStar)(
      const SimplexId &edgeId,
      const int &localStarId,
      SimplexId &starId) const override {

      starId = edgeStarData_[edgeId][localStarId];
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getEdgeStarNumber)(
      const SimplexId &edgeId) const override {
      return edgeStarData_.size(edgeId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeStars)() override {
      edgeStarData_.copyTo(edgeStarList_);
      return &edgeStarList_;
    }

    inline int getEdgeTriangle(const SimplexId &edgeId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      triangleId = -1;

      if constexpr(card==1)
        return -1;

      if(!hasPreconditionedEdgeTriangles())
        return -2;
#endif
      if constexpr(card==2)
        return getEdgeStar(edgeId, localTriangleId, triangleId);

      return getEdgeTriangleInternal(edgeId, localTriangleId, triangleId);
    }

    inline int getEdgeTriangleInternal(const SimplexId &edgeId,
                                       const int &localTriangleId,
                                       SimplexId &triangleId) const override {
      triangleId = edgeTriangleData_[edgeId][localTriangleId];
      return 0;
    }

    inline SimplexId
      getEdgeTriangleNumber(const SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if constexpr(card == 1)
        return -1;

      if(!hasPreconditionedEdgeTriangles())
        return -2;
#endif

      if constexpr(card == 2)
        return getEdgeStarNumber(edgeId);

      return getEdgeTriangleNumberInternal(edgeId);
    }
    inline SimplexId
      getEdgeTriangleNumberInternal(const SimplexId &edgeId) const override {
      return edgeTriangleData_.size(edgeId);
    }


    inline const std::vector<std::vector<SimplexId>> *
      getEdgeTriangles() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if constexpr(card==1)
        return nullptr;

      if(!hasPreconditionedEdgeTriangles())
        return nullptr;
#endif

      if constexpr(card==2)
        return getEdgeStars();

      return getEdgeTrianglesInternal();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getEdgeTrianglesInternal() override {
      edgeTriangleData_.copyTo(edgeTriangleList_);
      return &edgeTriangleList_;
    }

    inline int getEdgeVertex(const SimplexId &edgeId,
                                     const int &localVertexId,
                                     SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      vertexId = -1;
    
      if(!hasPreconditionedEdges())
        return -2;
#endif
      if constexpr (card == 1)
        return getCellVertex(edgeId, localVertexId, vertexId);

      return getEdgeVertexInternal(edgeId, localVertexId, vertexId);
    }

    inline int getEdgeVertexInternal(const SimplexId &edgeId,
                                     const int &localVertexId,
                                     SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId >= (SimplexId)edgeList_.size()))
        return -1;
      if((localVertexId != 0) && (localVertexId != 1))
        return -2;
#endif
      if(!localVertexId)
        vertexId = edgeList_[edgeId][0];
      else
        vertexId = edgeList_[edgeId][1];
      return 0;
    }

    inline SimplexId
      TTK_TRIANGULATION_INTERNAL(getNumberOfCells)() const override {
      return cellNumber_;
    }

    inline SimplexId getNumberOfEdges() const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedEdges())
        return -2;
#endif
      if constexpr (card==1)
        return getNumberOfCells();
    
      return getNumberOfEdgesInternal(); 
    }

    inline SimplexId getNumberOfEdgesInternal() const override {
      return edgeList_.size();
    }

    inline SimplexId getNumberOfTriangles() const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if constexpr(card==1)
        return -1;

      if(!hasPreconditionedTriangles())
        return -2;
#endif
      if constexpr(card==2)
        return getNumberOfCells();

      return getNumberOfTrianglesInternal();
    }

    inline SimplexId getNumberOfTrianglesInternal() const override {
      return triangleList_.size();
    }

    inline SimplexId
      TTK_TRIANGULATION_INTERNAL(getNumberOfVertices)() const override {
      return vertexNumber_;
    }

    inline const std::vector<std::array<SimplexId, 3>> *
      TTK_TRIANGULATION_INTERNAL(getTriangles)() override {
      return &triangleList_;
    }

    inline int getTriangleEdgeInternal(const SimplexId &triangleId,
                                       const int &localEdgeId,
                                       SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0)
         || (triangleId >= (SimplexId)triangleEdgeList_.size()))
        return -1;
      if((localEdgeId < 0) || (localEdgeId > 2))
        return -2;
#endif

      edgeId = triangleEdgeList_[triangleId][localEdgeId];

      return 0;
    }

    inline SimplexId getTriangleEdgeNumberInternal(
      const SimplexId &triangleId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0)
         || (triangleId >= (SimplexId)triangleEdgeList_.size()))
        return -1;
#endif

      return triangleEdgeList_[triangleId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getTriangleEdgesInternal() override {

      convertToVector(triangleEdgeList_, triangleEdgeVector_);
      return &triangleEdgeVector_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getTriangleLink)(
      const SimplexId &triangleId,
      const int &localLinkId,
      SimplexId &linkId) const override {

      linkId = triangleLinkData_[triangleId][localLinkId];
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getTriangleLinkNumber)(
      const SimplexId &triangleId) const override {
      return triangleLinkData_[triangleId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleLinks)() override {
      triangleLinkData_.copyTo(triangleLinkList_);
      return &triangleLinkList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getTriangleStar)(
      const SimplexId &triangleId,
      const int &localStarId,
      SimplexId &starId) const override {
      starId = triangleStarData_[triangleId][localStarId];
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getTriangleStarNumber)(
      const SimplexId &triangleId) const override {
      return triangleStarData_[triangleId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleStars)() override {
      triangleStarData_.copyTo(triangleStarList_);
      return &triangleStarList_;
    }

    inline int getTriangleVertex(const SimplexId &triangleId,
                                         const int &localVertexId,
                                         SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      vertexId = -1;

      if constexpr(card==1)
        return -1;

      if(!hasPreconditionedTriangles())
        return -2;
#endif
      if constexpr(card==2)
        return getCellVertex(triangleId, localVertexId, vertexId);

      return getTriangleVertexInternal(triangleId, localVertexId, vertexId);
    }

    inline int getTriangleVertexInternal(const SimplexId &triangleId,
                                         const int &localVertexId,
                                         SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0) || (triangleId >= (SimplexId)triangleList_.size()))
        return -1;
      if((localVertexId < 0)
         || (localVertexId >= (SimplexId)triangleList_[triangleId].size()))
        return -2;
#endif
      vertexId = triangleList_[triangleId][localVertexId];
      return 0;
    }

    inline int getVertexEdge(const SimplexId &vertexId,
                                     const int &localEdgeId,
                                     SimplexId &edgeId) const override {

#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      edgeId = -1;

      if(!hasPreconditionedVertexEdges())
        return -1;
#endif
      if constexpr(card==1)
        return getVertexStar(vertexId, localEdgeId, edgeId);

      return getVertexEdgeInternal(vertexId, localEdgeId, edgeId);
    }

    inline int getVertexEdgeInternal(const SimplexId &vertexId,
                                     const int &localEdgeId,
                                     SimplexId &edgeId) const override {
      edgeId = vertexEdgeData_[vertexId][localEdgeId];
      return 0;
    }

    inline SimplexId
      getVertexEdgeNumber(const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedVertexEdges())
        return -1;
#endif
      if constexpr(card == 1)
        return getVertexStarNumber(vertexId);

      return getVertexEdgeNumberInternal(vertexId);
    }

    inline SimplexId
      getVertexEdgeNumberInternal(const SimplexId &vertexId) const override {
      return vertexEdgeData_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *getVertexEdges() override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(!hasPreconditionedVertexEdges())
        return nullptr;
#endif
      if constexpr(card==1)
        return getVertexStars();

      return getVertexEdgesInternal();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexEdgesInternal() override {
      vertexEdgeData_.copyTo(vertexEdgeList_);
      return &vertexEdgeList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexLink)(
      const SimplexId &vertexId,
      const int &localLinkId,
      SimplexId &linkId) const override {

      linkId = vertexLinkData_[vertexId][localLinkId];
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexLinkNumber)(
      const SimplexId &vertexId) const override {
      return vertexLinkData_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexLinks)() override {
      vertexLinkData_.copyTo(vertexLinkList_);
      return &vertexLinkList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexNeighbor)(
      const SimplexId &vertexId,
      const int &localNeighborId,
      SimplexId &neighborId) const override {

      neighborId = vertexNeighborData_[vertexId][localNeighborId];
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexNeighborNumber)(
      const SimplexId &vertexId) const override {
      return vertexNeighborData_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexNeighbors)() override {
      vertexNeighborData_.copyTo(vertexNeighborList_);
      return &vertexNeighborList_;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexPoint)(
      const SimplexId &vertexId, float &x, float &y, float &z) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= vertexNumber_))
        return -1;
#endif

      if(doublePrecision_) {
        x = ((const double *)pointSet_)[3 * vertexId];
        y = ((const double *)pointSet_)[3 * vertexId + 1];
        z = ((const double *)pointSet_)[3 * vertexId + 2];
      } else {
        x = ((const float *)pointSet_)[3 * vertexId];
        y = ((const float *)pointSet_)[3 * vertexId + 1];
        z = ((const float *)pointSet_)[3 * vertexId + 2];
      }

      return 0;
    }

    inline int TTK_TRIANGULATION_INTERNAL(getVertexStar)(
      const SimplexId &vertexId,
      const int &localStarId,
      SimplexId &starId) const override {
      starId = vertexStarData_[vertexId][localStarId];
      return 0;
    }

    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexStarNumber)(
      const SimplexId &vertexId) const override {
      return vertexStarData_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexStars)() override {
      vertexStarData_.copyTo(vertexStarList_);
      return &vertexStarList_;
    }

    inline int getVertexTriangle(const SimplexId &vertexId,
                                         const int &localTriangleId,
                                         SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      // initialize output variable before early return
      triangleId = -1;
   
      if constexpr(card==1)
        return -1;
    
      if(!hasPreconditionedVertexTriangles())
        return -2;
#endif
      if constexpr(card==2)
        return getVertexStar(vertexId, localTriangleId, triangleId);

      return getVertexTriangleInternal(vertexId, localTriangleId, triangleId);
    }

    inline int getVertexTriangleInternal(const SimplexId &vertexId,
                                         const int &localTriangleId,
                                         SimplexId &triangleId) const override {
      triangleId = vertexTriangleData_[vertexId][localTriangleId];
      return 0;
    }

    inline SimplexId
      getVertexTriangleNumber(const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if constexpr(card==1)
        return -1;

      if(!hasPreconditionedVertexTriangles())
        return -2;
#endif

      if constexpr(card==2)
        return getVertexStarNumber(vertexId);

      return getVertexTriangleNumberInternal(vertexId);
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexTriangles() override {
#ifndef TTK_ENABLE_KAMIKAZE
        if constexpr (card==1)
          return nullptr;

        if(!hasPreconditionedVertexTriangles())
          return nullptr;
#endif
        if constexpr (card==2)
          return getVertexStars();

        return getVertexTrianglesInternal();
      }

    inline SimplexId getVertexTriangleNumberInternal(
      const SimplexId &vertexId) const override {
      return vertexTriangleData_[vertexId].size();
    }

    inline const std::vector<std::vector<SimplexId>> *
      getVertexTrianglesInternal() override {
      vertexTriangleData_.copyTo(vertexTriangleList_);
      return &vertexTriangleList_;
    }

    inline bool TTK_TRIANGULATION_INTERNAL(isEdgeOnBoundary)(
      const SimplexId &edgeId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((edgeId < 0) || (edgeId >= (SimplexId)boundaryEdges_.size()))
        return false;
#endif
      return boundaryEdges_[edgeId];
    }

    inline bool isEmpty() const override {
      return !vertexNumber_;
    }

    inline bool TTK_TRIANGULATION_INTERNAL(isTriangleOnBoundary)(
      const SimplexId &triangleId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((triangleId < 0)
         || (triangleId >= (SimplexId)boundaryTriangles_.size()))
        return false;
#endif
      return boundaryTriangles_[triangleId];
    }

    inline bool TTK_TRIANGULATION_INTERNAL(isVertexOnBoundary)(
      const SimplexId &vertexId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if((vertexId < 0) || (vertexId >= (SimplexId)boundaryVertices_.size()))
        return false;
#endif
      return boundaryVertices_[vertexId];
    }

    int preconditionBoundaryEdgesInternal() override;
    int preconditionBoundaryTrianglesInternal() override;
    int preconditionBoundaryVerticesInternal() override;

    //TODO precondition pas templatées
    int preconditionCellEdgesInternal() override;
    int preconditionCellNeighborsInternal() override;
    int preconditionCellTrianglesInternal() override;

    int preconditionEdgesInternal() override;
    int preconditionEdgeLinksInternal() override;
    int preconditionEdgeStarsInternal() override;
    int preconditionEdgeTrianglesInternal() override;

    int preconditionTrianglesInternal() override;
    int preconditionTriangleEdgesInternal() override;
    int preconditionTriangleLinksInternal() override;
    int preconditionTriangleStarsInternal() override;

    int preconditionVertexEdgesInternal() override;
    int preconditionVertexLinksInternal() override;
    int preconditionVertexNeighborsInternal() override;
    int preconditionVertexStarsInternal() override;
    int preconditionVertexTrianglesInternal() override;

    int preconditionManifoldInternal() override;

#ifdef TTK_CELL_ARRAY_NEW
    // Layout with connectivity + offset array (new)
    inline int setInputCells(const SimplexId &cellNumber,
                             const LongSimplexId *connectivity,
                             const LongSimplexId *offset) {

      // Cell Check
      {
        if(cellNumber > 0) {
          const auto &cellDimension = offset[1] - offset[0] - 1;

          if(cellDimension < 0 || cellDimension > 3) {
            this->printErr("Unable to create triangulation for cells of "
                           "dimension 4 or higher ("
                           + std::to_string(cellDimension) + ").");
            return -1;
          }

          bool error = false;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
          for(SimplexId i = 0; i < cellNumber; i++) {
            if(offset[i + 1] - offset[i] - 1 != cellDimension) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif // TTK_ENABLE_OPENMP
              error = true;
            }
          }

          if(error) {
            this->printErr("Unable to create triangulation for "
                           "inhomogeneous\ncell dimensions.");
            return -2;
          }
        }
      }

      if(cellNumber_)
        clear();

      cellNumber_ = cellNumber;

      cellArray_
        = std::make_shared<CellArray>(connectivity, offset, cellNumber);

      // TODO: ASSUME Regular Mesh Here to compute dimension!
      if(cellNumber) {
        maxCellDim_ = cellArray_->getCellVertexNumber(0) - 1;
      }
      return 0;
    }
#else
    // Flat layout with a single array (legacy & default one)
    inline int setInputCells(const SimplexId &cellNumber,
                             const LongSimplexId *cellArray) {
      if(cellNumber_)
        clear();

      cellNumber_ = cellNumber;

      if(cellNumber) {
        // assume regular mesh here to compute dimension
        cellArray_ = std::make_shared<CellArray>(
          cellArray, cellNumber, cellArray[0] - 1);
        maxCellDim_ = cellArray[0] - 1;
      }
      return 0;
    }
#endif

    inline int setInputPoints(const SimplexId &pointNumber,
                              const void *pointSet,
                              const bool &doublePrecision = false) {

      if(vertexNumber_)
        clear();

      vertexNumber_ = pointNumber;
      pointSet_ = pointSet;
      doublePrecision_ = doublePrecision;
      return 0;
    }

    /**
     * @brief Write internal state to disk
     *
     * Use a custom binary format for fast loading
     */
    int writeToFile(std::ofstream &stream) const;
    /**
     * @brief Write internal state to disk using an ASCII format
     */
    int writeToFileASCII(std::ofstream &stream) const;
    /**
     * @brief Read from disk into internal state
     *
     * Use a custom binary format for fast loading
     */
    int readFromFile(std::ifstream &stream);

#ifdef TTK_ENABLE_MPI

    inline void setCellsGlobalIds(const LongSimplexId *const cellGid) {
      this->cellGid_ = cellGid;
    }
    inline void setVertsGlobalIds(const LongSimplexId *array) {
      this->vertGid_ = array;
    }

    inline SimplexId
      getVertexGlobalIdInternal(const SimplexId lvid) const override {
      return this->vertGid_[lvid];
    }

    inline SimplexId
      getVertexLocalIdInternal(const SimplexId gvid) const override {
      const auto it{this->vertexGidToLid_.find(gvid)};
      if(it == this->vertexGidToLid_.end()) {
        return -1;
      }
      return it->second;
    }

    inline SimplexId
      getCellGlobalIdInternal(const SimplexId lcid) const override {
      return this->cellGid_[lcid];
    }

    inline SimplexId
      getCellLocalIdInternal(const SimplexId gcid) const override {
      const auto it{this->cellGidToLid_.find(gcid)};
#ifndef TTK_ENABLE_KAMIKAZE
      if(it == this->cellGidToLid_.end()) {
        return -1;
      }
#endif // TTK_ENABLE_KAMIKAZE
      return it->second;
    }

inline SimplexId getTriangleGlobalId(const SimplexId ltid) const override {
      const auto dim{this->getDimensionality()};
#ifndef TTK_ENABLE_KAMIKAZE
      if constexpr(card != 3 && card != 2) {
        this->printErr("Only 2D and  3D datasets are supported");
        return -1;
      }
      if(!this->hasPreconditionedDistributedEdges_) {
        this->printErr("TriangleGlobalId query without pre-process!");
        this->printErr(
          "Please call preconditionDistributedTriangles() in a pre-process.");
        return -1;
      }
      if(ltid < 0 || ltid >= this->getNumberOfTriangles()) {
        return -1;
      }
#endif // TTK_ENABLE_KAMIKAZE
      if(!ttk::isRunningWithMPI()) {
        return ltid;
      }
      if constexpr(card == 3) {
        return this->getTriangleGlobalIdInternal(ltid);
      } else if constexpr(card == 2) {
        return this->getCellGlobalIdInternal(ltid);
      }
      return -1;
    }

inline SimplexId getEdgeGlobalId(const SimplexId leid) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if constexpr(card != 1 && card != 2 && card != 3) {
        this->printErr("Only 1D, 2D and 3D datasets are supported");
        return -1;
      }
      if(!this->hasPreconditionedDistributedEdges_) {
        this->printErr("EdgeGlobalId query without pre-process!");
        this->printErr(
          "Please call preconditionDistributedEdges() in a pre-process.");
        return -1;
      }
      if(leid < 0 || leid >= this->getNumberOfEdges()) {
        return -1;
      }
#endif // TTK_ENABLE_KAMIKAZE
      if(!ttk::isRunningWithMPI()) {
        return leid;
      }
      if constexpr(card == 2 || card == 3) {
        return this->getEdgeGlobalIdInternal(leid);
      } else if constexpr(card == 1) {
        return this->getCellGlobalIdInternal(leid);
      }
      return -1;
    }

    inline SimplexId
      getEdgeGlobalIdInternal(const SimplexId leid) const override {
      return this->edgeLidToGid_[leid];
    }

    virtual inline SimplexId getEdgeLocalId(const SimplexId geid) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(card != 1 && card != 2 && card != 3) {
        this->printErr("Only 1D, 2D and 3D datasets are supported");
        return -1;
      }
      if(!this->hasPreconditionedDistributedEdges_) {
        this->printErr("EdgeLocalId query without pre-process!");
        this->printErr(
          "Please call preconditionDistributedEdges() in a pre-process.");
        return -1;
      }
#endif // TTK_ENABLE_KAMIKAZE
      if(!ttk::isRunningWithMPI()) {
        return geid;
      }
      if constexpr(card == 2 || card == 3) {
        return this->getEdgeLocalIdInternal(geid);
      } else if constexpr(card == 1) {
        return this->getCellLocalIdInternal(geid);
      }
      return -1;
    }


    inline SimplexId
      getEdgeLocalIdInternal(const SimplexId geid) const override {
      const auto it = this->edgeGidToLid_.find(geid);
#ifndef TTK_ENABLE_KAMIKAZE
      if(it == this->edgeGidToLid_.end()) {
        return -1;
      }
#endif // TTK_ENABLE_KAMIKAZE
      return it->second;
    }

    inline SimplexId getTriangleLocalId(const SimplexId gtid) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if constexpr(card != 3 && card != 2) {
        this->printErr("Only 2D and 3D datasets are supported");
        return -1;
      }
      if(!this->hasPreconditionedDistributedEdges_) {
        this->printErr("TriangleLocalId query without pre-process!");
        this->printErr(
          "Please call preconditionDistributedTriangles() in a pre-process.");
        return -1;
      }
#endif // TTK_ENABLE_KAMIKAZE
      if(!ttk::isRunningWithMPI()) {
        return gtid;
      }
      if constexpr(card == 3) {
        return this->getTriangleLocalIdInternal(gtid);
      } else if constexpr (card == 2) {
        return this->getCellLocalIdInternal(gtid);
      }
      return -1;
    }

    inline SimplexId
      getTriangleGlobalIdInternal(const SimplexId ltid) const override {
      return this->triangleLidToGid_[ltid];
    }

    inline SimplexId
      getTriangleLocalIdInternal(const SimplexId gtid) const override {
      const auto it = this->triangleGidToLid_.find(gtid);
#ifndef TTK_ENABLE_KAMIKAZE
      if(it == this->triangleGidToLid_.end()) {
        return -1;
      }
#endif // TTK_ENABLE_KAMIKAZE
      return it->second;
    }

    inline int setVertexRankArray(const int *rankArray) override {
      vertexRankArray_.resize(vertexNumber_);
      std::copy(rankArray, rankArray + vertexNumber_, vertexRankArray_.begin());
      return 0;
    }

    inline int setCellRankArray(const int *rankArray) override {
      cellRankArray_.resize(cellNumber_);
      std::copy(rankArray, rankArray + cellNumber_, cellRankArray_.begin());
      return 0;
    }

    inline int getVertexRankInternal(const SimplexId lvid) const override {
      return this->vertexRankArray_[lvid];
    }

    inline std::unordered_map<SimplexId, SimplexId> &getVertexGlobalIdMap() {
      return this->vertexGidToLid_;
    }

    inline void setBoundingBox(const double *const bBox) {
      this->boundingBox_
        = {bBox[0], bBox[1], bBox[2], bBox[3], bBox[4], bBox[5]};
    }

    inline int getTriangleIncenter(const SimplexId triangleId,
                                   float incenter[3]) const {
      std::array<SimplexId, 3> vertexId{};
      if constexpr(card == 2) {
        getCellVertex(triangleId, 0, vertexId[0]);
        getCellVertex(triangleId, 1, vertexId[1]);
        getCellVertex(triangleId, 2, vertexId[2]);
      } else if constexpr(card == 3) {
        getTriangleVertex(triangleId, 0, vertexId[0]);
        getTriangleVertex(triangleId, 1, vertexId[1]);
        getTriangleVertex(triangleId, 2, vertexId[2]);
      }

      std::array<float, 9> p{};
      getVertexPoint(vertexId[0], p[0], p[1], p[2]);
      getVertexPoint(vertexId[1], p[3], p[4], p[5]);
      getVertexPoint(vertexId[2], p[6], p[7], p[8]);

      std::array<float, 3> d{};
      d[0] = Geometry::distance(&p[3], &p[6]);
      d[1] = Geometry::distance(&p[0], &p[6]);
      d[2] = Geometry::distance(&p[0], &p[3]);
      const float sum = d[0] + d[1] + d[2];

      d[0] = d[0] / sum;
      d[1] = d[1] / sum;
      d[2] = d[2] / sum;

      incenter[0] = d[0] * p[0] + d[1] * p[3] + d[2] * p[6];
      incenter[1] = d[0] * p[1] + d[1] * p[4] + d[2] * p[7];
      incenter[2] = d[0] * p[2] + d[1] * p[5] + d[2] * p[8];

      return 0;
    }

    /* TODO dimension dimensionalité de la trig ?!
    inline int getCellIncenter(const SimplexId cellid,
                               const int dim,
                               float incenter[3]) const override {
      if constexpr(card == 0) {
        getVertexPoint(cellid, incenter[0], incenter[1], incenter[2]);
      }
      else if constexpr(card == 1) {
        getEdgeIncenter(cellid, incenter);
      }
      else if constexpr(card == 3) {
        getTriangleIncenter(cellid, incenter);
      }
      else if constexpr(card == 3) {
        getTetraIncenter(cellid, incenter);
      }
      return 0;
    } */

    

  protected:
    template <typename Func0, typename Func1, typename Func2>
    int exchangeDistributedInternal(const Func0 &getGlobalSimplexId,
                                    const Func1 &storeGlobalSimplexId,
                                    const Func2 &iterCond,
                                    const int nSimplicesPerCell);

    int preconditionDistributedCellRanges();
    size_t
      computeCellRangeOffsets(std::vector<size_t> &nSimplicesPerRange) const;

    int preconditionDistributedCells() override;
    int preconditionExchangeGhostCells() override;
    int preconditionDistributedEdges() override;
    int preconditionDistributedVertices() override;
    int preconditionExchangeGhostVertices() override;
    int preconditionDistributedTriangles() override;
    int preconditionVertexRankArray();
    int preconditionCellRankArray();
    int preconditionEdgeRankArray() override;
    int preconditionTriangleRankArray() override;

    // range of (local) cells owned by the current rank that have
    // contiguous global ids (to label edges & triangles)
    struct CellRange {
      // rank-local range id
      size_t id;
      // range beginning (global cell id)
      size_t begin;
      // range end (inclusive, global cell id)
      size_t end;
      // owner rank
      size_t rank;

      static inline MPI_Datatype getMPIType() {
        MPI_Datatype res{};
        const auto cellRangeSize = sizeof(CellRange) / sizeof(size_t);
        MPI_Type_contiguous(cellRangeSize, ttk::getMPIType(size_t{}), &res);
        return res;
      }
    };
    // cell ranges per rank
    std::vector<CellRange> localCellRanges_{};
    // cell ranges from all ranks (gathered on rank 0)
    std::vector<CellRange> gatheredCellRanges_{};
    // number of CellRanges per rank
    std::vector<int> nRangesPerRank_{};

    // "GlobalCellIds" from "Generate Global Ids"
    const LongSimplexId *cellGid_{};
    // "GlobalPointIds" from "Generate Global Ids"
    const LongSimplexId *vertGid_{};

    // inverse of vertGid_
    std::unordered_map<SimplexId, SimplexId> vertexGidToLid_{};
    // inverse of cellGid_
    std::unordered_map<SimplexId, SimplexId> cellGidToLid_{};

    std::vector<SimplexId> edgeLidToGid_{};
    std::unordered_map<SimplexId, SimplexId> edgeGidToLid_{};
    std::vector<SimplexId> triangleLidToGid_{};
    std::unordered_map<SimplexId, SimplexId> triangleGidToLid_{};

    std::array<double, 6> boundingBox_{};

    std::vector<int> vertexRankArray_{};
    std::vector<int> cellRankArray_{};
    std::vector<int> edgeRankArray_{};
    std::vector<int> triangleRankArray_{};

#endif // TTK_ENABLE_MPI

  private:
    bool doublePrecision_;
    SimplexId cellNumber_, vertexNumber_;
    const void *pointSet_;
    int maxCellDim_;
    std::shared_ptr<CellArray> cellArray_;

    FlatJaggedArray vertexNeighborData_{};
    FlatJaggedArray cellNeighborData_{};
    FlatJaggedArray vertexEdgeData_{};
    FlatJaggedArray vertexTriangleData_{};
    FlatJaggedArray edgeTriangleData_{};
    FlatJaggedArray vertexStarData_{};
    FlatJaggedArray edgeStarData_{};
    FlatJaggedArray triangleStarData_{};
    FlatJaggedArray vertexLinkData_{};
    FlatJaggedArray edgeLinkData_{};
    FlatJaggedArray triangleLinkData_{};

    // Char array that identifies the file format.
    static const char *magicBytes_;
    // Current version of the file format. To be incremented at every
    // breaking change to keep backward compatibility.
    static const unsigned long formatVersion_;
  };
  //
  //

} // namespace ttk
