// Todo peut-être le templater ?
/// \ingroup base
/// \class ttk::MultiresTopology
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date 2022.
///
/// \brief TTK processing package for progressive Topological Data Analysis
///
/// This package introduces a multiresolution hierarchical representation of the
/// data.
///
/// It allows the definition of efficient algorithms for approximate or
/// progressive computation in TDA.
//
/// It is used in ttk::ProgressiveTopology for the definition of
/// efficient progressive algorithms for the computation of critical points and
/// extremum-saddle persistence diagrams.
///
/// It is also used in ttk::ApproximateTopology for the approximate computation
/// of the extremum-saddle persistence diagrams with guarantees.
///
/// \b Related \b publications \n
/// "A Progressive Approach to Scalar Field Topology" \n
/// Jules Vidal, Pierre Guillou, Julien Tierny\n
/// IEEE Transactions on Visualization and Computer Graphics, 2021
///
/// "Fast Approximation of Persistence Diagrams with Guarantees" \n
/// Jules Vidal, Julien Tierny\n
/// IEEE Symposium on Large Data Visualization and Analysis (LDAV), 2021
///
///
/// \sa ProgressiveTopology
/// \sa ApproximateTopology

#pragma once

// base code includes
#include <DynamicTree.h>
#include <ImplicitTriangulation.h>
#include <MultiresTriangulation.h>

#include <limits>
#include <tuple>

namespace ttk {

  using triplet = std::tuple<ttk::SimplexId, ttk::SimplexId, ttk::SimplexId>;
  using polarity = unsigned char;

  class MultiresTopology : public Debug {

  public:
    struct PersistencePair {
      /** first (lower) vertex id */
      ttk::SimplexId birth{};
      /** second (higher) vertex id */
      ttk::SimplexId death{};
      /** pair type (min-saddle: 0, saddle-saddle: 1, saddle-max: 2) */
      ttk::SimplexId pairType{};

      PersistencePair() = default;
      PersistencePair(const SimplexId b,
                      const SimplexId d,
                      const SimplexId pType)
        : birth{b}, death{d}, pairType{pType} {
      }
    };

    MultiresTopology() {
      this->setDebugMsgPrefix("MultiresTopology");
    }

    inline void setupTriangulation(AbstractTriangulation *const data) {
      int card = data->getDimensionality();
      if (card == 0) {
        triangulation0_ = static_cast<ttk::ImplicitTriangulation<0>*>(data);
        multiresTriangulation0_.setTriangulation(triangulation0_);
        dimensionality_ = 0;
      }
      else if (card == 1) {
        triangulation1_ = static_cast<ttk::ImplicitTriangulation<1>*>(data);
        multiresTriangulation1_.setTriangulation(triangulation1_);
        dimensionality_ = 1;
      }
      else if (card == 2) {
        triangulation2_ = static_cast<ttk::ImplicitTriangulation<2>*>(data);
        multiresTriangulation2_.setTriangulation(triangulation2_);
        dimensionality_ = 2;
      }
      else if (card == 3) {
        triangulation3_ = static_cast<ttk::ImplicitTriangulation<3>*>(data);
        multiresTriangulation3_.setTriangulation(triangulation3_);
        dimensionality_ = 3;
      }
    }

    virtual void setStartingDecimationLevel(int data) {
      startingDecimationLevel_ = std::max(data, 0);
    }
    virtual void setStoppingDecimationLevel(int data) {
      stoppingDecimationLevel_ = std::max(data, 0);
    }
    virtual void setPreallocateMemory(const bool b) {
      this->preallocateMemory_ = b;
    }
    inline int getStoppingDecimationLevel() {
      return this->stoppingDecimationLevel_;
    }

    void setStartingResolutionLevel(int rl) {
      if (dimensionality_ == 0) {
        this->setStartingDecimationLevel(multiresTriangulation0_.RL_to_DL(rl));
      }
      else if (dimensionality_ == 1) {
        this->setStartingDecimationLevel(multiresTriangulation1_.RL_to_DL(rl));
      }
      else if (dimensionality_ == 2) {
        this->setStartingDecimationLevel(multiresTriangulation2_.RL_to_DL(rl));
      }
      else if (dimensionality_ == 3) {
        this->setStartingDecimationLevel(multiresTriangulation3_.RL_to_DL(rl));
      }
    }

    void setStoppingResolutionLevel(int rl) {
      if (dimensionality_ == 0) {
        this->setStoppingDecimationLevel(multiresTriangulation0_.RL_to_DL(rl));
      }
      else if (dimensionality_ == 1) {
        this->setStoppingDecimationLevel(multiresTriangulation1_.RL_to_DL(rl));
      }
      else if (dimensionality_ == 2) {
        this->setStoppingDecimationLevel(multiresTriangulation2_.RL_to_DL(rl));
      }
      else if (dimensionality_ == 3) {
        this->setStoppingDecimationLevel(multiresTriangulation3_.RL_to_DL(rl));
      }
    }


  protected:
    int dimensionality_ = -1;
    // maximum link size in 3D
    static const size_t nLink_ = 27;

    using VLBoundaryType
      = std::array<std::vector<std::pair<SimplexId, SimplexId>>, nLink_>;

    void buildVertexLinkByBoundary(const SimplexId vertexId,
                                   VLBoundaryType &vlbt) const;



    template <size_t card>
    void buildVertexLinkByBoundary(const SimplexId vertexId,
                                   VLBoundaryType &vlbt) const;

    char getCriticalTypeFromLink(
      const std::vector<std::pair<polarity, polarity>> &vlp,
      DynamicTree &link) const;

    void getValencesFromLink(
      const SimplexId vertexId,
      const std::vector<std::pair<polarity, polarity>> &vlp,
      DynamicTree &link,
      std::vector<polarity> &toPropagateMin,
      std::vector<polarity> &toPropagateMax,
      std::vector<std::vector<SimplexId>> &saddleCCMin,
      std::vector<std::vector<SimplexId>> &saddleCCMax) const;

    void getTripletsFromSaddles(
      const SimplexId vertexId,
      std::vector<triplet> &triplets,
      const std::vector<std::vector<SimplexId>> &vertexReps) const;

    void updateLinkPolarityPonctual(
      std::vector<std::pair<polarity, polarity>> &vlp) const;

    inline void setDecimationLevel(int decimationLevel) {
      if (this->dimensionality_ == 0) {
        multiresTriangulation0_.setDecimationLevel(decimationLevel);
      }
      else if (this->dimensionality_ == 1) {
        multiresTriangulation1_.setDecimationLevel(decimationLevel);
      }
      else if (this->dimensionality_ == 2) {
        multiresTriangulation2_.setDecimationLevel(decimationLevel);
      }
      else if (this->dimensionality_ == 3) {
        multiresTriangulation3_.setDecimationLevel(decimationLevel);
      }
    }

    inline int getVertexNumber(void) const {
      if (this->dimensionality_ == 0) {
        multiresTriangulation0_.getVertexNumber();
      }
      else if (this->dimensionality_ == 1) {
        multiresTriangulation1_.getVertexNumber();
      }
      else if (this->dimensionality_ == 2) {
        multiresTriangulation2_.getVertexNumber();
      }
      else if (this->dimensionality_ == 3) {
        multiresTriangulation3_.getVertexNumber();
      }
      
      this->printErr("Dimensionality should be between 0 and 3 inclusive.");
      return -1;
    }

    inline int DL_to_RL(int dl) {
      if (this->dimensionality_ == 0) {
        return multiresTriangulation0_.DL_to_RL(dl);
      }
      else if (this->dimensionality_ == 1) {
        return multiresTriangulation1_.DL_to_RL(dl);
      }
      else if (this->dimensionality_ == 2) {
        return multiresTriangulation2_.DL_to_RL(dl);
      }
      else if (this->dimensionality_ == 3) {
        return multiresTriangulation3_.DL_to_RL(dl);
      }

      this->printErr("Dimensionality should be between 0 and 3 inclusive.");
      return -1;
    }

    inline int getVertexNeighbor(const SimplexId &vertexId,
        const int &localNeighborId,
        SimplexId &neighborId) const {
      if (this->dimensionality_ == 0) {
        return multiresTriangulation0_.getVertexNeighbor(vertexId, localNeighborId, neighborId);
      }
      else if (this->dimensionality_ == 1) {
        return multiresTriangulation1_.getVertexNeighbor(vertexId, localNeighborId, neighborId);
      }
      else if (this->dimensionality_ == 2) {
        return multiresTriangulation2_.getVertexNeighbor(vertexId, localNeighborId, neighborId);
      }
      else if (this->dimensionality_ == 3) {
        return multiresTriangulation3_.getVertexNeighbor(vertexId, localNeighborId, neighborId);
      }

      this->printErr("Dimensionality should be between 0 and 3 inclusive.");
      return -1;
    }

    inline int getDecimatedVertexNumber() const {
      if (this->dimensionality_ == 0) {
        return multiresTriangulation0_.getDecimatedVertexNumber();
      }
      else if (this->dimensionality_ == 1) {
        return multiresTriangulation1_.getDecimatedVertexNumber();
      }
      else if (this->dimensionality_ == 2) {
        return multiresTriangulation2_.getDecimatedVertexNumber();
      }
      else if (this->dimensionality_ == 3) {
        return multiresTriangulation3_.getDecimatedVertexNumber();
      }
      
      this->printErr("Dimensionality should be between 0 and 3 inclusive.");
      return -1;
    }

    inline SimplexId localToGlobalVertexId(const SimplexId localId) const {
      if (this->dimensionality_ == 0) {
        return multiresTriangulation0_.localToGlobalVertexId(localId);
      }
      else if (this->dimensionality_ == 1) {
        return multiresTriangulation1_.localToGlobalVertexId(localId);
      }
      else if (this->dimensionality_ == 2) {
        return multiresTriangulation2_.localToGlobalVertexId(localId);
      }
      else if (this->dimensionality_ == 3) {
        return multiresTriangulation3_.localToGlobalVertexId(localId);
      }
      
      this->printErr("Dimensionality should be between 0 and 3 inclusive.");
      return -1;
    }


    inline SimplexId getVertexNeighborNumber(const SimplexId &vertexId) const {
      if (this->dimensionality_ == 0) {
        return multiresTriangulation0_.getVertexNeighborNumber(vertexId);
      }
      else if (this->dimensionality_ == 1) {
        return multiresTriangulation1_.getVertexNeighborNumber(vertexId);
      }
      else if (this->dimensionality_ == 2) {
        return multiresTriangulation2_.getVertexNeighborNumber(vertexId);
      }
      else if (this->dimensionality_ == 3) {
        return multiresTriangulation3_.getVertexNeighborNumber(vertexId);
      }

      this->printErr("Dimensionality should be between 0 and 3 inclusive.");
      return -1;
    }

    inline SimplexId getVertexBoundaryIndex(const SimplexId simplexId) const {
      if (this->dimensionality_ == 0) {
        return multiresTriangulation0_.getVertexBoundaryIndex(simplexId);
      }
      else if (this->dimensionality_ == 1) {
        return multiresTriangulation1_.getVertexBoundaryIndex(simplexId);
      }
      else if (this->dimensionality_ == 2) {
        return multiresTriangulation2_.getVertexBoundaryIndex(simplexId);
      }
      else if (this->dimensionality_ == 3) {
        return multiresTriangulation3_.getVertexBoundaryIndex(simplexId);
      }

      this->printErr("Dimensionality should be between 0 and 3 inclusive.");
      return -1;
    }

    inline bool isVertexOnTriangulationBoundary(const SimplexId lvid) const {
      if (this->dimensionality_ == 0) {
        return triangulation0_->isVertexOnBoundary(lvid);
      }
      else if (this->dimensionality_ == 1) {
        return triangulation1_->isVertexOnBoundary(lvid);
      }
      else if (this->dimensionality_ == 2) {
        return triangulation2_->isVertexOnBoundary(lvid);
      }
      else if (this->dimensionality_ == 3) {
        return triangulation3_->isVertexOnBoundary(lvid);
      }
      
      this->printErr("Dimensionality should be between 0 and 3 inclusive.");
      return false;
    }


    inline void getVertexNeighborAtDecimation(const SimplexId &vertexId,   
                                            const int &localNeighborId,    
                                            SimplexId &neighborId,
                                            int decimation) const {
      if (this->dimensionality_ == 0) {
        multiresTriangulation0_.getVertexNeighborAtDecimation(vertexId, localNeighborId, neighborId, decimation);
      }
      else if (this->dimensionality_ == 1) {
        multiresTriangulation1_.getVertexNeighborAtDecimation(vertexId, localNeighborId, neighborId, decimation);
      }
      else if (this->dimensionality_ == 2) {
        multiresTriangulation2_.getVertexNeighborAtDecimation(vertexId, localNeighborId, neighborId, decimation);
      }
      else if (this->dimensionality_ == 3) {
        multiresTriangulation3_.getVertexNeighborAtDecimation(vertexId, localNeighborId, neighborId, decimation);
      }
    }


    /*
 
    inline void tata(int toto) {
      if (this->dimensionality_ == 0) {
        multiresTriangulation0_.tata(toto);
      }
      else if (this->dimensionality_ == 1) {
        multiresTriangulation1_.tata(toto);
      }
      else if (this->dimensionality_ == 2) {
        multiresTriangulation2_.tata(toto);
      }
      else if (this->dimensionality_ == 3) {
        multiresTriangulation3_.tata(toto);
      }
    }

*/ 

    //getTriangulation ?

    inline void findBoundaryRepresentatives(std::vector<SimplexId> &boundaryRepresentatives) {
      if (this->dimensionality_ == 0) {
        multiresTriangulation0_.findBoundaryRepresentatives(boundaryRepresentatives);
      }
      else if (this->dimensionality_ == 1) {
        multiresTriangulation1_.findBoundaryRepresentatives(boundaryRepresentatives);
      }
      else if (this->dimensionality_ == 2) {
        multiresTriangulation1_.findBoundaryRepresentatives(boundaryRepresentatives);
      }
      else if (this->dimensionality_ == 3) {
        multiresTriangulation1_.findBoundaryRepresentatives(boundaryRepresentatives);
      }
    }


    std::string resolutionInfoString();

    ImplicitTriangulation<0> *triangulation0_{};
    ImplicitTriangulation<1> *triangulation1_{};
    ImplicitTriangulation<2> *triangulation2_{};
    ImplicitTriangulation<3> *triangulation3_{};

    MultiresTriangulation<0> multiresTriangulation0_{};
    MultiresTriangulation<1> multiresTriangulation1_{};
    MultiresTriangulation<2> multiresTriangulation2_{};
    MultiresTriangulation<3> multiresTriangulation3_{};

    // store the two global extrema extracted from the whole dataset vertices
    // sorting operation
    mutable SimplexId globalMax_{}, globalMin_{};

    int decimationLevel_{};
    int startingDecimationLevel_{};
    int stoppingDecimationLevel_{};
    bool preallocateMemory_{true};
    std::vector<PersistencePair> CTDiagram_{};
  };
} // namespace ttk
