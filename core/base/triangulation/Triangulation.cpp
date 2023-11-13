#include <Triangulation.h>

using namespace std;
using namespace ttk;

Triangulation::Triangulation() : abstractTriangulation_{nullptr} {
  debugLevel_ = 0; // overrides the global debug level.
  gridDimensions_ = {-1, -1, -1};
  hasPeriodicBoundaries_ = false;
}

Triangulation::Triangulation(const Triangulation &rhs)
  : AbstractTriangulation(rhs), abstractTriangulation_{nullptr},
    explicitTriangulation0_{rhs.explicitTriangulation0_},
    explicitTriangulation1_{rhs.explicitTriangulation1_},
    explicitTriangulation2_{rhs.explicitTriangulation2_},
    explicitTriangulation3_{rhs.explicitTriangulation3_},
    implicitTriangulation0_{rhs.implicitTriangulation0_},
    implicitTriangulation1_{rhs.implicitTriangulation1_},
    implicitTriangulation2_{rhs.implicitTriangulation2_},
    implicitTriangulation3_{rhs.implicitTriangulation3_},
    periodicImplicitTriangulation0_{rhs.periodicImplicitTriangulation0_},
    periodicImplicitTriangulation1_{rhs.periodicImplicitTriangulation1_},
    periodicImplicitTriangulation2_{rhs.periodicImplicitTriangulation2_},
    periodicImplicitTriangulation3_{rhs.periodicImplicitTriangulation3_},
    compactTriangulation0_{rhs.compactTriangulation0_},
    compactTriangulation1_{rhs.compactTriangulation1_},
    compactTriangulation2_{rhs.compactTriangulation2_},
    compactTriangulation3_{rhs.compactTriangulation3_} {

  gridDimensions_ = rhs.gridDimensions_;
  hasPeriodicBoundaries_ = rhs.hasPeriodicBoundaries_;

  switch(rhs.getType()) {
    case Type::EXPLICIT:
      if (rhs.getDimensionality() == 0) {
        this->abstractTriangulation_ = &this->explicitTriangulation0_;
      }
      else if (rhs.getDimensionality() == 1) {
        this->abstractTriangulation_ = &this->explicitTriangulation1_;
      }
      else if (rhs.getDimensionality() == 2) {
        this->abstractTriangulation_ = &this->explicitTriangulation2_;
      }
      else if (rhs.getDimensionality() == 3) {
        this->abstractTriangulation_ = &this->explicitTriangulation3_;
      }
      else {
        this->printErr("Error: dimensionality should be between 0 and 3.");
      }
      break;
    case Type::COMPACT:
      if (rhs.getDimensionality() == 0) {
        this->abstractTriangulation_ = &this->compactTriangulation0_;
      }
      else if (rhs.getDimensionality() == 1) {
        this->abstractTriangulation_ = &this->compactTriangulation1_;
      }
      else if (rhs.getDimensionality() == 2) {
        this->abstractTriangulation_ = &this->compactTriangulation2_;
      }
      else if (rhs.getDimensionality() == 3) {
        this->abstractTriangulation_ = &this->compactTriangulation3_;
      }
      else {
        this->printErr("Error: dimensionality should be between 0 and 3.");
      }
      break;
    case Type::IMPLICIT:
      if (rhs.getDimensionality() == 0) {
        this->abstractTriangulation_ = &this->implicitTriangulation0_;
      }
      else if (rhs.getDimensionality() == 1) {
        this->abstractTriangulation_ = &this->implicitTriangulation1_;
      }
      else if (rhs.getDimensionality() == 2) {
        this->abstractTriangulation_ = &this->implicitTriangulation2_;
      }
      else if (rhs.getDimensionality() == 3) {
        this->abstractTriangulation_ = &this->implicitTriangulation3_;
      }
      else {
        this->printErr("Error: dimensionality should be between 0 and 3.");
      }
      break;
    case Type::HYBRID_IMPLICIT:
      if (rhs.getDimensionality() == 0) {
        this->abstractTriangulation_ = &this->implicitPreconditionsTriangulation0_;
      }
      else if (rhs.getDimensionality() == 1) {
        this->abstractTriangulation_ = &this->implicitPreconditionsTriangulation1_;
      }
      else if (rhs.getDimensionality() == 2) {
        this->abstractTriangulation_ = &this->implicitPreconditionsTriangulation2_;
      }
      else if (rhs.getDimensionality() == 3) {
        this->abstractTriangulation_ = &this->implicitPreconditionsTriangulation3_;
      }
      else {
        this->printErr("Error: dimensionality should be between 0 and 3.");
      }
      break;
    case Type::PERIODIC:
      if (rhs.getDimensionality() == 0) {
        this->abstractTriangulation_ = &this->periodicImplicitTriangulation0_;
      }
      else if (rhs.getDimensionality() == 1) {
        this->abstractTriangulation_ = &this->periodicImplicitTriangulation1_;
      }
      else if (rhs.getDimensionality() == 2) {
        this->abstractTriangulation_ = &this->periodicImplicitTriangulation2_;
      }
      else if (rhs.getDimensionality() == 3) {
        this->abstractTriangulation_ = &this->periodicImplicitTriangulation3_;
      }
      else {
        this->printErr("Error: dimensionality should be between 0 and 3.");
      }
      break;
    case Type::HYBRID_PERIODIC:
      if (rhs.getDimensionality() == 0) {
        this->abstractTriangulation_ = &this->periodicPreconditionsTriangulation0_;
      }
      else if (rhs.getDimensionality() == 1) {
        this->abstractTriangulation_ = &this->periodicPreconditionsTriangulation1_;
      }
      else if (rhs.getDimensionality() == 2) {
        this->abstractTriangulation_ = &this->periodicPreconditionsTriangulation2_;
      }
      else if (rhs.getDimensionality() == 3) {
        this->abstractTriangulation_ = &this->periodicPreconditionsTriangulation3_;
      }
      else {
        this->printErr("Error: dimensionality should be between 0 and 3.");
      }
      break;
  }
}

//TODO preconditions hybrid?
Triangulation::Triangulation(Triangulation &&rhs) noexcept
  : AbstractTriangulation(
    std::move(*static_cast<AbstractTriangulation *>(&rhs))),
    abstractTriangulation_{nullptr},
    explicitTriangulation0_{std::move(rhs.explicitTriangulation0_)},
    explicitTriangulation1_{std::move(rhs.explicitTriangulation1_)},
    explicitTriangulation2_{std::move(rhs.explicitTriangulation2_)},
    explicitTriangulation3_{std::move(rhs.explicitTriangulation3_)},
    implicitTriangulation0_{std::move(rhs.implicitTriangulation0_)},
    implicitTriangulation1_{std::move(rhs.implicitTriangulation1_)},
    implicitTriangulation2_{std::move(rhs.implicitTriangulation2_)},
    implicitTriangulation3_{std::move(rhs.implicitTriangulation3_)},
    periodicImplicitTriangulation0_{std::move(rhs.periodicImplicitTriangulation0_)},
    periodicImplicitTriangulation1_{std::move(rhs.periodicImplicitTriangulation1_)},
    periodicImplicitTriangulation2_{std::move(rhs.periodicImplicitTriangulation2_)},
    periodicImplicitTriangulation3_{std::move(rhs.periodicImplicitTriangulation3_)},
    compactTriangulation0_{std::move(rhs.compactTriangulation0_)},
    compactTriangulation1_{std::move(rhs.compactTriangulation1_)},
    compactTriangulation2_{std::move(rhs.compactTriangulation2_)},
    compactTriangulation3_{std::move(rhs.compactTriangulation3_)} {

  gridDimensions_ = rhs.gridDimensions_;
  hasPeriodicBoundaries_ = rhs.hasPeriodicBoundaries_;

#define repeatAssign4(a, b, c) \
  switch(c) { \
    case 0: \
    ##a = &this->##b0_;\
    case 1: \
  ##a = &this->##b1_; \
    case 2: \
  ##a = &this->##b2_; \
    case 3: \
  ##a = &this->##b3_; \

  switch(rhs.getType()) {
    case Type::EXPLICIT:
      if (rhs.abstractTriangulation_ == &rhs.explicitTriangulation0_) {
      this->abstractTriangulation_ = &this->explicitTriangulation0_;
    }
      else if (rhs.abstractTriangulation_ == &rhs.explicitTriangulation1_) {
      this->abstractTriangulation_ = &this->explicitTriangulation1_;
    }
      else if (rhs.abstractTriangulation_ == &rhs.explicitTriangulation2_) {
      this->abstractTriangulation_ = &this->explicitTriangulation2_;
    }
      else if (rhs.abstractTriangulation_ == &rhs.explicitTriangulation3_) {
      this->abstractTriangulation_ = &this->explicitTriangulation3_;
    }
      break;
    case Type::COMPACT:
      if (rhs.abstractTriangulation_ == &rhs.compactTriangulation0_) {
      this->abstractTriangulation_ = &this->compactTriangulation0_;
    }
      else if (rhs.abstractTriangulation_ == &rhs.compactTriangulation1_) {
      this->abstractTriangulation_ = &this->compactTriangulation1_;
    }
      else if (rhs.abstractTriangulation_ == &rhs.compactTriangulation2_) {
      this->abstractTriangulation_ = &this->compactTriangulation2_;
    }
      else if (rhs.abstractTriangulation_ == &rhs.compactTriangulation3_) {
      this->abstractTriangulation_ = &this->compactTriangulation3_;
    }
      break;
    case Type::IMPLICIT:
      if (rhs.abstractTriangulation_ == &rhs.implicitTriangulation0_) {
      this->abstractTriangulation_ = &this->implicitTriangulation0_;
    }
      else if (rhs.abstractTriangulation_ == &rhs.implicitTriangulation1_) {
      this->abstractTriangulation_ = &this->implicitTriangulation1_;
    }
      else if (rhs.abstractTriangulation_ == &rhs.implicitTriangulation2_) {
      this->abstractTriangulation_ = &this->implicitTriangulation2_;
    }
      else if (rhs.abstractTriangulation_ == &rhs.implicitTriangulation3_) {
      this->abstractTriangulation_ = &this->implicitTriangulation3_;
    }
      break;
    case Type::HYBRID_IMPLICIT:
      if (rhs.abstractTriangulation_ == &rhs.implicitPreconditionsTriangulation0_) {
      this->abstractTriangulation_ = &this->implicitPreconditionsTriangulation0_;
    }
      else if (rhs.abstractTriangulation_ == &rhs.implicitPreconditionsTriangulation1_) {
      this->abstractTriangulation_ = &this->implicitPreconditionsTriangulation1_;
    }
      else if (rhs.abstractTriangulation_ == &rhs.implicitPreconditionsTriangulation2_) {
      this->abstractTriangulation_ = &this->implicitPreconditionsTriangulation2_;
    }
      else if (rhs.abstractTriangulation_ == &rhs.implicitPreconditionsTriangulation3_) {
      this->abstractTriangulation_ = &this->implicitPreconditionsTriangulation3_;
    }
      break;
    case Type::PERIODIC:
      if (rhs.abstractTriangulation_ == &rhs.periodicImplicitTriangulation0_) {
      this->abstractTriangulation_ = &this->periodicImplicitTriangulation0_;
    }
      else if (rhs.abstractTriangulation_ == &rhs.periodicImplicitTriangulation1_) {
      this->abstractTriangulation_ = &this->periodicImplicitTriangulation1_;
    }
      else if (rhs.abstractTriangulation_ == &rhs.periodicImplicitTriangulation2_) {
      this->abstractTriangulation_ = &this->periodicImplicitTriangulation2_;
    }
      else if (rhs.abstractTriangulation_ == &rhs.periodicImplicitTriangulation3_) {
      this->abstractTriangulation_ = &this->periodicImplicitTriangulation3_;
    }
      break;
    case Type::HYBRID_PERIODIC:
      if (rhs.abstractTriangulation_ == &rhs.periodicPreconditionsTriangulation0_) {
      this->abstractTriangulation_ = &this->periodicPreconditionsTriangulation0_;
    }
      else if (rhs.abstractTriangulation_ == &rhs.periodicPreconditionsTriangulation1_) {
      this->abstractTriangulation_ = &this->periodicPreconditionsTriangulation1_;
    }
      else if (rhs.abstractTriangulation_ == &rhs.periodicPreconditionsTriangulation2_) {
      this->abstractTriangulation_ = &this->periodicPreconditionsTriangulation2_;
    }
      else if (rhs.abstractTriangulation_ == &rhs.periodicPreconditionsTriangulation3_) {
      this->abstractTriangulation_ = &this->periodicPreconditionsTriangulation3_;
    }
      break;
  }
  rhs.abstractTriangulation_ = nullptr;
}

Triangulation &Triangulation::operator=(const Triangulation &rhs) {
  if(this != &rhs) {
    AbstractTriangulation::operator=(rhs);
    gridDimensions_ = rhs.gridDimensions_;
    abstractTriangulation_ = nullptr;
    explicitTriangulation0_ = rhs.explicitTriangulation0_;
    explicitTriangulation1_ = rhs.explicitTriangulation1_;
    explicitTriangulation2_ = rhs.explicitTriangulation2_;
    explicitTriangulation3_ = rhs.explicitTriangulation3_;
    implicitTriangulation0_ = rhs.implicitTriangulation0_;
    implicitTriangulation1_ = rhs.implicitTriangulation1_;
    implicitTriangulation2_ = rhs.implicitTriangulation2_;
    implicitTriangulation3_ = rhs.implicitTriangulation3_;
    periodicImplicitTriangulation0_ = rhs.periodicImplicitTriangulation0_;
    periodicImplicitTriangulation1_ = rhs.periodicImplicitTriangulation1_;
    periodicImplicitTriangulation2_ = rhs.periodicImplicitTriangulation2_;
    periodicImplicitTriangulation3_ = rhs.periodicImplicitTriangulation3_;
    compactTriangulation0_ = rhs.compactTriangulation0_;
    compactTriangulation1_ = rhs.compactTriangulation1_;
    compactTriangulation2_ = rhs.compactTriangulation2_;
    compactTriangulation3_ = rhs.compactTriangulation3_;
    hasPeriodicBoundaries_ = rhs.hasPeriodicBoundaries_;

    switch(rhs.getType()) {
      case Type::EXPLICIT:
        if (rhs.abstractTriangulation_ == &rhs.explicitTriangulation0_) {
          this->abstractTriangulation_ = &this->explicitTriangulation0_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.explicitTriangulation1_) {
          this->abstractTriangulation_ = &this->explicitTriangulation1_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.explicitTriangulation2_) {
          this->abstractTriangulation_ = &this->explicitTriangulation2_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.explicitTriangulation3_) {
          this->abstractTriangulation_ = &this->explicitTriangulation3_;
        }
        break;
      case Type::COMPACT:
        if (rhs.abstractTriangulation_ == &rhs.compactTriangulation0_) {
          this->abstractTriangulation_ = &this->compactTriangulation0_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.compactTriangulation1_) {
          this->abstractTriangulation_ = &this->compactTriangulation1_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.compactTriangulation2_) {
          this->abstractTriangulation_ = &this->compactTriangulation2_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.compactTriangulation3_) {
          this->abstractTriangulation_ = &this->compactTriangulation3_;
        }
        break;
      case Type::IMPLICIT:
        if (rhs.abstractTriangulation_ == &rhs.implicitTriangulation0_) {
          this->abstractTriangulation_ = &this->implicitTriangulation0_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.implicitTriangulation1_) {
          this->abstractTriangulation_ = &this->implicitTriangulation1_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.implicitTriangulation2_) {
          this->abstractTriangulation_ = &this->implicitTriangulation2_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.implicitTriangulation3_) {
          this->abstractTriangulation_ = &this->implicitTriangulation3_;
        }
        break;
      case Type::HYBRID_IMPLICIT:
        if (rhs.abstractTriangulation_ == &rhs.implicitPreconditionsTriangulation0_) {
          this->abstractTriangulation_ = &this->implicitPreconditionsTriangulation0_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.implicitPreconditionsTriangulation1_) {
          this->abstractTriangulation_ = &this->implicitPreconditionsTriangulation1_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.implicitPreconditionsTriangulation2_) {
          this->abstractTriangulation_ = &this->implicitPreconditionsTriangulation2_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.implicitPreconditionsTriangulation3_) {
          this->abstractTriangulation_ = &this->implicitPreconditionsTriangulation3_;
        }
        break;
      case Type::PERIODIC:
        if (rhs.abstractTriangulation_ == &rhs.periodicImplicitTriangulation0_) {
          this->abstractTriangulation_ = &this->periodicImplicitTriangulation0_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.periodicImplicitTriangulation1_) {
          this->abstractTriangulation_ = &this->periodicImplicitTriangulation1_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.periodicImplicitTriangulation2_) {
          this->abstractTriangulation_ = &this->periodicImplicitTriangulation2_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.periodicImplicitTriangulation3_) {
          this->abstractTriangulation_ = &this->periodicImplicitTriangulation3_;
        }
        break;
      case Type::HYBRID_PERIODIC:
        if (rhs.abstractTriangulation_ == &rhs.periodicPreconditionsTriangulation0_) {
          this->abstractTriangulation_ = &this->periodicPreconditionsTriangulation0_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.periodicPreconditionsTriangulation1_) {
          this->abstractTriangulation_ = &this->periodicPreconditionsTriangulation1_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.periodicPreconditionsTriangulation2_) {
          this->abstractTriangulation_ = &this->periodicPreconditionsTriangulation2_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.periodicPreconditionsTriangulation3_) {
          this->abstractTriangulation_ = &this->periodicPreconditionsTriangulation3_;
        }
        break;
    }
  }
  return *this;
}

Triangulation &Triangulation::operator=(Triangulation &&rhs) noexcept {
  if(this != &rhs) {
    gridDimensions_ = rhs.gridDimensions_;
    abstractTriangulation_ = nullptr;
    explicitTriangulation0_ = std::move(rhs.explicitTriangulation0_);
    explicitTriangulation1_ = std::move(rhs.explicitTriangulation1_);
    explicitTriangulation2_ = std::move(rhs.explicitTriangulation2_);
    explicitTriangulation3_ = std::move(rhs.explicitTriangulation3_);
    implicitTriangulation0_ = std::move(rhs.implicitTriangulation0_);
    implicitTriangulation1_ = std::move(rhs.implicitTriangulation1_);
    implicitTriangulation2_ = std::move(rhs.implicitTriangulation2_);
    implicitTriangulation3_ = std::move(rhs.implicitTriangulation3_);
    periodicImplicitTriangulation0_
      = std::move(rhs.periodicImplicitTriangulation0_);
    periodicImplicitTriangulation1_
      = std::move(rhs.periodicImplicitTriangulation1_);
    periodicImplicitTriangulation2_
      = std::move(rhs.periodicImplicitTriangulation2_);
    periodicImplicitTriangulation3_
      = std::move(rhs.periodicImplicitTriangulation3_);
    compactTriangulation0_ = std::move(rhs.compactTriangulation0_);
    compactTriangulation1_ = std::move(rhs.compactTriangulation1_);
    compactTriangulation2_ = std::move(rhs.compactTriangulation2_);
    compactTriangulation3_ = std::move(rhs.compactTriangulation3_);
    hasPeriodicBoundaries_ = rhs.hasPeriodicBoundaries_;

    switch(rhs.getType()) {
      case Type::EXPLICIT:
        if (rhs.abstractTriangulation_ == &rhs.explicitTriangulation0_) {
          this->abstractTriangulation_ = &this->explicitTriangulation0_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.explicitTriangulation1_) {
          this->abstractTriangulation_ = &this->explicitTriangulation1_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.explicitTriangulation2_) {
          this->abstractTriangulation_ = &this->explicitTriangulation2_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.explicitTriangulation3_) {
          this->abstractTriangulation_ = &this->explicitTriangulation3_;
        }
        break;
      case Type::COMPACT:
        if (rhs.abstractTriangulation_ == &rhs.compactTriangulation0_) {
          this->abstractTriangulation_ = &this->compactTriangulation0_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.compactTriangulation1_) {
          this->abstractTriangulation_ = &this->compactTriangulation1_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.compactTriangulation2_) {
          this->abstractTriangulation_ = &this->compactTriangulation2_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.compactTriangulation3_) {
          this->abstractTriangulation_ = &this->compactTriangulation3_;
        }
        break;
      case Type::IMPLICIT:
        if (rhs.abstractTriangulation_ == &rhs.implicitTriangulation0_) {
          this->abstractTriangulation_ = &this->implicitTriangulation0_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.implicitTriangulation1_) {
          this->abstractTriangulation_ = &this->implicitTriangulation1_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.implicitTriangulation2_) {
          this->abstractTriangulation_ = &this->implicitTriangulation2_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.implicitTriangulation3_) {
          this->abstractTriangulation_ = &this->implicitTriangulation3_;
        }
        break;
      case Type::HYBRID_IMPLICIT:
        if (rhs.abstractTriangulation_ == &rhs.implicitPreconditionsTriangulation0_) {
          this->abstractTriangulation_ = &this->implicitPreconditionsTriangulation0_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.implicitPreconditionsTriangulation1_) {
          this->abstractTriangulation_ = &this->implicitPreconditionsTriangulation1_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.implicitPreconditionsTriangulation2_) {
          this->abstractTriangulation_ = &this->implicitPreconditionsTriangulation2_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.implicitPreconditionsTriangulation3_) {
          this->abstractTriangulation_ = &this->implicitPreconditionsTriangulation3_;
        }
        break;
      case Type::PERIODIC:
        if (rhs.abstractTriangulation_ == &rhs.periodicImplicitTriangulation0_) {
          this->abstractTriangulation_ = &this->periodicImplicitTriangulation0_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.periodicImplicitTriangulation1_) {
          this->abstractTriangulation_ = &this->periodicImplicitTriangulation1_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.periodicImplicitTriangulation2_) {
          this->abstractTriangulation_ = &this->periodicImplicitTriangulation2_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.periodicImplicitTriangulation3_) {
          this->abstractTriangulation_ = &this->periodicImplicitTriangulation3_;
        }
        break;
      case Type::HYBRID_PERIODIC:
        if (rhs.abstractTriangulation_ == &rhs.periodicPreconditionsTriangulation0_) {
          this->abstractTriangulation_ = &this->periodicPreconditionsTriangulation0_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.periodicPreconditionsTriangulation1_) {
          this->abstractTriangulation_ = &this->periodicPreconditionsTriangulation1_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.periodicPreconditionsTriangulation2_) {
          this->abstractTriangulation_ = &this->periodicPreconditionsTriangulation2_;
        }
        else if (rhs.abstractTriangulation_ == &rhs.periodicPreconditionsTriangulation3_) {
          this->abstractTriangulation_ = &this->periodicPreconditionsTriangulation3_;
        }
        break;
    }
			}
    AbstractTriangulation::operator=(std::move(rhs));

    return *this;
}

Triangulation::~Triangulation() = default;

void Triangulation::switchGrid(const bool usePeriodic,
                               const bool usePreconditions) {
  if(abstractTriangulation_ != nullptr
     && !isImplicitTriangulation()
     && !isHybridImplicitTriangulation()
     && !isPeriodicImplicitTriangulation()
     && !isHybridPeriodicTriangulation()) {
    return;
  }

  if(abstractTriangulation_ != nullptr) {
    // clear preconditions when switching triangulation type
    //if(abstractTriangulation_ == &implicitPreconditionsTriangulation_
      if(isHybridImplicitTriangulation()
       && (usePeriodic || !usePreconditions)) {
      implicitPreconditionsTriangulation0_.clear();
      implicitPreconditionsTriangulation1_.clear();
      implicitPreconditionsTriangulation2_.clear();
      implicitPreconditionsTriangulation3_.clear();
    //} else if(abstractTriangulation_ == &periodicPreconditionsTriangulation_
    } else if (isHybridPeriodicTriangulation()
              && (!usePeriodic || !usePreconditions)) {
      periodicPreconditionsTriangulation0_.clear();
      periodicPreconditionsTriangulation1_.clear();
      periodicPreconditionsTriangulation2_.clear();
      periodicPreconditionsTriangulation3_.clear();
    }
  }

  if(!usePeriodic && !usePreconditions) {
    if (isTriangulation0()) {
      abstractTriangulation_ = &implicitTriangulation0_;
      implicitTriangulation0_.preconditionVerticesAndCells();
    } else if (isTriangulation1()) {
      abstractTriangulation_ = &implicitTriangulation1_;
      implicitTriangulation1_.preconditionVerticesAndCells();
    } else if (isTriangulation2()) {
      abstractTriangulation_ = &implicitTriangulation2_;
      implicitTriangulation2_.preconditionVerticesAndCells();
    } else { //triangulation 3
      abstractTriangulation_ = &implicitTriangulation3_;
      implicitTriangulation3_.preconditionVerticesAndCells();
    }
  } else if(!usePeriodic && usePreconditions) {
    if (isTriangulation0()) {
      abstractTriangulation_ = &implicitPreconditionsTriangulation0_;
      implicitPreconditionsTriangulation0_.preconditionVerticesAndCells();
    } else if (isTriangulation1()) {
      abstractTriangulation_ = &implicitPreconditionsTriangulation1_;
      implicitPreconditionsTriangulation1_.preconditionVerticesAndCells();
    } else if (isTriangulation2()) {
      abstractTriangulation_ = &implicitPreconditionsTriangulation2_;
      implicitPreconditionsTriangulation2_.preconditionVerticesAndCells();
    } else { //triangulation 3
      abstractTriangulation_ = &implicitPreconditionsTriangulation3_;
      implicitPreconditionsTriangulation3_.preconditionVerticesAndCells();
    }
  } else if(usePeriodic && !usePreconditions) {
    if (isTriangulation0()) {
      abstractTriangulation_ = &periodicImplicitTriangulation0_;
      periodicImplicitTriangulation0_.preconditionVerticesAndCells();
    } else if (isTriangulation1()) {
      abstractTriangulation_ = &periodicImplicitTriangulation1_;
      periodicImplicitTriangulation1_.preconditionVerticesAndCells();
    } else if (isTriangulation2()) {
      abstractTriangulation_ = &periodicImplicitTriangulation2_;
      periodicImplicitTriangulation2_.preconditionVerticesAndCells();
    } else { //triangulation 3
      abstractTriangulation_ = &periodicImplicitTriangulation3_;
      periodicImplicitTriangulation3_.preconditionVerticesAndCells();
    }
  } else if(usePeriodic && usePreconditions) {
    if (isTriangulation0()) {
      abstractTriangulation_ = &periodicPreconditionsTriangulation0_;
      periodicPreconditionsTriangulation0_.preconditionVerticesAndCells();
    } else if (isTriangulation1()) {
      abstractTriangulation_ = &periodicPreconditionsTriangulation1_;
      periodicPreconditionsTriangulation1_.preconditionVerticesAndCells();
    } else if (isTriangulation2()) {
      abstractTriangulation_ = &periodicPreconditionsTriangulation2_;
      periodicPreconditionsTriangulation2_.preconditionVerticesAndCells();
    } else { //triangulation 3
      abstractTriangulation_ = &periodicPreconditionsTriangulation3_;
      periodicPreconditionsTriangulation3_.preconditionVerticesAndCells();
    }
  }
}

bool Triangulation::processImplicitStrategy(const STRATEGY strategy) const {

  if(strategy == STRATEGY::DEFAULT) {

#ifndef TTK_IMPLICIT_PRECONDITIONS_THRESHOLD
#define TTK_IMPLICIT_PRECONDITIONS_THRESHOLD 256 * 256 * 256
#endif // TTK_IMPLICIT_PRECONDITIONS_THRESHOLD

    const uint64_t threshold{TTK_IMPLICIT_PRECONDITIONS_THRESHOLD};
    const uint64_t nVerts
      = gridDimensions_[0] * gridDimensions_[1] * gridDimensions_[2];

    // disable preconditioning above TTK_IMPLICIT_PRECONDITIONS_THRESHOLD
    const auto doPreconditioning = nVerts <= threshold;

    if(!doPreconditioning) {

      // ensure that the warning message is printed at creation time
      int prevDebugLevel{-1};
      if(this->debugLevel_ < 2) {
        prevDebugLevel = this->debugLevel_;
        this->debugLevel_ = 2;
      }
      this->printWrn("Large grid detected (> " + std::to_string(threshold)
                     + " vertices)");
      this->printWrn("Defaulting to the fully implicit triangulation");
      if(prevDebugLevel != -1) {
        this->debugLevel_ = prevDebugLevel;
      }
    }
    return doPreconditioning;
  } else if(strategy == STRATEGY::WITH_PRECONDITIONS) {
    return true;
  }
  return false;
}
