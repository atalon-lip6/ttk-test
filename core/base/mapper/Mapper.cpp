#include <Mapper.h>

ttk::Mapper::Mapper() {
  this->setDebugMsgPrefix("Mapper");
}

int ttk::Mapper::updateNonCentroidsCoords(
  float *const outputPointsCoords, const std::vector<float> &pointsPrev) {
  const size_t nbPoint = compSpecialCoeffToSave_.size() / 3;
  if(pointsPrev.size() != 3 * nbPoint) {
    printErr("Error in updating the dilatation coefficient only. We want to "
             "update the coordinates of "
             + std::to_string(nbPoint)
             + " points but the data we saved to do so concerns "
             + std::to_string(compSpecialCoeffToSave_.size()) + " points.");
    return 1;
  }

  size_t dim = LowerDimension == LOWER_DIMENSION::LOWER_DIM_2D ? 2 : 3;

  for(size_t i = 0; i < nbPoint; i++) {
    if(i < 5) {
      std::cerr << i << " pPrev => " << pointsPrev[3 * i] << ","
                << pointsPrev[3 * i + 1] << "," << pointsPrev[3 * i + 2]
                << std::endl;
      std::cerr << " oh la qui voilà = " << compSpecialCoeffToSave_[3 * i + 0]
                << " - " << compSpecialCoeffToSave_[3 * i + 1] << std::endl;
    }
    for(size_t iDim = 0; iDim < dim; iDim++) {
      outputPointsCoords[3 * i + iDim]
        = pointsPrev[3 * i + iDim]
          + (DilatationCoeff - prevDilatationCoeff_)
              * compSpecialCoeffToSave_[3 * i + iDim];
    }
    if(i < 5)
      std::cerr << i << " => " << outputPointsCoords[3 * i] << ","
                << outputPointsCoords[3 * i + 1] << ","
                << outputPointsCoords[3 * i + 2] << std::endl;
  }
  prevDilatationCoeff_ = DilatationCoeff;

  return 0;
}

int ttk::Mapper::reduceMatrix(
  std::vector<std::vector<double>> &outputCoords,
  const Matrix &mat,
  const bool isDistanceMatrix,
  const ttk::DimensionReduction::METHOD method) const {

  if(isDistanceMatrix) {
    if(this->ReductionAlgo == REDUCTION_ALGO::LLE
       || this->ReductionAlgo == REDUCTION_ALGO::PCA) {
      this->printErr("LLE and PCA don't support distances matrices as input");
      return -1;
    }
  }

  ttk::DimensionReduction dimRed{};
  dimRed.setDebugLevel(0);
  dimRed.setThreadNumber(1);

  if(!dimRed.isPythonFound()) {
    this->printErr(
      "Missing Python modules, could not perform Dimension Reduction");
    return -2;
  }

  dimRed.setInputMethod(method);
  dimRed.setIsInputDistanceMatrix(isDistanceMatrix);
  dimRed.setInputNumberOfComponents(static_cast<int>(this->LowerDimension));

  const auto ret
    = dimRed.execute(outputCoords, nullptr, mat.data(), mat.nRows(), mat.nCols());
  if(ret != 0) {
    this->printErr("DimensionReduction returned error code "
                   + std::to_string(ret));
    return ret;
  }

  return 0;
}

void ttk::Mapper::computeCompCentroid(
  std::vector<SimplexId> &centroidId,
  const Matrix &distMat,
  const std::vector<std::set<SimplexId>> &connCompVertices,
  const int *const outputConnComp) const {

  centroidId.resize(connCompVertices.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < connCompVertices.size(); ++i) {
    const auto &cc{connCompVertices[i]};
    bool useExtendedVerts{true};
    for(const auto v : cc) {
      if(outputConnComp[v] == static_cast<int>(i)) {
        // connected component contains at least one vertex
        useExtendedVerts = false;
        break;
      }
    }

    // for each connected component, find the vertex that minimize the
    // square distance to every other component vertex
    std::pair<SimplexId, double> min{-1, std::numeric_limits<double>::max()};
    for(const auto v : cc) {
      if(!useExtendedVerts && outputConnComp[v] != static_cast<int>(i)) {
        continue;
      }

      double res{};
      for(const auto w : cc) {
        if(v == w) {
          continue;
        }
        res += distMat.get(v, w) * distMat.get(v, w);
      }
      if(res < min.second) {
        min = {v, res};
      }
    }
    centroidId[i] = min.first;
  }
}

void ttk::Mapper::extractSubDistMat(Matrix &subDistMat,
                                    const std::vector<SimplexId> &vertsId,
                                    const Matrix &distMat) const {

  const auto nConnComps{vertsId.size()};
  subDistMat.alloc(nConnComps, nConnComps);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nConnComps - 1; ++i) {
    for(size_t j = i + 1; j < nConnComps; ++j) {
      const auto dij = distMat.get(vertsId[i], vertsId[j]);
      subDistMat.get(i, j) = dij;
      subDistMat.get(j, i) = dij;
    }
  }
}
