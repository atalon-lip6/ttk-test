#include <MultiresTopology.h>

void ttk::MultiresTopology::getValencesFromLink(
  const SimplexId vertexId,
  const std::vector<std::pair<polarity, polarity>> &vlp,
  DynamicTree &link,
  std::vector<polarity> &toPropagateMin,
  std::vector<polarity> &toPropagateMax,
  std::vector<std::vector<SimplexId>> &saddleCCMin,
  std::vector<std::vector<SimplexId>> &saddleCCMax) const {

  const auto nbCC = link.getNbCC();

  SimplexId downValence = 0, upValence = 0;
  saddleCCMin[vertexId].clear();
  saddleCCMax[vertexId].clear();

  if(nbCC > 2) {
    std::vector<size_t> CCIds;
    CCIds.reserve(nbCC);
    link.retrieveNbCC(CCIds);
    for(size_t i = 0; i < CCIds.size(); i++) {
      const SimplexId neighbor = CCIds[i];
      const polarity isUpper = vlp[neighbor].first;
      if(isUpper) {
        saddleCCMax[vertexId].emplace_back(neighbor);
        upValence++;
      } else {
        saddleCCMin[vertexId].emplace_back(neighbor);
        downValence++;
      }
    }

    if(downValence > 1) {
      toPropagateMin[vertexId] = 255;
    } else {
      saddleCCMin[vertexId].clear();
      toPropagateMin[vertexId] = 0;
    }
    if(upValence > 1) {
      toPropagateMax[vertexId] = 255;
    } else {
      saddleCCMax[vertexId].clear();
      toPropagateMax[vertexId] = 0;
    }

  } else { // not a saddle
    toPropagateMax[vertexId] = 0;
    toPropagateMin[vertexId] = 0;
  }
}




template<> void ttk::MultiresTopology::buildVertexLinkByBoundary<0>(
  const SimplexId vertexId, VLBoundaryType &vlbt) const {

  const auto bid = multiresTriangulation0_.getVertexBoundaryIndex(vertexId);
  const auto nneigh = multiresTriangulation0_.getVertexNeighborNumber(vertexId);
  vlbt[bid].reserve(nneigh);

  for(SimplexId i = 0; i < nneigh; i++) {
    SimplexId n0 = 0;
    multiresTriangulation0_.getVertexNeighbor(vertexId, i, n0);
    for(SimplexId j = i + 1; j < nneigh; j++) {
      SimplexId n1 = 0;
      multiresTriangulation0_.getVertexNeighbor(vertexId, j, n1);
      if(multiresTriangulation0_.areVerticesNeighbors(n0, n1)) {
        vlbt[bid].emplace_back(i, j);
      }
    }
  }
}

template<> void ttk::MultiresTopology::buildVertexLinkByBoundary<1>(
  const SimplexId vertexId, VLBoundaryType &vlbt) const {

  const auto bid = multiresTriangulation1_.getVertexBoundaryIndex(vertexId);
  const auto nneigh = multiresTriangulation1_.getVertexNeighborNumber(vertexId);
  vlbt[bid].reserve(nneigh);

  for(SimplexId i = 0; i < nneigh; i++) {
    SimplexId n0 = 0;
    multiresTriangulation1_.getVertexNeighbor(vertexId, i, n0);
    for(SimplexId j = i + 1; j < nneigh; j++) {
      SimplexId n1 = 0;
      multiresTriangulation1_.getVertexNeighbor(vertexId, j, n1);
      if(multiresTriangulation1_.areVerticesNeighbors(n0, n1)) {
        vlbt[bid].emplace_back(i, j);
      }
    }
  }
}

template<> void ttk::MultiresTopology::buildVertexLinkByBoundary<2>(
  const SimplexId vertexId, VLBoundaryType &vlbt) const {

  const auto bid = multiresTriangulation2_.getVertexBoundaryIndex(vertexId);
  const auto nneigh = multiresTriangulation2_.getVertexNeighborNumber(vertexId);
  vlbt[bid].reserve(nneigh);

  for(SimplexId i = 0; i < nneigh; i++) {
    SimplexId n0 = 0;
    multiresTriangulation2_.getVertexNeighbor(vertexId, i, n0);
    for(SimplexId j = i + 1; j < nneigh; j++) {
      SimplexId n1 = 0;
      multiresTriangulation2_.getVertexNeighbor(vertexId, j, n1);
      if(multiresTriangulation2_.areVerticesNeighbors(n0, n1)) {
        vlbt[bid].emplace_back(i, j);
      }
    }
  }
}

template<> void ttk::MultiresTopology::buildVertexLinkByBoundary<3>(
  const SimplexId vertexId, VLBoundaryType &vlbt) const {

  const auto bid = multiresTriangulation3_.getVertexBoundaryIndex(vertexId);
  const auto nneigh = multiresTriangulation3_.getVertexNeighborNumber(vertexId);
  vlbt[bid].reserve(nneigh);

  for(SimplexId i = 0; i < nneigh; i++) {
    SimplexId n0 = 0;
    multiresTriangulation3_.getVertexNeighbor(vertexId, i, n0);
    for(SimplexId j = i + 1; j < nneigh; j++) {
      SimplexId n1 = 0;
      multiresTriangulation3_.getVertexNeighbor(vertexId, j, n1);
      if(multiresTriangulation3_.areVerticesNeighbors(n0, n1)) {
        vlbt[bid].emplace_back(i, j);
      }
    }
  }
}

void ttk::MultiresTopology::buildVertexLinkByBoundary(
  const SimplexId vertexId, VLBoundaryType &vlbt) const {
  if (dimensionality_ == 0) {
    buildVertexLinkByBoundary<0>(vertexId, vlbt);
  }
  else if (dimensionality_ == 1) {
    buildVertexLinkByBoundary<1>(vertexId, vlbt);
  }
  else if (dimensionality_ == 2) {
    buildVertexLinkByBoundary<2>(vertexId, vlbt);
  }
  else if (dimensionality_ == 3) {
    buildVertexLinkByBoundary<3>(vertexId, vlbt);
  }
}



void ttk::MultiresTopology::getTripletsFromSaddles(
  const SimplexId vertexId,
  std::vector<triplet> &triplets,
  const std::vector<std::vector<SimplexId>> &vertexReps) const {

  const auto &reps = vertexReps[vertexId];
  const SimplexId m = reps[0];
#ifndef TTK_ENABLE_KAMIKAZE
  const auto &repsm = vertexReps[m];
  if(m == -1 || repsm.empty() || repsm[0] != m) {
    this->printErr("HERE PROBLEM");
  }
#endif // TTK_ENABLE_KAMIKAZE
  for(size_t i = 1; i < reps.size(); i++) {
    const SimplexId n = reps[i];
#ifndef TTK_ENABLE_KAMIKAZE
    const auto &repsn = vertexReps[n];
    if(n == -1 || repsn.empty() || repsn[0] != n) {
      this->printErr("HERE2 PROBLEM");
    }
#endif // TTK_ENABLE_KAMIKAZE
    triplets.emplace_back(vertexId, m, n);
  }
}

char ttk::MultiresTopology::getCriticalTypeFromLink(
  const std::vector<std::pair<polarity, polarity>> &vlp,
  DynamicTree &link) const {

  const auto nbCC = link.getNbCC();

  int const dimensionality = this->dimensionality_;//multiresTriangulation_.getDimensionality();
  SimplexId downValence = 0, upValence = 0;

  std::vector<size_t> CCIds;
  CCIds.reserve(nbCC);
  link.retrieveNbCC(CCIds);
  for(size_t i = 0; i < CCIds.size(); i++) {
    const SimplexId neighbor = CCIds[i];
    const polarity isUpper = vlp[neighbor].first;
    if(isUpper) {
      upValence++;
    } else {
      downValence++;
    }
  }

  if(downValence == -1 && upValence == -1) {
    return -1;
  } else if(downValence == 0 && upValence == 1) {
    return static_cast<char>(CriticalType::Local_minimum);
  } else if(downValence == 1 && upValence == 0) {
    return static_cast<char>(CriticalType::Local_maximum);
  } else if(downValence == 1 && upValence == 1) {
    // regular point
    return static_cast<char>(CriticalType::Regular);
  } else {
    // saddles
    if(dimensionality == 2) {
      if((downValence == 2 && upValence == 1)
         || (downValence == 1 && upValence == 2)
         || (downValence == 2 && upValence == 2)) {
        // regular saddle
        return static_cast<char>(CriticalType::Saddle1);
      } else {
        // monkey saddle, saddle + extremum
        return static_cast<char>(CriticalType::Degenerate);
        // NOTE: you may have multi-saddles on the boundary in that
        // configuration
        // to make this computation 100% correct, one would need to
        // disambiguate boundary from interior vertices
      }
    } else if(dimensionality == 3) {
      if(downValence == 2 && upValence == 1) {
        return static_cast<char>(CriticalType::Saddle1);
      } else if(downValence == 1 && upValence == 2) {
        return static_cast<char>(CriticalType::Saddle2);
      } else {
        // monkey saddle, saddle + extremum
        return static_cast<char>(CriticalType::Degenerate);
        // NOTE: we may have a similar effect in 3D (TODO)
      }
    }
  }

  // -2: regular points
  return static_cast<char>(CriticalType::Regular);
}

std::string ttk::MultiresTopology::resolutionInfoString() {
  std::stringstream res;
  if (dimensionality_ == 0) {
    res << "Resolution level "
      << multiresTriangulation0_.DL_to_RL(decimationLevel_);
  }
  else if (dimensionality_ == 1) {
    res << "Resolution level "
      << multiresTriangulation1_.DL_to_RL(decimationLevel_);
  }
  else if (dimensionality_ == 2) {
    res << "Resolution level "
      << multiresTriangulation2_.DL_to_RL(decimationLevel_);
  }
  else if (dimensionality_ == 3) {
    res << "Resolution level "
      << multiresTriangulation3_.DL_to_RL(decimationLevel_);
  }
  else {
    this->printErr("Error: dimensionality should be between 0 and 3 inclusive.");
    return "";
  }
  if(decimationLevel_ == 0) {
    res << " (final)";
  }
  return res.str();
}

void ttk::MultiresTopology::updateLinkPolarityPonctual(
  std::vector<std::pair<polarity, polarity>> &vlp) const {

  for(size_t i = 0; i < vlp.size(); i++) {
    if(vlp[i].second) {
      vlp[i].first = ~vlp[i].first;
      vlp[i].second = 0;
    }
  }
}
