#include "BoundarySystem.h"

namespace hfox{

void BoundarySystem::setBoundaryCondition(BoundaryModel * bModel, const std::set<int> * faces){

  std::tuple< BoundaryModel*, const std::set<int>* > bPair {bModel, faces};
  boundaryList.push_back(bPair);

};//setBoundaryCondition

};//hfox
