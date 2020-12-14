#include "BoundarySystem.h"

namespace hfox{

void BoundarySystem::setBoundaryCondition(FEModel * bModel, const std::set<int> * faces){

  std::tuple< FEModel*, const std::set<int>* > bPair = {bModel, faces};
  boundaryList.push_back(bPair);

};//setBoundaryCondition

};//hfox
