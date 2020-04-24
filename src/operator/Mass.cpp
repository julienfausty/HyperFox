#include "Mass.h"

namespace hfox{

void Mass::assemble(const std::vector< double > & dV, const std::vector< EMatrix > & invJacobians){
  if(!allocated){
    throw(ErrorHandle("Mass", "assemble", "cannot assemble before allocating."));
  }
  const std::vector< std::vector<double> > * shapes = refEl->getIPShapeFunctions();
  int nNodes = refEl->getNumNodes();
  EMatrix phiphi(refEl->getNumIPs(), nNodes*(nNodes+1)/2);
  const std::vector<double> * locShapes;
  int index = 0;
  for(int i = 0; i < refEl->getNumIPs(); i++){
    locShapes = &(shapes->at(i));
    index = 0;
    for(int j = 0; j < nNodes; j++){
      for(int k = j; k < nNodes; k++){
        phiphi(i, index) = locShapes->at(j)*locShapes->at(k);
        index += 1;
      }
    }
  }
  EVector buff = phiphi.transpose() * EMap<const EVector>(dV.data(), dV.size());
  index = 0;
  for(int j = 0; j < nNodes; j++){
    for(int k = j; k < nNodes; k++){
      op(j, k) = buff[index];
      if(j != k){
        op(k, j) = buff[index];
      }
      index += 1;
    }
  }
  if(nDOFsPerNode > 1){
    multiplyDOFs();
  }
};//assemble

}//hfox
