#include "Reaction.h"

namespace hfox{

void Reaction::calcReaction(const std::vector< std::vector<double> > & nodes){
  if(reactionFunc == NULL){
    throw(ErrorHandle("Reaction", "calcReaction", "must set a reaction function before calculating the reaction."));
  }
  reaction.resize(refEl->getNumIPs(), 0.0);
  const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
  std::vector<double> realNode(nodes[0].size());
  EMatrix matNodes(realNode.size(), refEl->getNumNodes());
  for(int i = 0; i < refEl->getNumNodes(); i++){
    matNodes.col(i) = EMap<const EVector>(nodes[i].data(), nodes[i].size());
  }
  for(int i = 0; i < refEl->getNumIPs(); i++){
    EMap<EVector> mapNode(realNode.data(), realNode.size());
    EMap<const EVector> shape(ipShapes->at(i).data(), ipShapes->at(i).size());
    mapNode = matNodes * shape;
    reaction[i] = reactionFunc(realNode);
  }
};//calcReaction

void Reaction::assemble(const std::vector< double > & dV, const std::vector< EMatrix > & invJacobians){
  if(!allocated){
    throw(ErrorHandle("Reaction", "assemble", "the reraction operator must be allocated before beign assembled"));
  }
  if(reaction.size() == 0){
    throw(ErrorHandle("Reaction", "assemble", "the reaction must be calculated in the element before beign assembled"));
  }
  const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
  std::vector<double> locMeasure(refEl->getNumIPs(), 0);
  //pointwise multiplication
  EMap<EVector>(locMeasure.data(), locMeasure.size()) = EMap<const EVector>(reaction.data(), reaction.size()).array()*EMap<const EVector>(dV.data(), reaction.size()).array();
  Mass::assemble(locMeasure, invJacobians);
};//assemble

};
