#include "Source.h"

namespace hfox{

void Source::calcSource(const std::vector< std::vector<double> > & nodes){
  if(sourceFunc == NULL){
    throw(ErrorHandle("Source", "calcSource", "must set a dource function before calculating the source."));
  }
  source.resize(refEl->getNumIPs(), 0.0);
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
    source[i] = sourceFunc(realNode);
  }
};//calcSource

void Source::assemble(const std::vector< double > & dV, const std::vector< EMatrix > & invJacobians){
  if(!allocated){
    throw(ErrorHandle("Source", "assemble", "the source operator must be allocated before beign assembled"));
  }
  if(source.size() == 0){
    throw(ErrorHandle("Source", "assemble", "the source must be calculated in the element before beign assembled"));
  }
  const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
  EVector locMeasure(refEl->getNumIPs());
  //pointwise multiplication
  locMeasure = EMap<const EVector>(source.data(), source.size()).array()*EMap<const EVector>(dV.data(), source.size()).array();
  EMatrix shapes(refEl->getNumNodes(), refEl->getNumIPs());
  for(int i = 0; i < refEl->getNumIPs(); i++){
    shapes.col(i) = EMap<const EVector>(ipShapes->at(i).data(), ipShapes->at(i).size());
  }
  op = shapes * locMeasure;
  if(nDOFsPerNode > 1){
    EVector buff = op.block(0,0,refEl->getNumNodes(), 1);
    for(int i = 0; buff.size(); i++){
      for(int j = 0; j < nDOFsPerNode; j++){
        op(i*nDOFsPerNode + j, 1) = buff(i, 1);
      }
    }
  }
};//assemble

}
