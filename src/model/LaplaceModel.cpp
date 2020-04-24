#include "LaplaceModel.h"

namespace hfox{

void LaplaceModel::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator it;
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "DiffusionTensor"){
      fieldMap[it->first] = &(it->second);
    }
  }
  fieldSet = 1;
};//setFieldMap


void LaplaceModel::allocate(int nDOFsPerNode){
  assembly.matrix = Add;
  assembly.rhs = None;
  initializeOperators();
  localMatrix = EMatrix::Zero(nDOFsPerNode*refEl->getNumNodes(), nDOFsPerNode*refEl->getNumNodes());
  localRHS = EVector::Zero(nDOFsPerNode*refEl->getNumNodes());
  for(std::map<std::string, Operator*>::iterator it = operatorMap.begin(); it != operatorMap.end(); it++){
    (it->second)->allocate(nDOFsPerNode);
  }
  allocated = 1;
};//allocate

void LaplaceModel::initializeOperators(){
  operatorMap["Diffusion"] = new Diffusion(refEl);
};//initializeOperators

void LaplaceModel::computeLocalMatrix(){
  if(nodeSet and allocated){
    jacobians = Operator::calcJacobians(*elementNodes, refEl);
    invJacobians = Operator::calcInvJacobians(jacobians);
    dV = Operator::calcMeasure(Operator::calcDetJacobians(jacobians), refEl);
    std::map<std::string, const std::vector<double> * >::iterator it = fieldMap.find("DiffusionTensor");
    if(it != fieldMap.end()){
      int dimTensor = std::sqrt(it->second->size()/(refEl->getNumNodes()));
      std::vector< EMatrix > diffTensors(refEl->getNumNodes(), EMatrix::Zero(dimTensor, dimTensor));
      int s = dimTensor*dimTensor;
      for(int i = 0; i < refEl->getNumNodes(); i++){
        diffTensors[i] = Eigen::Map<const EMatrix>(it->second->data() + i*s, dimTensor, dimTensor);
      }
      ((Diffusion*)operatorMap["Diffusion"])->setDiffTensor(diffTensors);
    }
    for(std::map<std::string, Operator*>::iterator it = operatorMap.begin(); it != operatorMap.end(); it++){
      (it->second)->assemble(dV, invJacobians);
    }
    localMatrix = *(operatorMap["Diffusion"]->getMatrix());
  } else {
    throw(ErrorHandle("LaplaceModel", "computeLocalMatrix", 
          "the nodes should be set and the object allocated before computing."));
  }
};//computeLocalMatrix

void LaplaceModel::computeLocalRHS(){
  if(!(nodeSet and allocated)){
    throw(ErrorHandle("LaplaceModel", "computeLocalRHS", 
          "the nodes should be set and the object allocated before computing."));
  }
};//computeLocalRHS


}//hfox
