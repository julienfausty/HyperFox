#include "DiffusionSource.h"

namespace hfox{

void DiffusionSource::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator itfm;
  itfm = fm->find("DiffusionTensor");
  if(itfm != fm->end()){
    fieldMap["DiffusionTensor"] = &(itfm->second);
  }
  if(timeScheme != NULL){
    timeScheme->setFieldMap(fm);
  }
  fieldSet = 1;
};//setFieldMap

void DiffusionSource::allocate(int nDOFsPerNode){
  assembly.matrix = Add;
  assembly.rhs = Add;
  if(timeScheme != NULL){
    timeScheme->allocate(nDOFsPerNode);
  }
  initializeOperators();
  int n = nDOFsPerNode * (refEl->getNumNodes());
  localMatrix = EMatrix::Zero(n, n);
  localRHS = EVector::Zero(n);
  for(std::map<std::string, Operator*>::iterator it = operatorMap.begin(); it != operatorMap.end(); it++){
    (it->second)->allocate(nDOFsPerNode);
  }
  allocated = 1;
};//allocate

void DiffusionSource::setSourceFunction(std::function<double(const std::vector<double>&)> s){
  if(!allocated){
    throw(ErrorHandle("DiffusionSource", "setSourceFunction", "the model must be allocated before setting the source function"));
  }
  ((Source*)operatorMap["Source"])->setSourceFunction(s); 
};//setSourceFunction

void DiffusionSource::initializeOperators(){
  if(operatorMap.find("Source") == operatorMap.end()){
    operatorMap["Source"] = new Source(refEl);
  }
  if(operatorMap.find("Diffusion") == operatorMap.end()){
    operatorMap["Diffusion"] = new Diffusion(refEl);
  }
};//initializeOperators

std::vector<EMatrix> DiffusionSource::parseDiffusionVals() const{
  const std::vector<double> * diffVals = fieldMap.at("DiffusionTensor");
  int nNodes = refEl->getNumNodes();
  int dimDiffVal = diffVals->size()/nNodes;
  int dimMat = std::sqrt(dimDiffVal);
  int dimSpace = refEl->getDimension();
  std::vector<EMatrix> res(nNodes, EMatrix::Identity(dimSpace, dimSpace));
  if(dimMat == 1){
    for(int i = 0; i < nNodes; i++){
      res[i] *= (diffVals->at(i));
    }
  } else if(dimMat == dimSpace){
    for(int i = 0; i < nNodes; i++){
      res[i] = EMap<const EMatrix>(diffVals->data() + i*dimDiffVal, dimMat, dimMat);
    }
  } else{
    throw(ErrorHandle("DiffusionSource", "parseDiffusionVals", "the dimension of the diffusion tensor vales are not correct, they should be either scalar or tensor of the dimension of the reference element"));
  }
  return res;
};//parseDiffusionVals

void DiffusionSource::computeLocalMatrix(){
  if(!(nodeSet and fieldSet)){
    throw("DiffusionSource", "computeLocalMatrix", "the nodes and the fields should be set before computing");
  }
  jacobians = Operator::calcJacobians(*elementNodes, refEl);
  invJacobians = Operator::calcInvJacobians(jacobians);
  dV = Operator::calcMeasure(Operator::calcDetJacobians(jacobians), refEl);
  if(fieldMap.find("DiffusionTensor") != fieldMap.end()){
    ((Diffusion*)operatorMap["Diffusion"])->setDiffTensor(parseDiffusionVals());
  }
  operatorMap["Diffusion"]->assemble(dV, invJacobians);
  localMatrix = *(operatorMap["Diffusion"]->getMatrix());
};//computeLocalMatrix

void DiffusionSource::computeLocalRHS(){
  if(!(nodeSet and fieldSet)){
    throw("DiffusionSource", "computeLocalRHS", "the nodes and the fields should be set before computing");
  }
  ((Source*)operatorMap["Source"])->calcSource(*elementNodes);
  operatorMap["Source"]->assemble(dV, invJacobians);
  localRHS = (EVector) *(operatorMap["Source"]->getMatrix());
};//computeLocalRHS


};
