#include "HDGDiffusionSource.h"

namespace hfox{

void HDGDiffusionSource::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator itfm;
  itfm = fm->find("DiffusionTensor");
  if(itfm != fm->end()){
    fieldMap["DiffusionTensor"] = &(itfm->second);
  }
  if(timeScheme != NULL){
    timeScheme->setFieldMap(fm);
  }
  HDGModel::setFieldMap(fm);
};//setFieldMap

void HDGDiffusionSource::allocate(int nDOFsPerNodeUser){
  assembly.matrix = Add;
  assembly.rhs = Add;
  if(timeScheme != NULL){
    timeScheme->allocate(nDOFsPerNodeUser);
  }
  HDGModel::allocate(nDOFsPerNodeUser);
};//allocate

void HDGDiffusionSource::setSourceFunction(std::function<double(const std::vector<double>&)> s){
  if(!allocated){
    throw(ErrorHandle("HDGDiffusionSource", "setSourceFunction", "the model must be allocated before setting the source function"));
  }
  ((Source*)operatorMap["Source"])->setSourceFunction(s); 
};//setSourceFunction

void HDGDiffusionSource::initializeOperators(){
  if(operatorMap.find("Source") == operatorMap.end()){
    operatorMap["Source"] = new Source(refEl);
  }
  if(operatorMap.find("Diffusion") == operatorMap.end()){
    operatorMap["Diffusion"] = new HDGDiffusion(refEl);
  }
  HDGModel::initializeOperators();
};//initializeOperators

void HDGDiffusionSource::computeLocalMatrix(){
  if(!(nodeSet and fieldSet)){
    throw("HDGDiffusionSource", "computeLocalMatrix", "the nodes and the fields should be set before computing");
  }
  computeElementJacobians();
  ((HDGBase*)operatorMap["Base"])->setTau(*(fieldMap["Tau"]));
  ((HDGBase*)operatorMap["Base"])->calcNormals(*elementNodes, jacobians);
  ((HDGDiffusion*)operatorMap["Diffusion"])->setFromBase(((HDGBase*)operatorMap["Base"])->getNormals());
  if(fieldMap.find("DiffusionTensor") != fieldMap.end()){
    ((HDGDiffusion*)operatorMap["Diffusion"])->setDiffusionTensor(parseDiffusionVals());
  }
  operatorMap["Base"]->assemble(dV, invJacobians);
  operatorMap["Diffusion"]->assemble(dV, invJacobians);
  localMatrix = *(operatorMap["Base"]->getMatrix());
  localMatrix += *(operatorMap["Diffusion"]->getMatrix());
};//computeLocalMatrix

std::vector<EMatrix> HDGDiffusionSource::parseDiffusionVals() const{
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
    throw(ErrorHandle("HDGDiffusionSource", "parseDiffusionVals", "the dimension of the diffusion tensor vales are not correct, they should be either scalar or tensor of the dimension of the reference element"));
  }
  return res;
};//parseDiffusionVals

void HDGDiffusionSource::computeLocalRHS(){
  ((Source*)operatorMap["Source"])->calcSource(*elementNodes);
  operatorMap["Source"]->assemble(dV, invJacobians);
  localRHS.segment(0, operatorMap["Source"]->getMatrix()->rows()) = (EVector) *(operatorMap["Source"]->getMatrix());
};//computeLocalRHS

}//hfox
