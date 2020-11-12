#include "HDGnGammaModel.h"

namespace hfox{

namespace nGamma{

void HDGnGammaModel::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator itfm;
  itfm = fm->find("D");
  if(itfm != fm->end()){
    fieldMap["D"] = &(itfm->second);
  }
  itfm = fm->find("G");
  if(itfm != fm->end()){
    fieldMap["G"] = &(itfm->second);
  }
  itfm = fm->find("b");
  if(itfm != fm->end()){
    fieldMap["b"] = &(itfm->second);
  } else {
    throw(ErrorHandle("HDGnGammaModel", "setFieldMap", "must provide a b vector field to the nGamma model"));
  }
  itfm = fm->find("Solution");
  if(itfm != fm->end()){
    fieldMap["Solution"] = &(itfm->second);
  } else {
    throw(ErrorHandle("HDGnGammaModel", "setFieldMap", "must provide a Solution field for the Newton-Raphson iterations"));
  }
  if(itfm->second.size() != 2*(refEl->getNumNodes())){
    throw(ErrorHandle("HDGnGammaModel", "setFieldMap", "the Solution must have 2 DOFs per node n and Gamma"));
  }
  itfm = fm->find("Trace");
  if(itfm != fm->end()){
    fieldMap["Trace"] = &(itfm->second);
  } else {
    throw(ErrorHandle("HDGnGammaModel", "setFieldMap", "must provide a Trace field for the Newton-Raphson iterations"));
  }
  if(itfm->second.size() != 2*(refEl->getFaceElement()->getNumNodes())*(refEl->getNumFaces())){
    throw(ErrorHandle("HDGnGammaModel", "setFieldMap", "the Trace must have 2 DOFs per node"));
  }
  if(timeScheme != NULL){
    timeScheme->setFieldMap(fm);
  }
  HDGModel::setFieldMap(fm);
};//setFieldMap

void HDGnGammaModel::allocate(int nDOFsPerNodeUser){
  if(nDOFsPerNodeUser != 2){
    throw(ErrorHandle("HDGnGammaModel", "allocate", "the number of DOFs per node must be equal to 2: n and Gamma"));
  }
  assembly.matrix = Add;
  assembly.rhs = Add;
  if(timeScheme != NULL){
    timeScheme->allocate(nDOFsPerNodeUser);
  }
  initializeOperators();
  nDOFsPNode = nDOFsPerNodeUser;
  int n = nDOFsPNode * (refEl->getNumNodes()*(refEl->getDimension() + 1) 
      + (refEl->getFaceElement()->getNumNodes())*(refEl->getNumFaces()));
  localMatrix = EMatrix::Zero(n, n);
  localRHS = EVector::Zero(n);
  for(std::map<std::string, Operator*>::iterator it = operatorMap.begin(); it != operatorMap.end(); it++){
    if(it->first.find("Convection") != std::string::npos){
      (it->second)->allocate(1);
    } else if(it->first.find("Diffusion") != std::string::npos){
      (it->second)->allocate(1);
    } else if(it->first.find("Source") != std::string::npos){
      (it->second)->allocate(1);
    } else {
      (it->second)->allocate(nDOFsPNode);
    }
  }
  allocated = 1;
};//allocate

void HDGnGammaModel::initializeOperators(){
  for(int i = 0; i < 4; i++){
    if(operatorMap.find("Convection_" + std::to_string(i)) == operatorMap.end()){
      operatorMap["Convection_" + std::to_string(i)] = new HDGConvection(refEl);
    }
  }
  for(int i = 0; i < 2; i++){
    if(operatorMap.find("Diffusion_" + std::to_string(i)) == operatorMap.end()){
      operatorMap["Diffusion_" + std::to_string(i)] = new HDGDiffusion(refEl);
    }
  }
  for(int i = 0; i < 2; i++){
    if(operatorMap.find("Source_" + std::to_string(i)) == operatorMap.end()){
      operatorMap["Source_" + std::to_string(i)] = new Source(refEl);
    }
  }
  HDGModel::initializeOperators();
};//initializeOperators

void HDGnGammaModel::setSourceFunction(std::function<double(const std::vector<double>&, int)> s){
  if(!allocated){
    throw(ErrorHandle("HDGnGammaModel", "setSourceFunction", "the model must be allocated before setting the source function"));
  }
  for(int i = 0; i < refEl->getDimension(); i++){
    ((Source*)operatorMap["Source_" + std::to_string(i)])->setSourceFunction([i,s](const std::vector<double> & x){return s(x, i);});
  }
  sourceSet = 1;
};//setSourceFunction

void HDGnGammaModel::computeLocalMatrix(){
  if(!(nodeSet and fieldSet)){
    throw(ErrorHandle("HDGnGammaModel", "computeLocalMatrix", "the nodes and the fields should be set before computing"));
  }
  computeElementJacobians();
  ((HDGBase*)operatorMap["Base"])->setTau(*(fieldMap["Tau"]));
  ((HDGBase*)operatorMap["Base"])->calcNormals(*elementNodes, jacobians);
  operatorMap["Base"]->assemble(dV, invJacobians);
  localMatrix = *(operatorMap["Base"]->getMatrix());
  std::vector<EVector> bs = parseBVals();
  std::vector<EVector> vs(refEl->getNumNodes(), EVector::Zero(2));
};//computeLocalMatrix

void HDGnGammaModel::computeLocalRHS(){
  if(!(nodeSet and fieldSet)){
    throw(ErrorHandle("HDGnGammaModel", "computeLocalRHS", "the nodes and the fields should be set before computing"));
  }
  localRHS = EVector::Zero(localRHS.size());
};//computeLocalRHS

std::vector<EVector> HDGnGammaModel::parseBVals(){

};//parseBVals

std::vector<EMatrix> HDGnGammaModel::parseDVals(){

};//parseDVals


std::vector<EMatrix> HDGnGammaModel::parseGVals(){

};//parseGVals

std::vector<EVector> HDGnGammaModel::parseSolutionVals(){

};//parseSolutionVals

};//nGamma

};//hfox
