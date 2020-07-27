#include "HDGConvectionDiffusionReactionSource.h"

namespace hfox{

void HDGConvectionDiffusionReactionSource::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  bool velFound = 0; 
  bool diffFound = 0;
  std::map<std::string, std::vector<double> >::const_iterator itfm;
  itfm = fm->find("DiffusionTensor");
  if(itfm != fm->end()){
    fieldMap["DiffusionTensor"] = &(itfm->second);
    diffFound = 1;
  }
  itfm = fm->find("Velocity");
  if(itfm != fm->end()){
    fieldMap["Velocity"] = &(itfm->second);
    velFound = 1;
  }
  if(!velFound and !diffFound){
    throw(ErrorHandle("HDGConvectionDiffusionReactionSource", "setFieldMap", "must provide at least either a Velocity field or a DiffusionTensor field"));
  }
  if(timeScheme != NULL){
    timeScheme->setFieldMap(fm);
  }
  HDGModel::setFieldMap(fm);
};//setFieldMap

void HDGConvectionDiffusionReactionSource::allocate(int nDOFsPerNodeUser){
  assembly.matrix = Add;
  assembly.rhs = Add;
  if(timeScheme != NULL){
    timeScheme->allocate(nDOFsPerNodeUser);
  }
  HDGModel::allocate(nDOFsPerNodeUser);
};//allocate

void HDGConvectionDiffusionReactionSource::setSourceFunction(std::function<double(const std::vector<double>&)> s){
  if(!allocated){
    throw(ErrorHandle("HDGConvectionDiffusionReactionSource", "setSourceFunction", "the model must be allocated before setting the source function"));
  }
  ((Source*)operatorMap["Source"])->setSourceFunction(s);
  sourceSet = 1;
};//setSourceFunction

void HDGConvectionDiffusionReactionSource::setReactionFunction(std::function<double(const std::vector<double>&)> r){
  if(!allocated){
    throw(ErrorHandle("HDGConvectionDiffusionReactionSource", "setReactionFunction", "the model must be allocated before setting the reaction function"));
  }
  ((Reaction*)operatorMap["Reaction"])->setReactionFunction(r);
  reactionSet = 1;
};//setSourceFunction

void HDGConvectionDiffusionReactionSource::initializeOperators(){
  if(operatorMap.find("Source") == operatorMap.end()){
    operatorMap["Source"] = new Source(refEl);
  }
  if(operatorMap.find("Convection") == operatorMap.end()){
    operatorMap["Convection"] = new HDGConvection(refEl);
  }
  if(operatorMap.find("Diffusion") == operatorMap.end()){
    operatorMap["Diffusion"] = new HDGDiffusion(refEl);
  }
  if(operatorMap.find("Reaction") == operatorMap.end()){
    operatorMap["Reaction"] = new Reaction(refEl);
  }
  HDGModel::initializeOperators();
};//initializeOperators

void HDGConvectionDiffusionReactionSource::computeLocalMatrix(){
  if(!(nodeSet and fieldSet)){
    throw(ErrorHandle("HDGConvectionDiffusionReactionSource", "computeLocalMatrix", "the nodes and the fields should be set before computing"));
  }
  computeElementJacobians();
  ((HDGBase*)operatorMap["Base"])->setTau(*(fieldMap["Tau"]));
  ((HDGBase*)operatorMap["Base"])->calcNormals(*elementNodes, jacobians);
  operatorMap["Base"]->assemble(dV, invJacobians);
  localMatrix = *(operatorMap["Base"]->getMatrix());
  if(fieldMap.find("Velocity") != fieldMap.end()){
    ((HDGConvection*)operatorMap["Convection"])->setVelocity(parseVelocityVals());
    ((HDGConvection*)operatorMap["Convection"])->setFromBase(((HDGBase*)operatorMap["Base"])->getNormals());
    operatorMap["Convection"]->assemble(dV, invJacobians);
    localMatrix += *(operatorMap["Convection"]->getMatrix());
  }
  if(fieldMap.find("DiffusionTensor") != fieldMap.end()){
    ((HDGDiffusion*)operatorMap["Diffusion"])->setDiffusionTensor(parseDiffusionVals());
    ((HDGDiffusion*)operatorMap["Diffusion"])->setFromBase(((HDGBase*)operatorMap["Base"])->getNormals());
    operatorMap["Diffusion"]->assemble(dV, invJacobians);
    localMatrix += *(operatorMap["Diffusion"]->getMatrix());
  }
  if(reactionSet){
    ((Reaction*)operatorMap["Reaction"])->calcReaction(*elementNodes);
    operatorMap["Reaction"]->assemble(dV, invJacobians);
    const EMatrix * rM = operatorMap["Reaction"]->getMatrix();
    localMatrix.block(0, 0, rM->rows(), rM->cols()) += *(rM);
  }
};//computeLocalMatrix

void HDGConvectionDiffusionReactionSource::computeLocalRHS(){
  if(!(nodeSet and fieldSet)){
    throw("HDGConvectionDiffusionReactionSource", "computeLocalRHS", "the nodes and the fields should be set before computing");
  }
  localRHS = EVector::Zero(localRHS.size());
  if(sourceSet){
    ((Source*)operatorMap["Source"])->calcSource(*elementNodes);
    operatorMap["Source"]->assemble(dV, invJacobians);
    localRHS.segment(0, operatorMap["Source"]->getMatrix()->rows()) = (EVector) *(operatorMap["Source"]->getMatrix());
  }
};//computeLocalRHS


std::vector<EMatrix> HDGConvectionDiffusionReactionSource::parseDiffusionVals() const{
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
    throw(ErrorHandle("HDGConvectionDiffusionReactionSource", "parseDiffusionVals", "the dimension of the diffusion tensor vales are not correct, they should be either scalar or tensor of the dimension of the reference element"));
  }
  return res;
};//parseDiffusionVals

std::vector<EVector> HDGConvectionDiffusionReactionSource::parseVelocityVals() const{
  const std::vector<double> * velVals = fieldMap.at("Velocity");
  int nNodes = refEl->getNumNodes();
  int dimVelVal = velVals->size()/nNodes;
  int dimSpace = refEl->getDimension();
  if(dimVelVal != dimSpace){
    throw(ErrorHandle("HDGConvectionDiffusionReactionSource", "parseVelocityVals", "the dimension of the velocity vector does not correspond to the dimension of the reference element"));
  }
  std::vector<EVector> res(nNodes, EVector::Zero(dimSpace));
  for(int i = 0; i < nNodes; i++){
    res[i] = EMap<const EVector>(velVals->data() + i*dimSpace, dimSpace);
  }
  return res;
};//parseVelocityVals

};
