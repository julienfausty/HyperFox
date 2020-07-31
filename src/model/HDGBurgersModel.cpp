#include "HDGBurgersModel.h"

namespace hfox{

void HDGBurgersModel::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator itfm;
  itfm = fm->find("DiffusionTensor");
  if(itfm != fm->end()){
    fieldMap["DiffusionTensor"] = &(itfm->second);
  }
  itfm = fm->find("BufferSolution");
  if(itfm != fm->end()){
    fieldMap["BufferSolution"] = &(itfm->second);
  } else {
    throw(ErrorHandle("HDGBurgersModel", "setFieldMap", "must provide a BufferSolution field for the Newton-Raphson iterations"));
  }
  if(itfm->second.size() != refEl->getDimension()*(refEl->getNumNodes())){
    throw(ErrorHandle("HDGBurgersModel", "setFieldMap", "the BufferSolution must have the same number of DOFs per node as spatial dimensions for the Burgers equation"));
  }
  if(timeScheme != NULL){
    timeScheme->setFieldMap(fm);
  }
  HDGModel::setFieldMap(fm);
};//setFieldMap

void HDGBurgersModel::allocate(int nDOFsPerNodeUser){
  if(nDOFsPerNodeUser != refEl->getDimension()){
    throw(ErrorHandle("HDGBurgersModel", "allocate", "the number of DOFs per node must be equal to the number of spatial dimensions for the Burgers equation"));
  }
  assembly.matrix = Add;
  assembly.rhs = Add;
  if(timeScheme != NULL){
    timeScheme->allocate(nDOFsPerNodeUser);
  }
  HDGModel::allocate(nDOFsPerNodeUser);
};//allocate


void HDGBurgersModel::initializeOperators(){
  if(operatorMap.find("Convection") == operatorMap.end()){
    operatorMap["Convection"] = new HDGConvection(refEl);
  }
  if(operatorMap.find("Diffusion") == operatorMap.end()){
    operatorMap["Diffusion"] = new HDGDiffusion(refEl);
  }
  HDGModel::initializeOperators();
};//initializeOperators


void HDGBurgersModel::computeLocalMatrix(){
  if(!(nodeSet and fieldSet)){
    throw(ErrorHandle("HDGBurgersModel", "computeLocalMatrix", "the nodes and the fields should be set before computing"));
  }
  computeElementJacobians();
  ((HDGBase*)operatorMap["Base"])->setTau(*(fieldMap["Tau"]));
  ((HDGBase*)operatorMap["Base"])->calcNormals(*elementNodes, jacobians);
  operatorMap["Base"]->assemble(dV, invJacobians);
  localMatrix = *(operatorMap["Base"]->getMatrix());
  ((HDGConvection*)operatorMap["Convection"])->setVelocity(parseSolutionVals());
  ((HDGConvection*)operatorMap["Convection"])->setFromBase(((HDGBase*)operatorMap["Base"])->getNormals());
  operatorMap["Convection"]->assemble(dV, invJacobians);
  localMatrix += *(operatorMap["Convection"]->getMatrix());
  localMatrix.block(0, 0, refEl->getNumNodes()*nDOFsPNode, refEl->getNumNodes()*nDOFsPNode) += operatorMap["Convection"]->getMatrix()->block(0, 0, refEl->getNumNodes()*nDOFsPNode, refEl->getNumNodes()*nDOFsPNode).transpose();
  if(fieldMap.find("DiffusionTensor") != fieldMap.end()){
    ((HDGDiffusion*)operatorMap["Diffusion"])->setDiffusionTensor(parseDiffusionVals());
    ((HDGDiffusion*)operatorMap["Diffusion"])->setFromBase(((HDGBase*)operatorMap["Base"])->getNormals());
    operatorMap["Diffusion"]->assemble(dV, invJacobians);
    localMatrix += *(operatorMap["Diffusion"]->getMatrix());
  }
};//computeLocalMatrix

void HDGBurgersModel::computeLocalRHS(){
  if(!(nodeSet and fieldSet)){
    throw("HDGBurgersModel", "computeLocalRHS", "the nodes and the fields should be set before computing");
  }
  localRHS = EVector::Zero(localRHS.size());
  localRHS.segment(0, refEl->getNumNodes()*nDOFsPNode) += operatorMap["Convection"]->getMatrix()->block(0, 0, refEl->getNumNodes()*nDOFsPNode, refEl->getNumNodes()*nDOFsPNode) * EMap<const EVector>(fieldMap["BufferSolution"]->data(), fieldMap["BufferSolution"]->size());
};//computeLocalRHS


std::vector<EMatrix> HDGBurgersModel::parseDiffusionVals() const{
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
    throw(ErrorHandle("HDGBurgersModel", "parseDiffusionVals", "the dimension of the diffusion tensor vales are not correct, they should be either scalar or tensor of the dimension of the reference element"));
  }
  return res;
};//parseDiffusionVals

std::vector<EVector> HDGBurgersModel::parseSolutionVals() const{
  const std::vector<double> * velVals = fieldMap.at("BufferSolution");
  int nNodes = refEl->getNumNodes();
  int dimVelVal = velVals->size()/nNodes;
  int dimSpace = refEl->getDimension();
  if(dimVelVal != dimSpace){
    throw(ErrorHandle("HDGBurgersModel", "parseSolutionVals", "the dimension of the solution buffer vector does not correspond to the dimension of the reference element"));
  }
  std::vector<EVector> res(nNodes, EVector::Zero(dimSpace));
  for(int i = 0; i < nNodes; i++){
    res[i] = EMap<const EVector>(velVals->data() + i*dimSpace, dimSpace);
  }
  return res;
};//parseSolutionVals

};//hfox