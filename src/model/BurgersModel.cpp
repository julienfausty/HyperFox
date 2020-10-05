#include "BurgersModel.h"

namespace hfox{

void BurgersModel::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator itfm;
  itfm = fm->find("DiffusionTensor");
  if(itfm != fm->end()){
    fieldMap["DiffusionTensor"] = &(itfm->second);
  }
  itfm = fm->find("BufferSolution");
  if(itfm != fm->end()){
    fieldMap["BufferSolution"] = &(itfm->second);
  } else {
    throw(ErrorHandle("BurgersModel", "setFieldMap", "must provide a BufferSolution field for the Newton-Raphson iterations"));
  }
  if(itfm->second.size() != refEl->getDimension()*(refEl->getNumNodes())){
    throw(ErrorHandle("BurgersModel", "setFieldMap", "the BufferSolution must have the same number of DOFs per node as spatial dimensions for the Burgers equation"));
  }
  if(timeScheme != NULL){
    timeScheme->setFieldMap(fm);
  }
  fieldSet = 1;
};//setFieldMap


void BurgersModel::allocate(int nDOFsPerNode){
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
    if(it->first.find("Source") == std::string::npos){
      (it->second)->allocate(nDOFsPerNode);
    } else {
      (it->second)->allocate(1);
    }
  }
  allocated = 1;
};//allocate

void BurgersModel::initializeOperators(){
  if(operatorMap.find("Convection") == operatorMap.end()){
    operatorMap["Convection"] = new UNabU(refEl);
  }
  if(operatorMap.find("Diffusion") == operatorMap.end()){
    operatorMap["Diffusion"] = new Diffusion(refEl);
  }
  for(int i = 0; i < refEl->getDimension(); i++){
    if(operatorMap.find("Source_" + std::to_string(i)) == operatorMap.end()){
      operatorMap["Source_" + std::to_string(i)] = new Source(refEl);
    }
  }
};//initializeOperators

void BurgersModel::setSourceFunction(std::function<double(const std::vector<double>&, int)> s){
  if(!allocated){
    throw(ErrorHandle("BurgersModel", "setSourceFunction", "the model must be allocated before setting the source function"));
  }
  for(int i = 0; i < refEl->getDimension(); i++){
    ((Source*)operatorMap["Source_" + std::to_string(i)])->setSourceFunction([i,s](const std::vector<double> & x){return s(x, i);});
  }
  sourceSet = 1;
};//setSourceFunction


void BurgersModel::computeLocalMatrix(){
  if(!(nodeSet and fieldSet)){
    throw(ErrorHandle("BurgersModel", "computeLocalMatrix", "the nodes and the fields should be set before computing"));
  }
  jacobians = Operator::calcJacobians(*elementNodes, refEl);
  invJacobians = Operator::calcInvJacobians(jacobians);
  dV = Operator::calcMeasure(Operator::calcDetJacobians(jacobians), refEl);
  ((UNabU*)operatorMap["Convection"])->setSolution(parseSolutionVals());
  operatorMap["Convection"]->assemble(dV, invJacobians);
  localMatrix = *(operatorMap["Convection"]->getMatrix());
  if(fieldMap.find("DiffusionTensor") != fieldMap.end()){
    ((Diffusion*)operatorMap["Diffusion"])->setDiffTensor(parseDiffusionVals());
    operatorMap["Diffusion"]->assemble(dV, invJacobians);
    localMatrix += *(operatorMap["Diffusion"]->getMatrix());
  }
};//computeLocalMatrix

void BurgersModel::computeLocalRHS(){
  if(!(nodeSet and fieldSet)){
    throw("HDGBurgersModel", "computeLocalRHS", "the nodes and the fields should be set before computing");
  }
  localRHS = *(((UNabU*)operatorMap["Convection"])->getRHS());
  if(sourceSet){
    int dim = refEl->getDimension();
    for(int i = 0; i < dim; i++){
      ((Source*)operatorMap["Source_" + std::to_string(i)])->calcSource(*elementNodes);
      operatorMap["Source_" + std::to_string(i)]->assemble(dV, invJacobians);
      for(int iN = 0; iN < refEl->getNumNodes(); iN++){
        localRHS[iN*dim + i] += (*(operatorMap["Source_"+std::to_string(i)]->getMatrix()))(iN, 0);
      }
    }
  }
};//computeLocalRHS

std::vector<EMatrix> BurgersModel::parseDiffusionVals() const{
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
    throw(ErrorHandle("BurgersModel", "parseDiffusionVals", "the dimension of the diffusion tensor vales are not correct, they should be either scalar or tensor of the dimension of the reference element"));
  }
  return res;
};//parseDiffusionVals

std::vector<EVector> BurgersModel::parseSolutionVals() const{
  const std::vector<double> * velVals = fieldMap.at("BufferSolution");
  int nNodes = refEl->getNumNodes();
  int dimVelVal = velVals->size()/nNodes;
  int dimSpace = refEl->getDimension();
  if(dimVelVal != dimSpace){
    throw(ErrorHandle("BurgersModel", "parseSolutionVals", "the dimension of the solution buffer vector does not correspond to the dimension of the reference element"));
  }
  std::vector<EVector> res(nNodes, EVector::Zero(dimSpace));
  for(int i = 0; i < nNodes; i++){
    res[i] = EMap<const EVector>(velVals->data() + i*dimSpace, dimSpace);
  }
  return res;
};//parseSolutionVals

};//hfox
