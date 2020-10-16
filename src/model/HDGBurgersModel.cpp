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
  itfm = fm->find("Trace");
  if(itfm != fm->end()){
    fieldMap["Trace"] = &(itfm->second);
  } else {
    throw(ErrorHandle("HDGBurgersModel", "setFieldMap", "must provide a Trace field for the Newton-Raphson iterations"));
  }
  if(itfm->second.size() != refEl->getDimension()*(refEl->getFaceElement()->getNumNodes())*(refEl->getNumFaces())){
    throw(ErrorHandle("HDGBurgersModel", "setFieldMap", "the Trace must have the same number of DOFs per node as spatial dimensions for the Burgers equation"));
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
  initializeOperators();
  nDOFsPNode = nDOFsPerNodeUser;
  int n = nDOFsPNode * (refEl->getNumNodes()*(refEl->getDimension() + 1) 
      + (refEl->getFaceElement()->getNumNodes())*(refEl->getNumFaces()));
  localMatrix = EMatrix::Zero(n, n);
  localRHS = EVector::Zero(n);
  for(std::map<std::string, Operator*>::iterator it = operatorMap.begin(); it != operatorMap.end(); it++){
    if(it->first.find("Source") == std::string::npos){
      (it->second)->allocate(nDOFsPNode);
    } else {
      (it->second)->allocate(1);
    }
  }
  allocated = 1;
};//allocate


void HDGBurgersModel::initializeOperators(){
  if(operatorMap.find("Convection") == operatorMap.end()){
    operatorMap["Convection"] = new HDGUNabU(refEl);
  }
  if(operatorMap.find("Diffusion") == operatorMap.end()){
    operatorMap["Diffusion"] = new HDGDiffusion(refEl);
  }
  for(int i = 0; i < refEl->getDimension(); i++){
    if(operatorMap.find("Source_" + std::to_string(i)) == operatorMap.end()){
      operatorMap["Source_" + std::to_string(i)] = new Source(refEl);
    }
  }
  HDGModel::initializeOperators();
};//initializeOperators

void HDGBurgersModel::setSourceFunction(std::function<double(const std::vector<double>&, int)> s){
  if(!allocated){
    throw(ErrorHandle("BurgersModel", "setSourceFunction", "the model must be allocated before setting the source function"));
  }
  for(int i = 0; i < refEl->getDimension(); i++){
    ((Source*)operatorMap["Source_" + std::to_string(i)])->setSourceFunction([i,s](const std::vector<double> & x){return s(x, i);});
  }
  sourceSet = 1;
};//setSourceFunction


void HDGBurgersModel::computeLocalMatrix(){
  if(!(nodeSet and fieldSet)){
    throw(ErrorHandle("HDGBurgersModel", "computeLocalMatrix", "the nodes and the fields should be set before computing"));
  }
  computeElementJacobians();
  ((HDGBase*)operatorMap["Base"])->setTau(*(fieldMap["Tau"]));
  ((HDGBase*)operatorMap["Base"])->calcNormals(*elementNodes, jacobians);
  operatorMap["Base"]->assemble(dV, invJacobians);
  localMatrix = *(operatorMap["Base"]->getMatrix());
  //std::cout << "BaseMat:\n" << *(operatorMap["Base"]->getMatrix()) << std::endl;
  ((HDGUNabU*)operatorMap["Convection"])->setTrace(parseTraceVals());
  ((HDGUNabU*)operatorMap["Convection"])->setSolution(parseSolutionVals());
  ((HDGUNabU*)operatorMap["Convection"])->setFromBase(((HDGBase*)operatorMap["Base"])->getNormals());
  operatorMap["Convection"]->assemble(dV, invJacobians);
  localMatrix += *(operatorMap["Convection"]->getMatrix());
  //std::cout << "ConvectionMat:\n" << *(operatorMap["Convection"]->getMatrix()) << std::endl;
  if(fieldMap.find("DiffusionTensor") != fieldMap.end()){
    ((HDGDiffusion*)operatorMap["Diffusion"])->setDiffusionTensor(parseDiffusionVals());
    ((HDGDiffusion*)operatorMap["Diffusion"])->setFromBase(((HDGBase*)operatorMap["Base"])->getNormals());
    operatorMap["Diffusion"]->assemble(dV, invJacobians);
    //std::cout << "DiffusionMat:\n" << *(operatorMap["Diffusion"]->getMatrix()) << std::endl;
    localMatrix += *(operatorMap["Diffusion"]->getMatrix());
  }
};//computeLocalMatrix

void HDGBurgersModel::computeLocalRHS(){
  if(!(nodeSet and fieldSet)){
    throw("HDGBurgersModel", "computeLocalRHS", "the nodes and the fields should be set before computing");
  }
  //localRHS = EVector::Zero(localRHS.size());
  localRHS = *(((HDGUNabU*)operatorMap["Convection"])->getRHS()); 
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


std::vector<EVector> HDGBurgersModel::parseTraceVals() const{
  const std::vector<double> * traceVals = fieldMap.at("Trace");
  int nNodes = refEl->getFaceElement()->getNumNodes() * (refEl->getNumFaces());
  int dimTraceVal = traceVals->size()/nNodes;
  int dimSpace = refEl->getDimension();
  if(dimTraceVal != dimSpace){
    throw(ErrorHandle("HDGBurgersModel", "parseTraceVals", "the dimension of the Trace vector does not correspond to the dimension of the reference element"));
  }
  std::vector<EVector> res(nNodes, EVector::Zero(dimSpace));
  for(int i = 0; i < nNodes; i++){
    res[i] = EMap<const EVector>(traceVals->data() + i*dimSpace, dimSpace);
  }
  return res;
};//parseTraceVals

};//hfox
