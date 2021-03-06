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
  itfm = fm->find("BufferSolution");
  if(itfm != fm->end()){
    fieldMap["BufferSolution"] = &(itfm->second);
  } else {
    throw(ErrorHandle("HDGnGammaModel", "setFieldMap", "must provide a BufferSolution field for the Newton-Raphson iterations"));
  }
  if(itfm->second.size() != 2*(refEl->getNumNodes())){
    throw(ErrorHandle("HDGnGammaModel", "setFieldMap", "the BufferSolution must have 2 DOFs per node n and Gamma"));
  }
  //itfm = fm->find("Trace");
  //if(itfm != fm->end()){
    //fieldMap["Trace"] = &(itfm->second);
  //} else {
    //throw(ErrorHandle("HDGnGammaModel", "setFieldMap", "must provide a BufferTrace field for the Newton-Raphson iterations"));
  //}
  //if(itfm->second.size() != 2*(refEl->getFaceElement()->getNumNodes())*(refEl->getNumFaces())){
    //throw(ErrorHandle("HDGnGammaModel", "setFieldMap", "the BufferTrace must have 2 DOFs per node"));
  //}
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
  std::vector<EVector> sols = parseSolutionVals();
  std::vector< std::vector<EVector> > vs(4, std::vector<EVector>(refEl->getNumNodes(), EVector::Zero(2)));
  //calculating velocities for terms
  double U1OverU0;
  for(int i = 0; i < refEl->getNumNodes(); i++){
    vs[1][i] = bs[i];
    U1OverU0 = sols[i][1]/sols[i][0];
    vs[2][i] = (std::pow(params.soundSpeed, 2) - std::pow(U1OverU0, 2))*bs[i];
    vs[3][i] = U1OverU0*2.0*bs[i];
  }

  const std::vector<EVector> * ns = ((HDGBase*)operatorMap["Base"])->getNormals();
  Operator * buffOp;
  //convection terms
  for(int i = 0; i < 4; i++){
    buffOp = operatorMap["Convection_"+std::to_string(i)];
    ((HDGConvection*)buffOp)->setFromBase(ns);
    ((HDGConvection*)buffOp)->setVelocity(vs[i]);
    buffOp->assemble(dV, invJacobians);
  }
  //Diffusion terms
  if(fieldMap.find("D") != fieldMap.end()){
    buffOp = operatorMap["Diffusion_0"];
    ((HDGDiffusion*)buffOp)->setFromBase(ns);
    ((HDGDiffusion*)buffOp)->setDiffusionTensor(parseDVals());
    buffOp->assemble(dV, invJacobians);
  }

  if(fieldMap.find("G") != fieldMap.end()){
    buffOp = operatorMap["Diffusion_1"];
    ((HDGDiffusion*)buffOp)->setFromBase(ns);
    ((HDGDiffusion*)buffOp)->setDiffusionTensor(parseGVals());
    buffOp->assemble(dV, invJacobians);
  }

  //combine terms into matrix
  for(int i = 0; i < localMatrix.rows()/2; i++){
    for(int j = 0; j < localMatrix.cols()/2; j++){
      localMatrix(i*2, j*2) += (*(operatorMap["Diffusion_0"]->getMatrix()))(i, j) + (*(operatorMap["Convection_0"]->getMatrix()))(i, j);
      localMatrix(i*2, j*2 + 1) += (*(operatorMap["Convection_1"]->getMatrix()))(i, j);
      localMatrix(i*2 + 1, j*2) += (*(operatorMap["Convection_2"]->getMatrix()))(i, j);
      localMatrix(i*2 + 1, j*2 + 1) += (*(operatorMap["Diffusion_1"]->getMatrix()))(i, j) + (*(operatorMap["Convection_3"]->getMatrix()))(i, j);
    }
  }
};//computeLocalMatrix

void HDGnGammaModel::computeLocalRHS(){
  if(!(nodeSet and fieldSet)){
    throw(ErrorHandle("HDGnGammaModel", "computeLocalRHS", "the nodes and the fields should be set before computing"));
  }
  localRHS = EVector::Zero(localRHS.size());
  if(sourceSet){
    for(int i = 0; i < 2; i++){
      ((Source*)operatorMap["Source_" + std::to_string(i)])->calcSource(*elementNodes);
      operatorMap["Source_" + std::to_string(i)]->assemble(dV, invJacobians);
      for(int iN = 0; iN < refEl->getNumNodes(); iN++){
        localRHS[iN*2 + i] += (*(operatorMap["Source_"+std::to_string(i)]->getMatrix()))(iN, 0);
      }
    }
  }
};//computeLocalRHS

std::vector<EVector> HDGnGammaModel::parseBVals(){
  const std::vector<double> * velVals = fieldMap.at("b");
  int nNodes = refEl->getNumNodes();
  int dimVelVal = velVals->size()/nNodes;
  if(dimVelVal != 2){
    throw(ErrorHandle("HDGnGammaModel", "parseBVals", "the dimension of the b vector is not 2"));
  }
  std::vector<EVector> res(nNodes, EVector::Zero(2));
  for(int i = 0; i < nNodes; i++){
    res[i] = EMap<const EVector>(velVals->data() + i*2, 2);
  }
  return res;
};//parseBVals

std::vector<EMatrix> HDGnGammaModel::parseDVals(){
  const std::vector<double> * diffVals = fieldMap.at("D");
  int nNodes = refEl->getNumNodes();
  int dimDiffVal = diffVals->size()/nNodes;
  int dimMat = std::sqrt(dimDiffVal);
  std::vector<EMatrix> res(nNodes, EMatrix::Identity(2, 2));
  if(dimMat == 1){
    for(int i = 0; i < nNodes; i++){
      res[i] *= (diffVals->at(i));
    }
  } else if(dimMat == 2){
    for(int i = 0; i < nNodes; i++){
      res[i] = EMap<const EMatrix>(diffVals->data() + i*dimDiffVal, dimMat, dimMat);
    }
  } else{
    throw(ErrorHandle("HDGnGammaModel", "parseDVals", "the dimension of the diffusion tensor vales are not correct, they should be either scalar or tensor of the dimension of 2"));
  }
  return res;
};//parseDVals

std::vector<EMatrix> HDGnGammaModel::parseGVals(){
  const std::vector<double> * diffVals = fieldMap.at("G");
  int nNodes = refEl->getNumNodes();
  int dimDiffVal = diffVals->size()/nNodes;
  int dimMat = std::sqrt(dimDiffVal);
  std::vector<EMatrix> res(nNodes, EMatrix::Identity(2, 2));
  if(dimMat == 1){
    for(int i = 0; i < nNodes; i++){
      res[i] *= (diffVals->at(i));
    }
  } else if(dimMat == 2){
    for(int i = 0; i < nNodes; i++){
      res[i] = EMap<const EMatrix>(diffVals->data() + i*dimDiffVal, dimMat, dimMat);
    }
  } else{
    throw(ErrorHandle("HDGnGammaModel", "parseGVals", "the dimension of the diffusion tensor vales are not correct, they should be either scalar or tensor of the dimension of 2"));
  }
  return res;
};//parseGVals

std::vector<EVector> HDGnGammaModel::parseSolutionVals(){
  const std::vector<double> * velVals = fieldMap.at("BufferSolution");
  int nNodes = refEl->getNumNodes();
  int dimVelVal = velVals->size()/nNodes;
  if(dimVelVal != 2){
    throw(ErrorHandle("HDGnGammaModel", "parseSolutionVals", "the dimension of the BufferSolution vector is not 2"));
  }
  std::vector<EVector> res(nNodes, EVector::Zero(2));
  for(int i = 0; i < nNodes; i++){
    res[i] = EMap<const EVector>(velVals->data() + i*2, 2);
  }
  return res;
};//parseSolutionVals

};//nGamma

};//hfox
