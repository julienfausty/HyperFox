#include "Transport.h"

namespace hfox{

void Transport::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator itfm;
  itfm = fm->find("Velocity");
  if(itfm == fm->end()){
    throw(ErrorHandle("Transport", "setFieldMap", "one must provide a Velocity field to use the Transport model."));
  }
  fieldMap["Velocity"] = &(itfm->second);
  if(timeScheme != NULL){
    timeScheme->setFieldMap(fm);
  }
  fieldSet = 1;
};//setFieldMap

void Transport::allocate(int nDOFsPerNode){
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

void Transport::initializeOperators(){
  if(operatorMap.find("Convection") == operatorMap.end()){
    operatorMap["Convection"] = new Convection(refEl);
  }
};//initializeOperators

std::vector<EVector> Transport::parseVelocityVals() const{
  const std::vector<double> * velVals = fieldMap.at("Velocity");
  int nNodes = refEl->getNumNodes();
  int dimVelVal = velVals->size()/nNodes;
  int dimSpace = refEl->getDimension();
  if(dimVelVal != dimSpace){
    throw(ErrorHandle("Transport", "parseVelocityVals", "the dimension of the velocity vector does not correspond to the dimension of the reference element"));
  }
  std::vector<EVector> res(nNodes, EVector::Zero(dimSpace));
  for(int i = 0; i < nNodes; i++){
    res[i] = EMap<const EVector>(velVals->data() + i*dimSpace, dimSpace);
  }
  return res;
};//parseVelocityVals

void Transport::computeLocalMatrix(){
  if(!(nodeSet and fieldSet)){
    throw(ErrorHandle("Transport", "computeLocalMatrix", "the nodes and the fields should be set before computing"));
  }
  jacobians = Operator::calcJacobians(*elementNodes, refEl);
  invJacobians = Operator::calcInvJacobians(jacobians);
  dV = Operator::calcMeasure(Operator::calcDetJacobians(jacobians), refEl);
  ((Convection*)operatorMap["Convection"])->setVelocity(parseVelocityVals());
  operatorMap["Convection"]->assemble(dV, invJacobians);
  localMatrix = *(operatorMap["Convection"]->getMatrix());
};//computeLocalMatrix

void Transport::computeLocalRHS(){
  if(!(nodeSet and fieldSet)){
    throw("Transport", "computeLocalRHS", "the nodes and the fields should be set before computing");
  }
  localRHS = EVector::Zero(localRHS.size());
};//computeLocalRHS

}//hfox
