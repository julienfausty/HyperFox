#include "HDGTransport.h"

namespace hfox{

void HDGTransport::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator itfm;
  itfm = fm->find("Velocity");
  if(itfm == fm->end()){
    throw(ErrorHandle("HDGTransport", "setFieldMap", "one must provide a Velocity field to use the Transport model."));
  }
  fieldMap["Velocity"] = &(itfm->second);
  if(timeScheme != NULL){
    timeScheme->setFieldMap(fm);
  }
  HDGModel::setFieldMap(fm);
};//setFieldMap

void HDGTransport::allocate(int nDOFsPerNodeUser){
  assembly.matrix = Add;
  assembly.rhs = Add;
  if(timeScheme != NULL){
    timeScheme->allocate(nDOFsPerNodeUser);
  }
  HDGModel::allocate(nDOFsPerNodeUser);
};//allocate

void HDGTransport::initializeOperators(){
  if(operatorMap.find("Convection") == operatorMap.end()){
    operatorMap["Convection"] = new HDGConvection(refEl);
  }
  HDGModel::initializeOperators();
};//initializeOperators

std::vector<EVector> HDGTransport::parseVelocityVals() const{
  const std::vector<double> * velVals = fieldMap.at("Velocity");
  int nNodes = refEl->getNumNodes();
  int dimVelVal = velVals->size()/nNodes;
  int dimSpace = refEl->getDimension();
  if(dimVelVal != dimSpace){
    throw(ErrorHandle("HDGTransport", "parseVelocityVals", "the dimension of the velocity vector does not correspond to the dimension of the reference element"));
  }
  std::vector<EVector> res(nNodes, EVector::Zero(dimSpace));
  for(int i = 0; i < nNodes; i++){
    res[i] = EMap<const EVector>(velVals->data() + i*dimSpace, dimSpace);
  }
  return res;
};//parseVelocityVals

void HDGTransport::computeLocalMatrix(){
  if(!(nodeSet and fieldSet)){
    throw(ErrorHandle("HDGTransport", "computeLocalMatrix", "the nodes and the fields should be set before computing"));
  }
  computeElementJacobians();
  ((HDGBase*)operatorMap["Base"])->setTau(*(fieldMap["Tau"]));
  ((HDGBase*)operatorMap["Base"])->calcNormals(*elementNodes, jacobians);
  ((HDGConvection*)operatorMap["Convection"])->setVelocity(parseVelocityVals());
  operatorMap["Base"]->assemble(dV, invJacobians);
  operatorMap["Convection"]->assemble(dV, invJacobians);
  localMatrix = *(operatorMap["Base"]->getMatrix());
  localMatrix += *(operatorMap["Convection"]->getMatrix());
};//computeLocalMatrix

void HDGTransport::computeLocalRHS(){
  if(!(nodeSet and fieldSet)){
    throw("HDGTransport", "computeLocalRHS", "the nodes and the fields should be set before computing");
  }
  localRHS = EVector::Zero(localRHS.size());
};//computeLocalRHS

}//hfox
