#include "BohmModel.h"

namespace hfox{

namespace nGamma{

void BohmModel::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator it;
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "b"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("b") == fieldMap.end()){
    throw(ErrorHandle("BohmModel", "setFieldMap", 
          "need to give a field named b to the BohmModel"));
  }
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "SoundVelocity"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("SoundVelocity") == fieldMap.end()){
    throw(ErrorHandle("BohmModel", "setFieldMap", 
          "need to give a field named SoundVelocity to the BohmModel"));
  }
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "ExteriorNormals"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("ExteriorNormals") == fieldMap.end()){
    throw(ErrorHandle("BohmModel", "setFieldMap", 
          "need to give a field named ExteriorNormals to the BohmModel"));
  }
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "Trace"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("Trace") == fieldMap.end()){
    throw(ErrorHandle("BohmModel", "setFieldMap", 
          "need to give a field named Trace to the BohmModel"));
  }
  fieldSet = 1;
};//setFieldMap

void BohmModel::allocate(int nDOFsPerNode){
  assembly.matrix = Set;
  assembly.rhs = Set;
  initializeOperators();
  localMatrix = EMatrix::Zero(nDOFsPerNode*refEl->getNumNodes(), nDOFsPerNode*refEl->getNumNodes()*4);
  localRHS = EVector::Zero(nDOFsPerNode*refEl->getNumNodes());
  for(std::map<std::string, Operator*>::iterator it = operatorMap.begin(); it != operatorMap.end(); it++){
    (it->second)->allocate(1);
  }
  allocated = 1;
};//allocate

void BohmModel::initializeOperators(){
  if(operatorMap.find("Mass") == operatorMap.end()){
    operatorMap["Mass"] = new Mass(refEl);
  }
};//initializeOperators

void BohmModel::computeLocalMatrix(){
  if(!(fieldSet and allocated)){
    throw(ErrorHandle("BohmModel", "computeLocalMatrix", 
          "the fields should be set and the object allocated before computing."));
  }
  jacobians = Operator::calcJacobians(*elementNodes, refEl);
  invJacobians = Operator::calcInvJacobians(jacobians);
  dV = Operator::calcMeasure(Operator::calcDetJacobians(jacobians), refEl);
  operatorMap["Mass"]->assemble(dV, invJacobians);
  localMatrix = EMatrix::Zero(localMatrix.rows(), localMatrix.cols());
  int nNodes = refEl->getNumNodes();
  const std::vector< std::vector<double> > * shapes = refEl->getIPShapeFunctions();
  double dbuff;
  std::vector<EVector> normals(refEl->getNumIPs(), EVector::Zero(2));
  EMatrix normalMat = EMatrix::Zero(2, nNodes);
  for(int ndof = 0; ndof < nNodes; ndof++){
    normalMat.col(ndof) = EMap<const EVector>(fieldMap["ExteriorNormals"]->data() + ndof*2, 2);
  }
  for(int ip = 0; ip < refEl->getNumIPs(); ip++){
    normals[ip] = normalMat * EMap<const EVector>(shapes->at(ip).data(), shapes->at(ip).size());
  }
  //n condition
  for(int ip = 0; ip < refEl->getNumIPs(); ip++){
    for(int ndof = 0; ndof < nNodes; ndof++){
      for(int mdof = 0; mdof < nNodes; mdof++){
        dbuff = dV[ip]*(shapes->at(ip)[ndof])*(shapes->at(ip)[mdof]);
        for(int d = 0; d < 2; d++){
          localMatrix(mdof*2, nNodes + (ndof*2 + d)*2) += dbuff*normals[ip][d];
        }
      }
    }
  }
  //Gamma condition

};//computeLocalMatrix


void BohmModel::computeLocalRHS(){
  if(!(allocated and fieldSet)){
    throw(ErrorHandle("BohmModel", "computeLocalRHS", 
          "the fields should be set and the object allocated before computing."));
  }
};//computeLocalRHS

};//nGamma

};//hfox
