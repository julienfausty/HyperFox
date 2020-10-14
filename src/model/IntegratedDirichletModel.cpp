#include "IntegratedDirichletModel.h"

namespace hfox{

void IntegratedDirichletModel::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator it;
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "Dirichlet"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("Dirichlet") == fieldMap.end()){
    throw(ErrorHandle("DirichletModel", "setFieldMap", 
          "need to give a field named Dirichlet to the DirichletModel"));
  }
  fieldSet = 1;
};//setFieldMap

void IntegratedDirichletModel::allocate(int nDOFsPerNode){
  assembly.matrix = Set;
  assembly.rhs = Set;
  initializeOperators();
  localMatrix = EMatrix::Zero(nDOFsPerNode*refEl->getNumNodes(), nDOFsPerNode*refEl->getNumNodes());
  localRHS = EVector::Zero(nDOFsPerNode*refEl->getNumNodes());
  for(std::map<std::string, Operator*>::iterator it = operatorMap.begin(); it != operatorMap.end(); it++){
    (it->second)->allocate(nDOFsPerNode);
  }
  allocated = 1;
};//allocate

void IntegratedDirichletModel::initializeOperators(){
  if(operatorMap.find("Mass") == operatorMap.end()){
    operatorMap["Mass"] = new Mass(refEl);
  }
};//initializeOperators

void IntegratedDirichletModel::computeLocalMatrix(){
  if(!(fieldSet and allocated)){
    throw(ErrorHandle("IntegratedDirichletModel", "computeLocalMatrix", 
          "the fields should be set and the object allocated before computing."));
  }
  jacobians = Operator::calcJacobians(*elementNodes, refEl);
  invJacobians = Operator::calcInvJacobians(jacobians);
  dV = Operator::calcMeasure(Operator::calcDetJacobians(jacobians), refEl);
  operatorMap["Mass"]->assemble(dV, invJacobians);
  localMatrix = *(operatorMap["Mass"]->getMatrix());
};//computeLocalMatrix


void IntegratedDirichletModel::computeLocalRHS(){
  if(!(allocated and fieldSet)){
    throw(ErrorHandle("IntegratedDirichletModel", "computeLocalRHS", 
          "the fields should be set and the object allocated before computing."));
  }
  localRHS = localMatrix*Eigen::Map<const EVector>(fieldMap["Dirichlet"]->data(), fieldMap["Dirichlet"]->size());
};//computeLocalRHS

};//hfox
