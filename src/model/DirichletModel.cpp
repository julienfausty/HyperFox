#include "DirichletModel.h"

namespace hfox{

void DirichletModel::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
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

void DirichletModel::allocate(int nDOFsPerNode){
  assembly.matrix = Set;
  assembly.rhs = Set;
  localMatrix = EMatrix::Identity(nDOFsPerNode*refEl->getNumNodes(), nDOFsPerNode*refEl->getNumNodes());
  localRHS = EVector::Zero(nDOFsPerNode*refEl->getNumNodes());
  allocated = 1;
};//allocate

void DirichletModel::initializeOperators(){
};//initializeOperators

void DirichletModel::computeLocalMatrix(){
  if(!(fieldSet and allocated)){
    throw(ErrorHandle("DirichletModel", "computeLocalMatrix", 
          "the fields should be set and the object allocated before computing."));
  }
};//computeLocalMatrix


void DirichletModel::computeLocalRHS(){
  if(!(allocated and fieldSet)){
    throw(ErrorHandle("DirichletModel", "computeLocalRHS", 
          "the fields should be set and the object allocated before computing."));
  }
  localRHS = Eigen::Map<const EVector>(fieldMap["Dirichlet"]->data(), fieldMap["Dirichlet"]->size());
};//computeLocalRHS

}
