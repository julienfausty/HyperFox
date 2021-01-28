#include "HDGLaplaceModel.h"

namespace hfox{

void HDGLaplaceModel::allocate(int nDOFsPerNode){
  assembly.matrix = Add;
  assembly.rhs = Add;
  HDGModel::allocate(nDOFsPerNode);
};//allocate

void HDGLaplaceModel::initializeOperators(){
  if(operatorMap.find("Diffusion") == operatorMap.end()){
    operatorMap["Diffusion"] = new HDGDiffusion(refEl);
  }
  HDGModel::initializeOperators();
};//initializeOperators

void HDGLaplaceModel::computeLocalMatrix(){
  if(!nodeSet){
    throw(ErrorHandle("HDGLaplaceModel", "computeLocalMatrix", "the nodes should be set before computing matrix."));
  }
  computeElementJacobians();
  ((HDGBase*)operatorMap["Base"])->setTau(*(fieldMap["Tau"]));
  ((HDGBase*)operatorMap["Base"])->calcNormals(*elementNodes, jacobians);
  operatorMap["Base"]->assemble(dV, invJacobians);
  localMatrix = *(operatorMap["Base"]->getMatrix());
  ((HDGDiffusion*)operatorMap["Diffusion"])->setFromBase(((HDGBase*)operatorMap["Base"])->getNormals());
  operatorMap["Diffusion"]->assemble(dV, invJacobians);
  localMatrix += *(operatorMap["Diffusion"]->getMatrix());
};//computeLocalMatrix

void HDGLaplaceModel::computeLocalRHS(){
  if(!(nodeSet and allocated)){
    throw(ErrorHandle("HDGLaplaceModel", "computeLocalRHS", 
          "the nodes should be set and the object allocated before computing."));
  }
};//computeLocalRHS

}//hfox
