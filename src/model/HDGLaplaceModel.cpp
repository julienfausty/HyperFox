#include "HDGLaplaceModel.h"

namespace hfox{

void HDGLaplaceModel::allocate(int nDOFsPerNode){
  assembly.matrix = Add;
  assembly.rhs = None;
  HDGModel::allocate(nDOFsPerNode);
};//allocate

void HDGLaplaceModel::computeLocalMatrix(){
  if(!nodeSet){
    throw(ErrorHandle("HDGLaplaceModel", "computeLocalMatrix", "the nodes should be set before computing matrix."));
  }
  computeElementJacobians();
  ((HDGBase*)operatorMap["Base"])->setTau(*(fieldMap["Tau"]));
  ((HDGBase*)operatorMap["Base"])->calcNormals(*elementNodes, jacobians);
  operatorMap["Base"]->assemble(dV, invJacobians);
  localMatrix = *(operatorMap["Base"]->getMatrix());
  int dim = refEl->getDimension();
  int lenU = refEl->getNumNodes() * nDOFsPNode;
  int startQ = lenU;
  int lenQ = lenU * dim;

  for(int i = 0; i < lenU; i++){
    for(int j = 0; j < lenU; j++){
      for(int k = 0; k < dim; k++){
        localMatrix(i, startQ + j*dim + k) += localMatrix(startQ + i*dim + k, j);
      }
    }
  }

  int nNodesFace = refEl->getFaceElement()->getNumNodes();
  int nFaces = refEl->getNumFaces();
  int startL = lenU + lenQ;
  int lenL = nDOFsPNode * nFaces * nNodesFace;
  const std::vector< std::vector<int> > * nodeMap = refEl->getFaceNodes();
  int startFace = 0;
  for(int iFace = 0; iFace < nFaces; iFace++){
    startFace = iFace * nNodesFace * nDOFsPNode;
    for(int i = 0; i < nNodesFace; i++){
      for(int dof = 0; dof < nDOFsPNode; dof++){
        localMatrix.block((nodeMap->at(iFace))[i]*nDOFsPNode + dof, startQ, 1, lenQ) -= localMatrix.block(startL + startFace + i*nDOFsPNode + dof, startQ, 1, lenQ);
      }
    }
  }
};//computeLocalMatrix

void HDGLaplaceModel::computeLocalRHS(){
  if(!(nodeSet and allocated)){
    throw(ErrorHandle("HDGLaplaceModel", "computeLocalRHS", 
          "the nodes should be set and the object allocated before computing."));
  }
};//computeLocalRHS

}//hfox
