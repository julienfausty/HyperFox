#include "HDGLaplaceModel.h"

namespace hfox{

void HDGLaplaceModel::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator itfm;
  itfm = fm->find("Tau");
  if(itfm != fm->end()){
    fieldMap["Tau"] = &(itfm->second);
  } else {
    throw(ErrorHandle("HDGLaplaceModel", "setFieldMap", "the Tau field must be present in the field map."));
  }
  fieldSet = 1;
};//setFieldMap

void HDGLaplaceModel::allocate(int nDOFsPerNode){
  assembly.matrix = Add;
  assembly.rhs = None;
  initializeOperators();
  nDOFsPNode = nDOFsPerNode;
  int n = nDOFsPerNode * (refEl->getNumNodes()*(refEl->getDimension() + 1) 
      + (refEl->getFaceElement()->getNumNodes())*(refEl->getNumFaces()));
  localMatrix = EMatrix::Zero(n, n);
  localRHS = EVector::Zero(n);
  for(std::map<std::string, Operator*>::iterator it = operatorMap.begin(); it != operatorMap.end(); it++){
    (it->second)->allocate(nDOFsPerNode);
  }
  allocated = 1;
};//allocate

void HDGLaplaceModel::initializeOperators(){
  operatorMap["Base"] = new HDGBase(refEl);
};//initializeOperators

void HDGLaplaceModel::computeElementJacobians(){
  const ReferenceElement * fEl = refEl->getFaceElement();
  int nFaces = refEl->getNumFaces();
  int nNodesFace = fEl->getNumIPs();
  int nNodesEl = refEl->getNumIPs();
  int numJacs = nNodesEl + nNodesFace*nFaces;
  jacobians.resize(numJacs);
  invJacobians.resize(numJacs);
  dV.resize(numJacs, 0.0);
  std::vector<EMatrix> elMats = Operator::calcJacobians(*elementNodes, refEl);
  std::copy(elMats.begin(), elMats.end(), jacobians.begin());
  std::vector<double> eldV = Operator::calcMeasure(Operator::calcDetJacobians(elMats), refEl);
  std::copy(eldV.begin(), eldV.end(), dV.begin());
  elMats = Operator::calcInvJacobians(elMats);
  std::copy(elMats.begin(), elMats.end(), invJacobians.begin());
  std::vector<EMatrix> faceMats(nNodesFace);
  std::vector<double> facedV(nNodesFace, 0.0);
  std::vector< std::vector<double> > faceNodes(fEl->getNumNodes(), std::vector<double>(elementNodes->at(0).size(), 0.0));
  const std::vector< std::vector<int> > * faceNodeMap = refEl->getFaceNodes();
  int offset = 0;
  for(int iFace = 0; iFace < nFaces; iFace++){
    offset = nNodesEl + iFace*nNodesFace;
    for(int i = 0; i < fEl->getNumNodes(); i++){
      faceNodes[i] = elementNodes->at(faceNodeMap->at(iFace)[i]);
    }
    faceMats = Operator::calcJacobians(faceNodes, fEl);
    std::copy(faceMats.begin(), faceMats.end(), jacobians.begin() + offset);
    facedV = Operator::calcMeasure(Operator::calcDetJacobians(faceMats), fEl);
    std::copy(facedV.begin(), facedV.end(), dV.begin() + offset);
    faceMats = Operator::calcInvJacobians(faceMats);
    std::copy(faceMats.begin(), faceMats.end(), invJacobians.begin() + offset);
  }
};//computElementJacobians

void HDGLaplaceModel::computeLocalMatrix(){
  if(!nodeSet){
    throw(ErrorHandle("HDGLaplaceModel", "computeLocalMatrix", "the nodes should be set before computing matrix."));
  }
  computeElementJacobians();
  ((HDGBase*)operatorMap["Base"])->setTau(*(fieldMap["Tau"]));
  ((HDGBase*)operatorMap["Base"])->calcNormals(*elementNodes, jacobians);
  operatorMap["Base"]->assemble(dV, invJacobians);
  localMatrix = *(operatorMap["Base"]->getMatrix());
  int lenU = refEl->getNumNodes() * nDOFsPNode;
  int startQ = lenU;
  int lenQ = lenU * (refEl->getDimension());

  localMatrix.block(0, startQ, lenU, lenQ) += localMatrix.block(startQ, 0, lenQ, lenU).transpose();

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
