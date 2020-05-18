#include "HDGModel.h"

namespace hfox{

void HDGModel::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator itfm;
  itfm = fm->find("Tau");
  if(itfm != fm->end()){
    fieldMap["Tau"] = &(itfm->second);
  } else {
    throw(ErrorHandle("HDGModel", "setFieldMap", "the Tau field must be present in the field map."));
  }
  fieldSet = 1;
};//setFieldMap

void HDGModel::allocate(int nDOFsPerNode){
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

void HDGModel::initializeOperators(){
  if(operatorMap.find("Base") == operatorMap.end()){
    operatorMap["Base"] = new HDGBase(refEl);
  }
};//initializeOperators

void HDGModel::compute(){
  if(allocated){
    computeLocalMatrix();
    computeLocalRHS();
    if(timeScheme != NULL){
      int lenU = refEl->getNumNodes()*nDOFsPNode;
      EMatrix Su = localMatrix.block(0, 0, lenU, localMatrix.cols());
      EVector Fu = localRHS.segment(0, lenU);
      timeScheme->assemble(dV, invJacobians);
      timeScheme->apply(&Su, &Fu);
      localMatrix.block(0, 0, lenU, localMatrix.cols()) = Su;
      localRHS.segment(0, lenU) = Fu;
    }
  } else {
    throw(ErrorHandle("HDGModel", "compute", "the model must be allocated before computing."));
  }
};//compute

void HDGModel::computeElementJacobians(){
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

}//hfox
