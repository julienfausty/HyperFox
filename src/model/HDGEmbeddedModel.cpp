#include "HDGEmbeddedModel.h"

namespace hfox{

void HDGEmbeddedModel::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator itfm;
  itfm = fm->find("Tau");
  if(itfm != fm->end()){
    fieldMap["Tau"] = &(itfm->second);
  } else {
    throw(ErrorHandle("HDGEmbeddedModel", "setFieldMap", "the Tau field must be present in the field map."));
  }
  itfm = fm->find("Metric");
  if(itfm != fm->end()){
    fieldMap["Metric"] = &(itfm->second);
  } else {
    throw(ErrorHandle("HDGEmbeddedModel", "setFieldMap", "the Metric field must be present in the field map (computed using IP2NodesSolver)."));
  }
  fieldSet = 1;
  itfm = fm->find("Jacobian");
  if(itfm != fm->end()){
    fieldMap["Jacobian"] = &(itfm->second);
  } else {
    throw(ErrorHandle("HDGEmbeddedModel", "setFieldMap", "the Jacobian field must be present in the field map (computed using IP2NodesSolver)."));
  }
  fieldSet = 1;
};//setFieldMap

void HDGEmbeddedModel::allocate(int nDOFsPerNode){
  if(embeddingDim == 0){
    throw(ErrorHandle("HDGEmbeddedModel", "allocate", "the embedding dimension of the model should be set before allocating."));
  }
  initializeOperators();
  nDOFsPNode = nDOFsPerNode;
  int n = nDOFsPerNode * (refEl->getNumNodes()*(embeddingDim + 1) 
      + (refEl->getFaceElement()->getNumNodes())*(refEl->getNumFaces()));
  localMatrix = EMatrix::Zero(n, n);
  localRHS = EVector::Zero(n);
  for(std::map<std::string, Operator*>::iterator it = operatorMap.begin(); it != operatorMap.end(); it++){
    (it->second)->allocate(nDOFsPerNode);
  }
  allocated = 1;
};//allocate

void HDGEmbeddedModel::setEmbeddingDimension(int dim){
  int refDim = refEl->getDimension();
  if(dim < refDim){
    throw(ErrorHandle("HDGEmbeddedModel", "setEmbeddingDimension", "embedding dimension must be higher or equal to the reference space dimension."));
  }
  embeddingDim = dim;
};//setEmbeddingDimension

void HDGEmbeddedModel::initializeOperators(){
  if(operatorMap.find("Base") == operatorMap.end()){
    //operatorMap["Base"] = new HDGEmbeddedBase(refEl);
  }
};//initializeOperators

void HDGEmbeddedModel::compute(){
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
    throw(ErrorHandle("HDGEmbeddedModel", "compute", "the model must be allocated before computing."));
  }
};//compute

void HDGEmbeddedModel::computeElementGeometry(){
  int refDim = refEl->getDimension();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int nNodesEl = refEl->getNumNodes();
  int nIPsEl = refEl->getNumIPs();
  int nFaces = refEl->getNumFaces();
  int nNodesFc = fEl->getNumNodes();
  int nIPsFc = fEl->getNumIPs();
  //jacobians
  int nJacs = nIPsEl + nIPsFc*nFaces;
  jacobians.resize(nJacs);
  dV.resize(nJacs, 0.0);
  std::vector<EMatrix> elMats = Operator::calcJacobians(*elementNodes, refEl);
  std::copy(elMats.begin(), elMats.end(), jacobians.begin());
  std::vector<double> eldV = Operator::calcMeasure(Operator::calcDetJacobians(elMats), refEl);
  std::copy(eldV.begin(), eldV.end(), dV.begin());
  std::vector<EMatrix> faceMats(nIPsFc);
  std::vector<double> facedV(nIPsFc, 0.0);
  std::vector< std::vector<double> > faceNodes(nNodesFc, std::vector<double>(embeddingDim, 0.0));
  const std::vector< std::vector<int> > * faceNodeMap = refEl->getFaceNodes();
  int offset = 0;
  for(int iFace = 0; iFace < nFaces; iFace++){
    offset = nNodesEl + iFace*nIPsFc;
    for(int i = 0; i < nNodesFc; i++){
      faceNodes[i] = elementNodes->at(faceNodeMap->at(iFace)[i]);
    }
    faceMats = Operator::calcJacobians(faceNodes, fEl);
    std::copy(faceMats.begin(), faceMats.end(), jacobians.begin() + offset);
    facedV = Operator::calcMeasure(Operator::calcDetJacobians(faceMats), fEl);
    std::copy(facedV.begin(), facedV.end(), dV.begin() + offset);
  }
  //metric and inverses
  metricTensor.resize(nJacs);
  inverseMetricTensor.resize(nJacs);
  invJacobians.resize(nJacs);
  for(int ip = 0 ; ip < nJacs; ip++){
    metricTensor[ip] = jacobians[ip]*(jacobians[ip].transpose());
    inverseMetricTensor[ip] = metricTensor[ip].inverse();
    invJacobians[ip] = jacobians[ip].transpose()*inverseMetricTensor[ip];
  }
  //deriv metric and hessians
  derivMetricTensor.resize(nIPsEl, EMatrix::Zero(embeddingDim, refDim*refDim));
  int leng = std::pow(refDim, 2);
  int lenJac = refDim*embeddingDim;
  std::fill(derivMetricTensor.begin(), derivMetricTensor.end(), EMatrix::Zero(refDim, leng));
  hessians.resize(nIPsEl, std::vector<EMatrix>(embeddingDim, EMatrix::Zero(refDim, refDim)));
  std::fill(hessians.begin(), hessians.end(), std::vector<EMatrix>(embeddingDim, EMatrix::Zero(refDim, refDim)));
  const std::vector< std::vector< std::vector<double> > > * derivShapes = refEl->getIPDerivShapeFunctions();
  for(int ip = 0; ip < nIPsEl; ip++){
    const std::vector< std::vector<double> > * derivShape = &(derivShapes->at(ip));
    for(int iN = 0; iN < nNodesEl; iN++){
      derivMetricTensor[ip] += EMap<const EVector>(derivShape->at(iN).data(), refDim)*(EMap<const EVector>(fieldMap["Metric"]->data() + iN*leng, leng));
      EMap<const EMatrix> jac(fieldMap["Jacobian"]->data() + iN*lenJac, refDim, embeddingDim);
      for(int iD = 0; iD < embeddingDim; iD++){
        hessians[ip][iD] += EMap<const EVector>(derivShape->at(iN).data(), refDim)*(jac.col(iD).transpose());
      }
    }
  }
  //Christoffel Symbols
  //normals
};//computeElementGeometry


};//hfox
