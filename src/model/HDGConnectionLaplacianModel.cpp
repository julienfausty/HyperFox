#include "HDGConnectionLaplacianModel.h"

namespace hfox{


void HDGConnectionLaplacianModel::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator itfm;
  itfm = fm->find("DiffusionTensor");
  if(itfm != fm->end()){
    fieldMap["DiffusionTensor"] = &(itfm->second);
  }
  HDGEmbeddedModel::setFieldMap(fm);
};//setFieldMap

void HDGConnectionLaplacianModel::allocate(int nDOFsPerNode){
  assembly.matrix = Add;
  assembly.rhs = Add;
  HDGEmbeddedModel::allocate(nDOFsPerNode);
};//allocate

void HDGConnectionLaplacianModel::initializeOperators(){
  if(operatorMap.find("Source") == operatorMap.end()){
    operatorMap["Source"] = new Source(refEl);
  }
  HDGEmbeddedModel::initializeOperators();
};//initializeOperators

void HDGConnectionLaplacianModel::setSourceFunction(std::function<double(const std::vector<double>&)> s){
  if(!allocated){
    throw(ErrorHandle("HDGConnectionLaplacianModel", "setSourceFunction", "the model must be allocated before setting the source function"));
  }
  ((Source*)operatorMap["Source"])->setSourceFunction(s);
  sourceSet = 1;
};//setSourceFunction

void HDGConnectionLaplacianModel::parseDiffusionVals(){
  const std::vector<double> * diffVals = fieldMap.at("DiffusionTensor");
  int nNodes = refEl->getNumNodes();
  int dimDiffVal = diffVals->size()/nNodes;
  int dimMat = std::sqrt(dimDiffVal);
  std::vector<EMatrix> res(nNodes, EMatrix::Identity(embeddingDim, embeddingDim));
  if(dimMat == 1){
    for(int i = 0; i < nNodes; i++){
      res[i] *= (diffVals->at(i));
    }
  } else if(dimMat == embeddingDim){
    for(int i = 0; i < nNodes; i++){
      res[i] = EMap<const EMatrix>(diffVals->data() + i*dimDiffVal, dimMat, dimMat);
    }
  } else{
    throw(ErrorHandle("HDGConnectionLaplacianModel", "parseDiffusionVals", "the dimension of the diffusion tensor vales are not correct, they should be either scalar or tensor of the dimension of the embedding space"));
  }
  //get values at integration points
  int nIPsEl = refEl->getNumIPs();
  int nFaces = refEl->getNumFaces();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int nNodesFc = fEl->getNumNodes();
  int nIPsFc = fEl->getNumIPs();
  const std::vector< std::vector<double> > * shapes = refEl->getIPShapeFunctions();
  const std::vector< std::vector<int> > * faceIndexes = refEl->getFaceNodes();
  const std::vector< std::vector<double> > * fcShapes = fEl->getIPShapeFunctions();
  diffusionTensorVals.resize(nIPsEl + nFaces*nIPsFc);
  std::fill(diffusionTensorVals.begin(), diffusionTensorVals.end(), EMatrix::Zero(embeddingDim, embeddingDim));
  for(int ip = 0; ip < nIPsEl; ip++){
    const std::vector<double> * shape = &(shapes->at(ip));
    for(int iN = 0; iN < nNodes; iN++){
      diffusionTensorVals[ip] += shape->at(iN) * res[iN];
    }
  }
  for(int iFace = 0; iFace < nFaces; iFace++){
    for(int ip = 0; ip < nIPsFc; ip++){
      const std::vector<double> * fcShape = &(fcShapes->at(ip));
      for(int iN = 0; iN < nNodesFc; iN++){
        diffusionTensorVals[nIPsEl + iFace*nIPsFc + ip] += fcShape->at(iN)*res[faceIndexes->at(iFace)[iN]];
      }
    }
  }
};//parseDiffusionVals

void HDGConnectionLaplacianModel::computeLocalMatrix(){
  if(!(nodeSet and fieldSet)){
    throw(ErrorHandle("HDGConnectionLaplacianModel", "computeLocalMatrix", "the nodes and the fields should be set before computing"));
  }
  //set Matrix to null
  localMatrix = EMatrix::Zero(localMatrix.rows(), localMatrix.cols());
  computeElementGeometry();
  int refDim = refEl->getDimension();
  int nNodesEl = refEl->getNumNodes();
  int nIPsEl = refEl->getNumIPs();
  int nFaces = refEl->getNumFaces();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int nNodesFc = fEl->getNumNodes();
  int nIPsFc = fEl->getNumIPs();
  const std::vector< std::vector<double> > * elShapes = refEl->getIPShapeFunctions();
  const std::vector< std::vector< std::vector<double> > > * elDerivShapes = refEl->getIPDerivShapeFunctions();
  const std::vector< std::vector<double> > * fcShapes = fEl->getIPShapeFunctions();
  const std::vector< std::vector< std::vector<double> > > * fcDerivShapes = fEl->getIPDerivShapeFunctions();
  int lenU = nDOFsPNode * nNodesEl;
  int lenQ = lenU * embeddingDim;
  int lenL = nDOFsPNode * nFaces * nNodesFc;
  EMatrix tempMat;
  EVector tempVec;
  double tempVal, tempVal1;
  //bulk integrals
  for(int ip = 0; ip < nIPsEl; ip++){
    const std::vector<double> * elShape = &(elShapes->at(ip));
    const std::vector< std::vector<double> > * elDerivShape = &(elDerivShapes->at(ip));
    //Diffusive part
    //Squ
    for(int iN = 0; iN < nNodesEl; iN++){
      for(int jN = 0; jN < nNodesEl; jN++){
        tempVec = dV[ip]*elShape->at(iN)*diffusionTensorVals[ip]*invJacobians[ip]*EMap<const EVector>(elDerivShape->at(jN).data(), refDim);
        for(int iD = 0; iD < embeddingDim; iD++){
          for(int idof = 0; idof < nDOFsPNode; idof++){
            localMatrix(jN*nDOFsPNode + idof, lenU + (iN*embeddingDim + iD)*nDOFsPNode + idof) -= tempVec[iD];
          }
        }
      }
    }
    //gradient part
    int Qdim = embeddingDim*nDOFsPNode;
    for(int iN = 0; iN < nNodesEl; iN++){
      tempVal = dV[ip]*elShape->at(iN);
      for(int jN = 0; jN < nNodesEl; jN++){
        //Sqq
        localMatrix.block(lenU + jN*Qdim, lenU + iN*Qdim, Qdim, Qdim) = dV[ip]*elShape->at(iN)*elShape->at(jN)*EMatrix::Identity(Qdim, Qdim);
        for(int iD = 0; iD < embeddingDim; iD++){
          //here figure out the best assembly for curvature
          //tempVal1 = tempVal*invJacobian[ip].col(iD).dot(EMap<const EVector>(elDerivShape->at(jN).data(), refDim)*(elShape->at(ip)));
          for(int idof = 0; idof < nDOFsPNode; idof++){
            //localMatrix(lenU + jN*Qdim + iD*nDOFsPNode + idof, iN*nDOFsPNode + idof) -= tempVal;
          }
        }
      }
    }
  }
};//computeLocalMatrix

void HDGConnectionLaplacianModel::computeLocalRHS(){

};//computeLocalRHS

};
