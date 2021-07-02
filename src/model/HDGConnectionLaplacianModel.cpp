#include "HDGConnectionLaplacianModel.h"

namespace hfox{


void HDGConnectionLaplacianModel::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator itfm;
  itfm = fm->find("DiffusionTensor");
  if(itfm != fm->end()){
    fieldMap["DiffusionTensor"] = &(itfm->second);
  }
  if(timeScheme != NULL){
    timeScheme->setFieldMap(fm);
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
  parseTauValues();
  int refDim = refEl->getDimension();
  int nNodesEl = refEl->getNumNodes();
  int nIPsEl = refEl->getNumIPs();
  int nFaces = refEl->getNumFaces();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int nNodesFc = fEl->getNumNodes();
  int nIPsFc = fEl->getNumIPs();
  if(fieldMap.find("DiffusionTensor") != fieldMap.end()){
    parseDiffusionVals();
  } else {
    diffusionTensorVals.resize(nIPsEl + nFaces*nIPsFc, EMatrix::Identity(embeddingDim, embeddingDim));
    std::fill(diffusionTensorVals.begin(), diffusionTensorVals.end(), EMatrix::Identity(embeddingDim, embeddingDim));
  }
  const std::vector< std::vector<double> > * elShapes = refEl->getIPShapeFunctions();
  const std::vector< std::vector< std::vector<double> > > * elDerivShapes = refEl->getIPDerivShapeFunctions();
  const std::vector< std::vector<int> > * faceIndexes = refEl->getFaceNodes();
  const std::vector< std::vector<double> > * fcShapes = fEl->getIPShapeFunctions();
  int lenU = nDOFsPNode * nNodesEl;
  int lenQ = lenU * embeddingDim;
  int lenL = nDOFsPNode * nFaces * nNodesFc;
  int Qdim = embeddingDim*nDOFsPNode;
  EMatrix tempMat;
  EVector tempVec;
  double tempVal, tempVal1, tempVal2, tempVal3;
  //bulk integrals
  for(int ip = 0; ip < nIPsEl; ip++){
    const std::vector<double> * elShape = &(elShapes->at(ip));
    const std::vector< std::vector<double> > * elDerivShape = &(elDerivShapes->at(ip));
    //Diffusive part
    //Suq
    for(int iN = 0; iN < nNodesEl; iN++){
      tempMat = dV[ip]*(elShape->at(iN))*diffusionTensorVals[ip]*invJacobians[ip];
      for(int jN = 0; jN < nNodesEl; jN++){
        tempVec = tempMat*EMap<const EVector>(elDerivShape->at(jN).data(), refDim);
        for(int iD = 0; iD < embeddingDim; iD++){
          for(int idof = 0; idof < nDOFsPNode; idof++){
            localMatrix(jN*nDOFsPNode + idof, lenU + (iN*embeddingDim + iD)*nDOFsPNode + idof) += tempVec[iD];
          }
        }
      }
    }
    //gradient part
    for(int iN = 0; iN < nNodesEl; iN++){
      tempVal = dV[ip]*elShape->at(iN);
      for(int jN = 0; jN < nNodesEl; jN++){
        //Sqq
        localMatrix.block(lenU + jN*Qdim, lenU + iN*Qdim, Qdim, Qdim) += tempVal*elShape->at(jN)*EMatrix::Identity(Qdim, Qdim);
        //Squ
        for(int iD = 0; iD < embeddingDim; iD++){
          for(int kD = 0; kD < refDim; kD++){
            for(int lD = 0; lD < refDim; lD++){
              tempVal1 = 0;
              for(int rD = 0; rD < refDim; rD++){
                tempVal1 += christoffelSymbols[ip][rD](kD, lD)*jacobians[ip](rD, iD);
              }
              tempVal2 = tempVal*inverseMetricTensor[ip](kD, lD)*(jacobians[ip](lD, iD)*(elDerivShape->at(jN)[kD]) + (hessians[ip][iD](kD, lD) - tempVal1)*(elShape->at(jN)));
              for(int idof = 0; idof < nDOFsPNode; idof++){
                localMatrix(lenU + jN*Qdim + iD*nDOFsPNode + idof, iN*nDOFsPNode + idof) += tempVal2;
              }
            }
          }
        }
      }
    }
  }
  //face integrals
  for(int iFace = 0; iFace < nFaces; iFace++){
    const std::vector<int> * locFaceIndexes = &(faceIndexes->at(iFace));
    for(int ip = 0; ip < nIPsFc; ip++){
      const std::vector<double> * fcShape = &(fcShapes->at(ip));
      for(int iN = 0; iN < nNodesFc; iN++){
        tempVal = dV[nIPsEl + iFace*nIPsFc + ip]*(fcShape->at(iN));
        for(int jN = 0; jN < nNodesFc; jN++){
          tempVal1 = tempVal * (fcShape->at(jN));
          //tau terms
          for(int idof = 0; idof < nDOFsPNode; idof++){
            for(int jdof = 0; jdof < nDOFsPNode; jdof++){
              tempVal2 = taus[iFace*nIPsFc + ip](jdof, idof)*tempVal1;
              //Sll
              localMatrix(lenU + lenQ + (iFace*nNodesFc + jN)*nDOFsPNode + jdof, lenU + lenQ + (iFace*nNodesFc + iN)*nDOFsPNode + idof) -= tempVal2;
              //Slu
              localMatrix(lenU + lenQ + (iFace * nNodesFc + jN)*nDOFsPNode + jdof, (locFaceIndexes->at(iN))*nDOFsPNode + idof) += tempVal2;
              //Sul
              localMatrix((locFaceIndexes->at(jN))*nDOFsPNode + jdof, lenU + lenQ + (iFace* nNodesFc + iN)*nDOFsPNode + idof) -= tempVal2;
              //Suu
              localMatrix((locFaceIndexes->at(jN))*nDOFsPNode + jdof, (locFaceIndexes->at(iN))*nDOFsPNode + idof) += tempVal2;
            }
          }
          //q face terms
          for(int iD = 0; iD < embeddingDim; iD++){
            tempVal2 = tempVal1 * diffusionTensorVals[nIPsEl + iFace*nIPsFc + ip].col(iD).dot(normals[iFace*nIPsFc + ip]);
            tempVal3 = tempVal1 * normals[iFace*nIPsFc + ip][iD];
            for(int idof = 0; idof < nDOFsPNode; idof++){
              //Slq
              localMatrix(lenU + lenQ + (iFace * nNodesFc + jN)*nDOFsPNode + idof, lenU + Qdim*(locFaceIndexes->at(iN)) + iD*nDOFsPNode + idof) -= tempVal2;
              //Suq
              localMatrix(locFaceIndexes->at(jN)*nDOFsPNode + idof, lenU + Qdim*(locFaceIndexes->at(iN)) + iD*nDOFsPNode + idof) -= tempVal2;
              //Sql
              localMatrix(lenU + locFaceIndexes->at(jN)*Qdim + iD*nDOFsPNode + idof, lenU + lenQ + (iFace*nNodesFc + iN)*nDOFsPNode + idof) -= tempVal3;
            }
          }
        }
      }
    }
  }
};//computeLocalMatrix

void HDGConnectionLaplacianModel::computeLocalRHS(){
  if(!(nodeSet and fieldSet)){
    throw("HDGConnectionLaplacianModel", "computeLocalRHS", "the nodes and the fields should be set before computing");
  }
  localRHS = EVector::Zero(localRHS.size());
  if(sourceSet){
    ((Source*)operatorMap["Source"])->calcSource(*elementNodes);
    operatorMap["Source"]->assemble(dV, invJacobians);
    localRHS.segment(0, operatorMap["Source"]->getMatrix()->rows()) = (EVector) *(operatorMap["Source"]->getMatrix());
  }
};//computeLocalRHS

};
