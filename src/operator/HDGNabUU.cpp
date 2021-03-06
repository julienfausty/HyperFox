#include "HDGNabUU.h"

namespace hfox{

HDGNabUU::HDGNabUU(const ReferenceElement * re) : HDGNonLinearOperator(re){

};//constructor

HDGNabUU::~HDGNabUU(){

};//destructor

void HDGNabUU::allocate(int nDOFsPerNodeUser){
  HDGNonLinearOperator::allocate(nDOFsPerNodeUser);
};//allocate

void HDGNabUU::setFromBase(const std::vector<EVector> * ns){
  if(ns->at(0).size() != refEl->getDimension()){
    throw(ErrorHandle("HDGNabUU", "setFromBase", "the dimension of the normals must be equal to the dimension of the reference space"));
  }
  if(ns->size() != refEl->getNumFaces() * (refEl->getFaceElement()->getNumIPs())){
    throw(ErrorHandle("HDGNabUU", "setFromBase", "the number of normals should be the number of faces times the number of integration points per face"));
  }
  normals = ns;
};//setFromBase

void HDGNabUU::setSolution(const std::vector<EVector> & sol){
  int spaceDim = sol[0].size();
  solNodes.resize(refEl->getNumNodes(), EVector::Zero(spaceDim));
  for(int i = 0; i < refEl->getNumNodes(); i++){
    solNodes[i] = sol[i];
  }
  sols.resize(refEl->getNumIPs(), EVector::Zero(spaceDim));
  const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
  const std::vector< std::vector< std::vector<double> > > * ipDerivShapes = refEl->getIPDerivShapeFunctions();
  EMatrix solMat(spaceDim, refEl->getNumNodes());
  for(int i = 0; i < refEl->getNumNodes(); i++){
    solMat.col(i) = sol[i];
  }
  const std::vector<double> * pshape;
  const std::vector< std::vector<double> > * derivShape;
  for(int i = 0; i < refEl->getNumIPs(); i++){
    pshape = &(ipShapes->at(i));
    derivShape = &(ipDerivShapes->at(i));
    EMap<const EVector> shapes(pshape->data(), pshape->size());
    sols[i] = solMat*shapes;
  }
  int nFaces = refEl->getNumFaces();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int nfIPs = fEl->getNumIPs();
  int nNodesPFc = fEl->getNumNodes();
  const std::vector< std::vector<int> > * faceNodeMap = refEl->getFaceNodes();
  const std::vector< std::vector<double> > * fipShapes = fEl->getIPShapeFunctions();
  faceSols.resize(nFaces*nfIPs, EVector::Zero(spaceDim));
  EMatrix nodeSols(spaceDim, nNodesPFc);
  for(int iFace = 0; iFace < nFaces; iFace++){
    for(int i = 0; i < nNodesPFc; i++){
      nodeSols.col(i) = sol[faceNodeMap->at(iFace)[i]];
    }
    for(int i = 0; i < nfIPs; i++){
      EMap<const EVector> shape(fipShapes->at(i).data(), fipShapes->at(i).size());
      faceSols[iFace*nfIPs + i] = nodeSols*shape;
    }
  }
};//setSolution

void HDGNabUU::assemble(const std::vector<double> & dV, const std::vector<EMatrix> & invJacobians){
  if(!allocated){
    throw(ErrorHandle("HDGNabUU", "assemble", "the operator should be allocated before being assembled."));
  }
  if(sols.size() == 0){
    throw(ErrorHandle("HDGNabUU", "assemble", "the solution must be set before assembling"));
  }
  if(normals == NULL){
    throw(ErrorHandle("HDGNabUU", "assemble", "must set the normals from the base before assembling"));
  }
  op = EMatrix::Zero(op.rows(), op.cols());
  rhs = EVector::Zero(rhs.size());
  int spaceDim = refEl->getDimension();
  int nNodesEl = refEl->getNumNodes();
  int nIPsEl = refEl->getNumIPs();
  int nFaces = refEl->getNumFaces();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int nNodesPFc = fEl->getNumNodes();
  int nIPsFc = fEl->getNumIPs();
  const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
  const std::vector< std::vector< std::vector<double> > > * ipDerivShapes = refEl->getIPDerivShapeFunctions();
  const std::vector< std::vector<int> > * faceNodeMap = refEl->getFaceNodes();
  const std::vector< std::vector<double> > * fipShapes = fEl->getIPShapeFunctions();
  const std::vector<int> * faceNodes;
  const std::vector<double> * shapes;
  const std::vector< std::vector<double> > * derivShapes;
  int lenU = nNodesEl * nDOFsPerNode;
  int lenQ = spaceDim*lenU;
  int startL = lenU + lenQ;
  int lenL = nNodesPFc * nFaces * nDOFsPerNode;
  //start with face integrals
  EMatrix buffMat(spaceDim, spaceDim);
  EVector buffVec(spaceDim);
  EVector buffVec2(spaceDim);
  EVector buffVec3(spaceDim);
  double solNdV;
  double solShapedV;
  int offset;
  //face integrals
  for(int iFace = 0; iFace < nFaces; iFace++){
    faceNodes = &(faceNodeMap->at(iFace));
    for(int ip = 0; ip < nIPsFc; ip++){
      offset = iFace*nIPsFc + ip;
      shapes = &(fipShapes->at(ip));
      solNdV = (faceSols[offset].transpose()*(normals->at(offset)))(0, 0)*dV[nIPsEl + offset];
      for(int iN = 0; iN < nNodesPFc; iN++){
        for(int ndof = 0; ndof < nDOFsPerNode; ndof++){
          solShapedV = faceSols[offset][ndof]*(shapes->at(iN))*dV[nIPsEl + offset];
          for(int jN = 0; jN < nNodesPFc; jN++){
            op(startL + (iFace*nNodesPFc + iN)*nDOFsPerNode + ndof, (faceNodes->at(jN))*nDOFsPerNode + ndof) += solNdV * shapes->at(iN) * shapes->at(jN);
            for(int mdof = 0; mdof < nDOFsPerNode; mdof++){
              op(startL + (iFace*nNodesPFc + iN)*nDOFsPerNode + ndof, (faceNodes->at(jN))*nDOFsPerNode + mdof) += solShapedV * normals->at(offset)[mdof] * shapes->at(jN);
            }
          }
        }
      }
    }
  }
  //sum face integrals into first equations
  for(int iFace = 0; iFace < nFaces; iFace++){
    faceNodes = &(faceNodeMap->at(iFace));
    for(int iN = 0; iN < nNodesPFc; iN++){
      op.block(faceNodes->at(iN)*nDOFsPerNode, 0, nDOFsPerNode, lenU) += op.block(startL + (iFace*nNodesPFc + iN)*nDOFsPerNode, 0, nDOFsPerNode, lenU);
    }
  }
  //bulk integrals
  double solNabShapedV;
  for(int ip = 0; ip < nIPsEl; ip++){
    shapes = &(ipShapes->at(ip));
    derivShapes = &(ipDerivShapes->at(ip));
    for(int iN = 0; iN < nNodesEl; iN++){
      buffVec = (invJacobians[ip]*EMap<const EVector>(derivShapes->at(iN).data(), derivShapes->at(iN).size()))*dV[ip];
      solNabShapedV = (sols[ip].transpose()*buffVec)(0,0);
      for(int ndof = 0; ndof < nDOFsPerNode; ndof++){
        for(int jN = 0; jN < nNodesEl; jN++){
          op(iN*nDOFsPerNode + ndof, jN*nDOFsPerNode + ndof) -= solNabShapedV * shapes->at(jN);
          for(int mdof = 0; mdof < nDOFsPerNode; mdof++){
            op(iN*nDOFsPerNode + ndof, jN*nDOFsPerNode + mdof) -= buffVec[mdof]*sols[ip][ndof]*shapes->at(jN);
          }
        }
      }
    }
  }
  //remove the divergence components
  //double divSol0dV;
  //EMatrix derivShapesMat(spaceDim, nNodesEl);
  //for(int ip = 0; ip < nIPsEl; ip++){
    //shapes = &(ipShapes->at(ip));
    //derivShapes = &(ipDerivShapes->at(ip));
    //divSol0dV = 0;
    //derivShapesMat = EMatrix::Zero(spaceDim, nNodesEl);
    //for(int iN = 0; iN < nNodesEl; iN++){
      //derivShapesMat.col(iN) = invJacobians[ip]*EMap<const EVector>(derivShapes->at(iN).data(), derivShapes->at(iN).size());
      //divSol0dV += ((derivShapesMat.col(iN).transpose())*solNodes[iN])(0,0);
    //}
    //divSol0dV *= dV[ip];
    //for(int iN = 0; iN < nNodesEl; iN++){
      //buffVec = shapes->at(iN) * sols[ip] * dV[ip];
      //for(int ndof = 0; ndof < nDOFsPerNode; ndof++){
        //for(int jN = 0; jN < nNodesEl; jN++){
          //op(iN * nDOFsPerNode + ndof, jN*nDOFsPerNode + ndof) -= divSol0dV * shapes->at(iN) * shapes->at(jN);
          //for(int mdof = 0; mdof < nDOFsPerNode; mdof++){
            //op(iN * nDOFsPerNode + ndof, jN*nDOFsPerNode + mdof) -= buffVec[ndof] * derivShapesMat(mdof, jN);
          //}
        //}
      //}
    //}
  //}
  EVector solDiv2(lenU);
  for(int i = 0; i < nNodesEl; i++){
    solDiv2.segment(i*nDOFsPerNode, nDOFsPerNode) = solNodes[i]/2.0;
  }
  rhs.segment(0, lenU) += op.block(0, 0, lenU, lenU)*solDiv2;
  rhs.segment(startL, lenL) += op.block(startL, 0, lenL, lenU)*solDiv2;
};//assemble

};//hfox
