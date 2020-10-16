#include "HDGUNabU.h"

namespace hfox{

HDGUNabU::HDGUNabU(const ReferenceElement * re) : HDGNonLinearOperator(re){

};//constructor

HDGUNabU::~HDGUNabU(){

};//destructor

void HDGUNabU::allocate(int nDOFsPerNodeUser){
  HDGNonLinearOperator::allocate(nDOFsPerNodeUser);
};//allocate

void HDGUNabU::setFromBase(const std::vector<EVector> * ns){
  if(ns->at(0).size() != refEl->getDimension()){
    throw(ErrorHandle("HDGUNabU", "setFromBase", "the dimension of the normals must be equal to the dimension of the reference space"));
  }
  if(ns->size() != refEl->getNumFaces() * (refEl->getFaceElement()->getNumIPs())){
    throw(ErrorHandle("HDGUNabU", "setFromBase", "the number of normals should be the number of faces times the number of integration points per face"));
  }
  normals = ns;
};//setFromBase

void HDGUNabU::setSolution(const std::vector<EVector> & sol){
  int spaceDim = sol[0].size();
  std::fill(solNodes.begin(), solNodes.end(), EVector::Zero(spaceDim));
  solNodes.resize(refEl->getNumNodes(), EVector::Zero(spaceDim));
  for(int i = 0; i < refEl->getNumNodes(); i++){
    solNodes[i] = sol[i];
  }
  std::fill(sols.begin(), sols.end(), EVector::Zero(spaceDim));
  sols.resize(refEl->getNumIPs(), EVector::Zero(spaceDim));
  const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
  EMatrix solMat(spaceDim, refEl->getNumNodes());
  for(int i = 0; i < refEl->getNumNodes(); i++){
    solMat.col(i) = sol[i];
  }
  const std::vector<double> * pshape;
  for(int i = 0; i < refEl->getNumIPs(); i++){
    pshape = &(ipShapes->at(i));
    EMap<const EVector> shapes(pshape->data(), pshape->size());
    sols[i] = solMat*shapes;
  }
  int nFaces = refEl->getNumFaces();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int nfIPs = fEl->getNumIPs();
  int nNodesPFc = fEl->getNumNodes();
  const std::vector< std::vector<int> > * faceNodeMap = refEl->getFaceNodes();
  const std::vector< std::vector<double> > * fipShapes = fEl->getIPShapeFunctions();
  std::fill(faceSols.begin(), faceSols.end(), EVector::Zero(spaceDim));
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

void HDGUNabU::setTrace(const std::vector<EVector> & trace){
  int dim = trace[0].size();
  traceNodes = trace;
  int nFaces = refEl->getNumFaces();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int nNodesFc = fEl->getNumNodes();
  int nIPsFc = fEl->getNumIPs();
  const std::vector< std::vector<double> > * ipShapes = fEl->getIPShapeFunctions();
  const std::vector<double> * shapes;
  traces.resize(nFaces*nIPsFc, EVector::Zero(dim));
  std::fill(traces.begin(), traces.end(), EVector::Zero(dim));
  EMatrix traceFace(dim, nNodesFc);
  for(int iFace = 0; iFace < nFaces; iFace++){
    for(int iN = 0; iN < nNodesFc; iN++){
      traceFace.col(iN) = trace[iFace*nNodesFc + iN];
    }
    for(int ip = 0; ip < nIPsFc; ip++){
      shapes = &(ipShapes->at(ip));
      traces[iFace*nIPsFc + ip] = traceFace*EMap<const EVector>(shapes->data(), nNodesFc);
    }
  }
};//setTrace

void HDGUNabU::assemble(const std::vector<double> & dV, const std::vector<EMatrix> & invJacobians){
  if(!allocated){
    throw(ErrorHandle("HDGUNabU", "assemble", "the operator should be allocated before being assembled."));
  }
  if(sols.size() == 0){
    throw(ErrorHandle("HDGUNabU", "assemble", "the solution must be set before assembling"));
  }
  if(traces.size() == 0){
    throw(ErrorHandle("HDGUNabU", "assemble", "the trace must be set before assembling"));
  }
  if(normals == NULL){
    throw(ErrorHandle("HDGUNabU", "assemble", "must set the normals from the base before assembling"));
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
  //face integrals
  double faceSoldotNdV;
  int foffset;
  for(int iFace = 0; iFace < nFaces; iFace++){
    faceNodes = &(faceNodeMap->at(iFace));
    for(int ip = 0; ip < nIPsFc; ip++){
      shapes = &(fipShapes->at(ip));
      foffset = iFace*nIPsFc + ip;
      //faceSoldotNdV = (faceSols[foffset].transpose()*normals->at(foffset))(0,0)*dV[nIPsEl + foffset];
      faceSoldotNdV = (traces[foffset].transpose()*(normals->at(foffset)))(0,0)*dV[nIPsEl + foffset];
      for(int iN = 0; iN < nNodesPFc; iN++){
        for(int ndof = 0; ndof < nDOFsPerNode; ndof++){
          for(int jN = 0; jN < nNodesPFc; jN++){
            //op(startL + (iFace*nNodesPFc + iN)*nDOFsPerNode + ndof, faceNodes->at(jN)*nDOFsPerNode + ndof) +=
              //faceSoldotNdV*shapes->at(jN)*shapes->at(iN);
            //op(faceNodes->at(iN)*nDOFsPerNode + ndof, faceNodes->at(jN)*nDOFsPerNode + ndof) +=
              //faceSoldotNdV*shapes->at(jN)*shapes->at(iN);
            //op(faceNodes->at(iN)*nDOFsPerNode + ndof, startL + (iFace*nNodesPFc + jN)*nDOFsPerNode + ndof) +=
              //faceSoldotNdV*shapes->at(jN)*shapes->at(iN);
            op(startL + (iFace*nNodesPFc + iN)*nDOFsPerNode + ndof, startL + (iFace*nNodesPFc + jN)*nDOFsPerNode + ndof) +=
              faceSoldotNdV*(shapes->at(jN))*(shapes->at(iN));
            for(int mdof = 0; mdof < nDOFsPerNode; mdof++){
              //op(startL + (iFace*nNodesPFc + iN)*nDOFsPerNode + ndof, faceNodes->at(jN)*nDOFsPerNode + mdof) +=
                //faceSols[foffset][ndof]*shapes->at(jN)*shapes->at(iN)*normals->at(foffset)[mdof]*dV[nIPsEl + foffset];
              //op(faceNodes->at(iN)*nDOFsPerNode + ndof, faceNodes->at(jN)*nDOFsPerNode + mdof) +=
                //faceSols[foffset][ndof]*shapes->at(jN)*shapes->at(iN)*normals->at(foffset)[mdof]*dV[nIPsEl + foffset];
              //op(faceNodes->at(iN)*nDOFsPerNode + ndof, startL + (iFace*nNodesPFc + jN)*nDOFsPerNode + mdof) +=
                //traces[foffset][ndof]*shapes->at(jN)*shapes->at(iN)*normals->at(foffset)[mdof]*dV[nIPsEl + foffset];
              op(startL + (iFace*nNodesPFc + iN)*nDOFsPerNode + ndof, startL + (iFace*nNodesPFc + jN)*nDOFsPerNode + mdof) +=
                traces[foffset][ndof]*shapes->at(jN)*shapes->at(iN)*normals->at(foffset)[mdof]*dV[nIPsEl + foffset];
            }
          }
        }
      }
    }
  }
  //sum face contribution into first equations
  //for(int iFace = 0; iFace < nFaces; iFace++){
    //faceNodes = &(faceNodeMap->at(iFace));
    //for(int iN = 0; iN < nNodesPFc; iN++){
      //op.block(faceNodes->at(iN)*nDOFsPerNode, 0, nDOFsPerNode, lenU) += op.block(startL + (iFace*nNodesPFc + iN)*nDOFsPerNode, 0, nDOFsPerNode, lenU);
    //}
  //}
  for(int iFace = 0; iFace < nFaces; iFace++){
    faceNodes = &(faceNodeMap->at(iFace));
    for(int iN = 0; iN < nNodesPFc; iN++){
      op.block(faceNodes->at(iN)*nDOFsPerNode, startL, nDOFsPerNode, lenL) += op.block(startL + (iFace*nNodesPFc + iN)*nDOFsPerNode, startL, nDOFsPerNode, lenL);
    }
  }
  //bulk integrals
  EMatrix gradMat(spaceDim, nNodesEl);
  EMatrix gradSol(spaceDim, nDOFsPerNode);
  double divSol;
  for(int ip = 0; ip < nIPsEl; ip++){
    shapes = &(ipShapes->at(ip));
    derivShapes = &(ipDerivShapes->at(ip));
    gradSol = EMatrix::Zero(spaceDim, nDOFsPerNode);
    for(int iN = 0; iN < nNodesEl; iN++){
      gradMat.col(iN) = invJacobians[ip]*EMap<const EVector>(derivShapes->at(iN).data(), spaceDim);
      gradSol += gradMat.col(iN)*solNodes[iN].transpose();
    }
    divSol = gradSol.trace();
    for(int iN = 0; iN < nNodesEl; iN++){
      for(int ndof = 0; ndof < nDOFsPerNode; ndof++){
        for(int jN = 0; jN < nNodesEl; jN++){
          op(iN*nDOFsPerNode + ndof, jN*nDOFsPerNode + ndof) -= (divSol * shapes->at(iN) + 
            (sols[ip].transpose()*gradMat.col(iN))(0, 0))*dV[ip]*shapes->at(jN);
          //op(iN*nDOFsPerNode + ndof, jN*nDOFsPerNode + ndof) -= 
            //(sols[ip].transpose()*gradMat.col(iN))(0, 0)*dV[ip]*shapes->at(jN);
          for(int mdof = 0; mdof < nDOFsPerNode; mdof++){
            op(iN*nDOFsPerNode + ndof, jN*nDOFsPerNode + mdof) -= 
              sols[ip][ndof]*(gradMat(mdof, jN)*shapes->at(iN) + gradMat(mdof, iN)*shapes->at(jN))*dV[ip];
            //op(iN*nDOFsPerNode + ndof, jN*nDOFsPerNode + mdof) -= 
              //sols[ip][ndof]*(gradMat(mdof, iN)*shapes->at(jN))*dV[ip];
          }
        }
      }
    }
  }
  EVector solDiv2(lenU);
  EVector traceDiv2(lenL);
  for(int iN = 0 ; iN < nNodesEl; iN++){
    solDiv2.segment(iN*nDOFsPerNode, nDOFsPerNode) = solNodes[iN]/2.0;
  }
  for(int iN = 0 ; iN < nFaces*nNodesPFc; iN++){
    traceDiv2.segment(iN*nDOFsPerNode, nDOFsPerNode) = traceNodes[iN]/2.0;
  }
  rhs = EVector::Zero(rhs.size());
  rhs.segment(0, lenU) += op.block(0, 0, lenU, lenU)*solDiv2;
  rhs.segment(0, lenU) += op.block(0, startL, lenU, lenL)*traceDiv2;
  rhs.segment(startL, lenL) += op.block(startL, 0, lenL, lenU)*solDiv2;
  rhs.segment(startL, lenL) += op.block(startL, startL, lenL, lenL)*traceDiv2;
};//assemble

};//hfox
