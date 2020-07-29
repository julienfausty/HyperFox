#include "HDGBase.h"

namespace hfox{

HDGBase::HDGBase(const ReferenceElement * re) : HDGOperator(re){
  bulkMass = new Mass(re);
  faceMass = new Mass(re->getFaceElement());
  convection = new Convection(re);
};//constructor

HDGBase::~HDGBase(){
  delete bulkMass;
  delete faceMass;
  delete convection;
};//destructor

void HDGBase::allocate(int nDOFsPerNode){
  bulkMass->allocate(nDOFsPerNode*(refEl->getDimension()));
  convection->allocate(nDOFsPerNode);
  faceMass->allocate(nDOFsPerNode);
  HDGOperator::allocate(nDOFsPerNode);
};//allocate

void HDGBase::setTau(const std::vector<double> & userTaus){
  int sizeTau = nDOFsPerNode * nDOFsPerNode;
  int nFaces = refEl->getNumFaces();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int nNodesFace = fEl->getNumNodes();
  taus.resize(nFaces*(fEl->getNumIPs())*sizeTau, 0.0);
  const std::vector< std::vector<double> > * ipShapes = fEl->getIPShapeFunctions();
  for(int i = 0; i < nFaces; i++){
    EMap<const EMatrix> tau(userTaus.data() + i*nNodesFace*sizeTau, nNodesFace, sizeTau);
    for(int j = 0; j < fEl->getNumIPs(); j++){
      EMap<const EVector> shape((ipShapes->at(j)).data(), (ipShapes->at(j)).size());
      EMap<EMatrix>(taus.data() + (i*(fEl->getNumIPs()) + j)*sizeTau, 1, sizeTau) = shape.transpose() * tau;
    }
  }
};//setTau

void HDGBase::calcNormals(const std::vector< std::vector<double> > & nodes, const std::vector<EMatrix> & jacobians){
  int offset = refEl->getNumIPs();
  int nFaces = refEl->getNumFaces();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int nIPsFace = fEl->getNumIPs();
  normals.resize(nFaces*nIPsFace, EVector::Zero(nodes[0].size()));
  const std::vector< std::vector<int> > * faceNodeMap = refEl->getFaceNodes();
  int v0, vn;
  std::vector<int>::const_iterator itv;
  EVector testVec(nodes[0].size());
  for(int i = 0; i < nFaces; i++){
    v0 = faceNodeMap->at(i)[0];
    for(int k = 0; k < refEl->getNumNodes(); k++){
      itv = std::find((faceNodeMap->at(i)).begin(), (faceNodeMap->at(i)).end(), k);
      if(itv == (faceNodeMap->at(i)).end()){
        vn = k;
        break;
      }
    }
    testVec = EMap<const EVector>(nodes[vn].data(), nodes[vn].size()) - EMap<const EVector>(nodes[v0].data(), nodes[v0].size());
    for(int j = 0; j < nIPsFace; j++){
      int indexN = i*nIPsFace + j;
      EMatrix temp = jacobians[offset + indexN].fullPivLu().kernel();
      normals[indexN] = temp;
      normals[indexN].normalize();
      double prod = testVec.dot(normals[i*nIPsFace+j]);
      if(prod > 0){
        normals[indexN] *= -1;
      }
    }
  }
};//calcNormals

void HDGBase::assemble(const std::vector< double > & dV, const std::vector< EMatrix > & invJacobians){
  if(!allocated){
    throw(ErrorHandle("HDGBase", "assemble", "the operator should be allocated before being assembled."));
  }
  if(taus.size() == 0){
    throw(ErrorHandle("HDGBase", "assemble", "the stabilization parameter (tau) must be set before assembling"));
  }
  if(normals.size() == 0){
    throw(ErrorHandle("HDGBase", "assemble", "the normals must be calculated before assembling"));
  }
  op = EMatrix::Zero(op.rows(), op.cols());
  const ReferenceElement * fEl = refEl->getFaceElement();
  int nNodesEl = refEl->getNumNodes();
  int nIPsEl = refEl->getNumIPs();
  int dim = refEl->getDimension();
  int nFaces = refEl->getNumFaces();
  int nNodesFace = fEl->getNumNodes();
  int nIPsFace = fEl->getNumIPs();
  int startUblock = 0;
  int startQblock = nNodesEl;
  int startLblock = nNodesEl*(dim+1);
  int lenUblock = nNodesEl;
  int lenQblock = nNodesEl*dim;
  int lenLblock = nFaces*nNodesFace;
  int sizeTau = nDOFsPerNode * nDOFsPerNode;
  //Sqq
  bulkMass->assemble(dV, invJacobians);
  op.block(startQblock, startQblock, lenQblock, lenQblock) += *(bulkMass->getMatrix());
  //Squ
  for(int i = 0; i < dim; i++){
    std::vector<EVector> velocity(nNodesEl, EMatrix::Identity(dim, dim).col(i));
    convection->setVelocity(velocity);
    convection->assemble(dV, invJacobians);
    for(int j = 0; j < nNodesEl; j++){
      op.block(startQblock + j*dim + i, startUblock, 1, lenUblock) += (convection->getMatrix()->col(j)).transpose();
    }
  }
  std::vector<double> normaldV(nIPsFace, 0.0);
  std::vector<EMatrix> faceInvJacs(nIPsFace, EMatrix::Zero(dim, dim-1));
  std::vector<double> faceTaus(nIPsFace, 0.0);
  std::vector<EVector> faceNormals(nIPsFace, EVector::Zero(dim));
  const std::vector<int> * nodeMap;
  for(int iFace = 0; iFace < nFaces; iFace++){
    std::copy(invJacobians.begin() + nIPsEl + iFace*nIPsFace, invJacobians.begin() + nIPsEl + (iFace+1)*nIPsFace, faceInvJacs.begin());
    nodeMap = &(refEl->getFaceNodes()->at(iFace));
    int startDiag = startLblock + iFace*nNodesFace;
    //Sql
    int findex = iFace*nIPsFace;
    for(int i = 0; i < dim; i++){
      for(int j = 0; j < nIPsFace; j++){
        normaldV[j] = dV[nIPsEl + findex + j]*(normals[findex + j][i]);
      }
      faceMass->assemble(normaldV, faceInvJacs);
      for(int j = 0; j < nNodesFace; j++){
        op.block(startQblock + dim*(nodeMap->at(j)) + i, startDiag, 1, nNodesFace) -= (faceMass->getMatrix()->col(j)).transpose();
      }
    }
  }

  if(nDOFsPerNode > 1){
    multiplyDOFs();
  }
  double val;
  bool allZeros;
  for(int dofUp = 0; dofUp < nDOFsPerNode; dofUp++){
    for(int dofDown = 0; dofDown < nDOFsPerNode; dofDown++){
      for(int iFace = 0; iFace < nFaces; iFace++){
        std::copy(invJacobians.begin() + nIPsEl + iFace*nIPsFace, invJacobians.begin() + nIPsEl + (iFace+1)*nIPsFace, faceInvJacs.begin());
        nodeMap = &(refEl->getFaceNodes()->at(iFace));
        for(int i = 0; i < nIPsFace; i++){
          faceTaus[i] = taus[(iFace*nIPsFace + i)*sizeTau + dofUp*nDOFsPerNode + dofDown];
        }
        allZeros = 1;
        for(int i = 0; i < nIPsFace; i++){
          if(faceTaus[i] != 0){
            allZeros = 0;
            break;
          }
        }
        if(allZeros){
          continue;
        }
        for(int i = 0; i < nIPsFace; i++){
          faceTaus[i] *= dV[nIPsEl + iFace * nIPsFace + i];
        } 
        faceMass->assemble(faceTaus, faceInvJacs);
        for(int k = 0; k < nNodesFace; k++){
          for(int l = 0; l < nNodesFace; l++){
            val = (*(faceMass->getMatrix()))(k,l);
            //Sll
            op((startLblock + iFace*nNodesFace + k)*nDOFsPerNode + dofUp, (startLblock + iFace*nNodesFace + l)*nDOFsPerNode + dofDown) -= val;
            //Slu
            op((startLblock + iFace*nNodesFace + k)*nDOFsPerNode + dofUp, nodeMap->at(l)*nDOFsPerNode + dofDown) += val;
            //Suu
            op(nodeMap->at(k)*nDOFsPerNode + dofUp, nodeMap->at(l)*nDOFsPerNode + dofDown) += val;
            //Sul
            op(nodeMap->at(k)*nDOFsPerNode + dofUp, (startLblock + iFace*nNodesFace + l)*nDOFsPerNode + dofDown) -= val;
          }
        }
        //op.block(startDiag, startDiag, nNodesFace, nNodesFace) -= *(faceMass->getMatrix());
        ////Slu
        //for(int k = 0; k < nNodesFace; k++){
          //op.block(startDiag, nodeMap->at(k), nNodesFace, 1) += (*(faceMass->getMatrix())).col(k);
        //}
      }
    }
  }
  //Suu
  //for(int iFace = 0; iFace < nFaces; iFace++){
    //nodeMap = &(refEl->getFaceNodes()->at(iFace));
    //int startDiag = startLblock + iFace*nNodesFace;
    //for(int i = 0; i < nNodesFace; i++){
      //op.block(nodeMap->at(i), 0, 1, lenUblock) += op.block(startDiag + i, 0, 1, lenUblock);
    //}
  //}
  ////Sul
  //op.block(startUblock, startLblock, lenUblock, lenLblock) -= op.block(startLblock, startUblock, lenLblock, lenUblock).transpose();
};//assemble

}//hfox
