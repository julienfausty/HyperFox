#include "HDGBase.h"

namespace hfox{

HDGBase::HDGBase(const ReferenceElement * re) : HDGOperator(re){

};//constructor

HDGBase::~HDGBase(){

};//destructor

void HDGBase::allocate(int nDOFsPerNode){

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
    EMap<const EMatrix> tau(userTaus.data() + i*nNodesFace*sizeTau, sizeTau, nNodesFace);
    for(int j = 0; j < fEl->getNumIPs(); j++){
      EMap<const EVector> shape((ipShapes->at(j)).data(), (ipShapes->at(j)).size());
      EMap<EMatrix>(taus.data() + (i*(fEl->getNumIPs()) + j)*sizeTau, sizeTau, 1) = tau * shape;
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
  int startQblock = nNodesEl*nDOFsPerNode;
  int startLblock = nNodesEl*(dim+1)*nDOFsPerNode;
  int lenUblock = nNodesEl*nDOFsPerNode;
  int lenQblock = nNodesEl*dim*nDOFsPerNode;
  int lenLblock = nFaces*nNodesFace*nDOFsPerNode;
  int sizeTau = nDOFsPerNode * nDOFsPerNode;
  const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
  const std::vector< std::vector< std::vector<double> > > * ipDerivShapes = refEl->getIPDerivShapeFunctions();
  const std::vector< std::vector<int> > * faceNodeMap = refEl->getFaceNodes();
  const std::vector< std::vector<double> > * fipShapes = fEl->getIPShapeFunctions();
  const std::vector<int> * faceNodes;
  const std::vector<double> * shapes;
  const std::vector< std::vector<double> > * derivShapes;
  int offset;
  double val;
  EMatrix matMeasure(nDOFsPerNode, nDOFsPerNode);
  EVector vecMeasure(dim);
  EVector buffVec(dim);
  EMatrix shapeShape(nNodesFace, nNodesFace);
  //start with face integrals
  for(int ip = 0; ip < nIPsFace; ip++){
    shapes = &(fipShapes->at(ip));
    shapeShape = EMap<const EVector>(shapes->data(), nNodesFace)*EMap<const EVector>(shapes->data(), nNodesFace).transpose();
    for(int iFace = 0; iFace < nFaces; iFace++){
      offset = iFace*nIPsFace + ip;
      faceNodes = &(faceNodeMap->at(iFace));
      matMeasure = dV[nIPsEl + offset]*EMap<EMatrix>(taus.data() + offset*sizeTau, nDOFsPerNode, nDOFsPerNode);
      vecMeasure = dV[nIPsEl + offset]*(normals[offset]);
      for(int iN = 0; iN < nNodesFace; iN++){
        for(int ndof = 0; ndof < nDOFsPerNode; ndof++){
          for(int jN = 0; jN < nNodesFace; jN++){
            for(int mdof = 0; mdof < nDOFsPerNode; mdof++){
              val = matMeasure(ndof, mdof) * shapeShape(iN, jN);
              //Sll
              op(startLblock + (iFace*nNodesFace + iN)*nDOFsPerNode + ndof, startLblock + (iFace*nNodesFace + jN)*nDOFsPerNode + mdof) -= val;
              //Slu
              op(startLblock + (iFace*nNodesFace + iN)*nDOFsPerNode + ndof, faceNodes->at(jN)*nDOFsPerNode + mdof) += val;
              //Suu
              op(faceNodes->at(iN)*nDOFsPerNode + ndof, faceNodes->at(jN)*nDOFsPerNode + mdof) += val;
              //Sul
              op(faceNodes->at(iN)*nDOFsPerNode + ndof, startLblock + (iFace*nNodesFace + jN)*nDOFsPerNode + mdof) -= val;
            }
            buffVec = vecMeasure*shapeShape(iN, jN);
            for(int d = 0; d < dim; d++){
              //Sql
              op(startQblock + (faceNodes->at(iN)*dim + d)*nDOFsPerNode + ndof, startLblock + (iFace*nNodesFace + jN)*nDOFsPerNode + ndof) -= buffVec[d];
            }
          }
        }
      }
    }
  }
  //do bulk integrals
  shapeShape.resize(nNodesEl, nNodesEl);
  for(int ip = 0; ip < nIPsEl; ip++){
    shapes = &(ipShapes->at(ip));
    derivShapes = &(ipDerivShapes->at(ip));
    shapeShape = EMap<const EVector>(shapes->data(), nNodesEl)*EMap<const EVector>(shapes->data(), nNodesEl).transpose();
    for(int iN = 0; iN < nNodesEl; iN++){
      vecMeasure = (invJacobians[ip]*EMap<const EVector>(derivShapes->at(iN).data(), dim))*dV[ip];
      for(int d = 0; d < dim; d++){
        for(int ndof = 0; ndof < nDOFsPerNode; ndof++){
          for(int jN = 0; jN < nNodesEl; jN++){
            //Squ
            op(startQblock + (iN*dim + d)*nDOFsPerNode + ndof, jN*nDOFsPerNode + ndof) += vecMeasure[d]*shapes->at(jN);
            //Sqq
            op(startQblock + (iN*dim + d)*nDOFsPerNode + ndof, startQblock + (jN*dim + d)* nDOFsPerNode + ndof) += dV[ip]*shapeShape(iN, jN);
          }
        }
      }
    }
  }
};//assemble

}//hfox
