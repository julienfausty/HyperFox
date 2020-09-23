#include "UNabU.h"

namespace hfox{

void UNabU::setSolution(const std::vector<EVector> & sol){
  int spaceDim = sol[0].size();
  solNodes.resize(refEl->getNumNodes(), EVector::Zero(spaceDim));
  for(int i = 0; i < refEl->getNumNodes(); i++){
    solNodes[i] = sol[i];
  }
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
};//setSolution

void UNabU::assemble(const std::vector<double> & dV, const std::vector<EMatrix> & invJacobians){
  if(!allocated){
    throw(ErrorHandle("UNabU", "assemble", "the operator should be allocated before being assembled."));
  }
  if(sols.size() == 0){
    throw(ErrorHandle("UNabU", "assemble", "the solution must be set before assembling"));
  }
  op = EMatrix::Zero(op.rows(), op.cols());
  rhs = EVector::Zero(op.cols());
  int dim = refEl->getDimension();
  int nNodesEl = refEl->getNumNodes();
  int nIPsEl = refEl->getNumIPs();
  const std::vector< std::vector<double> > * ipShapes = refEl->getIPShapeFunctions();
  const std::vector< std::vector< std::vector<double> > > * ipDerivShapes = refEl->getIPDerivShapeFunctions();
  const std::vector<double> * shapes;
  const std::vector< std::vector<double> > * derivShapes;
  EMatrix gradSol(dim, nDOFsPerNode);
  EVector preGradOp(nDOFsPerNode);
  double buff;
  for(int ip = 0; ip < nIPsEl; ip++){
    shapes = &(ipShapes->at(ip));
    derivShapes = &(ipDerivShapes->at(ip));
    gradSol = EMatrix::Zero(dim, nDOFsPerNode);
    for(int iN =0; iN < nNodesEl; iN++){
      gradSol += invJacobians[ip]*EMap<const EVector>(derivShapes->at(iN).data(), dim)*solNodes[iN].transpose();
    }
    gradSol *= dV[ip];
    preGradOp = sols[ip]*dV[ip];
    for(int iN = 0; iN < nNodesEl; iN++){
      buff = (preGradOp.transpose()*invJacobians[ip]*EMap<const EVector>(derivShapes->at(iN).data(), dim))(0,0);
      for(int ndof = 0; ndof < nDOFsPerNode; ndof++){
        for(int jN = 0; jN < nNodesEl; jN++){
          op(jN*nDOFsPerNode + ndof, iN*nDOFsPerNode + ndof) += buff*(shapes->at(jN));
          for(int mdof = 0; mdof < nDOFsPerNode; mdof++){
            op(iN*nDOFsPerNode + ndof, jN*nDOFsPerNode + mdof) += gradSol(mdof, ndof)*shapes->at(iN)*shapes->at(jN);
          }
        }
      }
    }
  }
  EVector solDiv2(nDOFsPerNode*nNodesEl);
  for(int iN = 0; iN < nNodesEl; iN++){
    solDiv2.segment(iN*nDOFsPerNode, nDOFsPerNode) = solNodes[iN]/2.0;
  }
  rhs = op*solDiv2;
}

};//hfox
