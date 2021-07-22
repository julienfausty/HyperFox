#include "IP2NodeSolver.h"

namespace hfox{

void IP2NodeSolver::setField(Field * pField){
  if(*(pField->getFieldType()) != Cell){
    throw(ErrorHandle("IP2NodeSolver", "setField", "field must be a Cell field."));
  }
  myField = pField;
};//setField

void IP2NodeSolver::setFunction(std::function< void(std::vector<double>&, std::vector<double>*) > ipFunc){

  ipFuncNodes = ipFunc;
  ipFuncJacs = nullptr;
  ipFuncCell = nullptr;

};//setFunction

void IP2NodeSolver::setFunction(std::function< void(std::vector<EMatrix>&, std::vector<double>*) > ipFunc){

  ipFuncJacs = ipFunc;
  ipFuncNodes = nullptr;
  ipFuncCell = nullptr;

};//setFunction

void IP2NodeSolver::setFunction(std::function< void(int, std::vector<double>*) > ipFunc){

  ipFuncCell = ipFunc;
  ipFuncNodes = nullptr;
  ipFuncJacs = nullptr;

};//setFunction

void IP2NodeSolver::solve(){
  //Checks
  if(myField == NULL){
    throw(ErrorHandle("IP2NodeSolver", "solve", "must set field before attempting a solve."));
  }
  if((ipFuncCell == nullptr) and (ipFuncNodes == nullptr) and (ipFuncJacs == nullptr)){
    throw(ErrorHandle("IP2NodeSolver", "solve", "must set function before attempting a solve."));
  }
  const Partitioner * part = myMesh->getPartitioner();
  const ReferenceElement * refEl = myMesh->getReferenceElement();
  int nDOFsField = *(myField->getNumValsPerObj());
  int ndim = myMesh->getNodeSpaceDimension();
  int nNodesEl = refEl->getNumNodes();
  int nIPs = refEl->getNumIPs();
  int refDim = refEl->getDimension();
  std::vector<int> cell;
  std::vector<double> node(ndim, 0.0);
  std::vector< std::vector<double> > nodes(nNodesEl, std::vector<double>(ndim, 0.0));
  std::vector<double> ipNodes(nIPs*ndim, 0.0);
  std::vector<double> vals(nIPs*nDOFsField, 0.0);
  //integration helpers
  const std::vector< std::vector<double> > * shapes = refEl->getIPShapeFunctions();
  const std::vector< std::vector< std::vector<double> > > * derivShapes = refEl->getIPDerivShapeFunctions();
  EMatrix phiMat(nNodesEl, nIPs);
  for(int ip = 0; ip < nIPs; ip++){
    phiMat.col(ip) = EMap<const EVector>(shapes->at(ip).data(), nNodesEl);
  }
  EMatrix nodeMat(ndim, nNodesEl);
  std::vector<double> measure(nIPs, 0.0);
  std::vector<EMatrix> jacobians(nIPs, EMatrix::Zero(refDim, ndim));
  std::vector<EMatrix> invJacobians(nIPs, EMatrix::Zero(ndim, refDim));
  EVector rhs = EVector::Zero(nNodesEl*nDOFsField);
  myMass.allocate(nDOFsField);
  Eigen::HouseholderQR<EMatrix> invMass(nNodesEl * nDOFsField, nNodesEl * nDOFsField);
  for(int iEl = 0; iEl < myMesh->getNumberCells(); iEl++){
    myMesh->getCell(iEl, &cell);
    for(int iN = 0; iN < nNodesEl; iN++){
      int locN = cell[iN];
      if(part != NULL){
        locN = part->global2LocalNode(cell[iN]);
      }
      if(locN != -1){
        myMesh->getPoint(locN, &node);
      } else{
        myMesh->getGhostPoint(cell[iN], &node);
      }
      nodes[iN] = node;
      nodeMat.col(iN) = EMap<EVector>(node.data(), ndim);
    }
    jacobians = Operator::calcJacobians(nodes, refEl);
    invJacobians = Operator::calcInvJacobians(jacobians);
    measure = Operator::calcMeasure(Operator::calcDetJacobians(jacobians), refEl);
    myMass.assemble(measure, invJacobians);
    invMass.compute(*(myMass.getMatrix()));
    if(ipFuncNodes != nullptr){
      EMap<EMatrix>(ipNodes.data(), ndim, nIPs) = nodeMat*phiMat;
      ipFuncNodes(ipNodes, &vals);
    } else if(ipFuncJacs != nullptr){
      ipFuncJacs(jacobians, &vals);
    } else {
      ipFuncCell(iEl, &vals);
    }
    rhs = EVector::Zero(nNodesEl*nDOFsField);
    for(int ip = 0; ip < nIPs; ip++){
      EMap<EVector> val(vals.data() + nDOFsField*ip, nDOFsField);
      for(int iN = 0; iN < nNodesEl; iN++){
        rhs.segment(iN*nDOFsField, nDOFsField) += measure[ip] * phiMat(iN, ip) * val;
      }
    }
    EMap<EVector>(myField->getValues()->data() + iEl*nNodesEl*nDOFsField, nNodesEl*nDOFsField) = invMass.solve(rhs);
  }
};//solve

void IP2NodeSolver::jacobianVals(std::vector<EMatrix> & jacs, std::vector<double> * vals){
  vals->resize(jacs.size() * jacs[0].size(), 0.0);
  for(int i = 0; i < jacs.size(); i++){
    std::copy(jacs[i].data(), jacs[i].data() + jacs[i].size(), vals->begin() + i*jacs[i].size());
  }
};//jacobiansVals

void IP2NodeSolver::metricVals(std::vector<EMatrix> & jacs, std::vector<double> * vals){
  int refDim = jacs[0].rows();
  int leng = std::pow(refDim, 2);
  vals->resize(jacs.size() * leng, 0.0);
  EMatrix g = EMatrix::Zero(refDim, refDim);
  for(int i = 0; i < jacs.size(); i++){
    g = jacs[i]*jacs[i].transpose();
    std::copy(g.data(), g.data() + leng, vals->begin() + i*leng);
  }
};//metricVals

};//hfox
