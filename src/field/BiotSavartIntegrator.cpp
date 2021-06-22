#include "BiotSavartIntegrator.h"

namespace hfox{

BiotSavartIntegrator::BiotSavartIntegrator(Mesh * myMesh, int dim) : Integrator(myMesh, dim){
  setType(Cell);
};//constructor

void BiotSavartIntegrator::setCurrent(Field * myCurrent){
  pCurrent = myCurrent;
};//setCurrent

void BiotSavartIntegrator::setCoordinates(std::vector<double> * coords){
  pCoords = coords;
};//setCoordinates

void BiotSavartIntegrator::setDim(int dim){
  if(dim != 3){
    throw(ErrorHandle("BiotSavartIntegrator", "setDim", "dimension of BiotSavartIntegrator must be set to 3"));
  }
  dimIntegral = 3;
};//setDim

void BiotSavartIntegrator::evaluateIntegrand(int iEnt, std::vector<double> * ipVals){
  FieldType ftype = *(pCurrent->getFieldType());
  if(type == Cell and ftype == Face){
    throw(ErrorHandle("BiotSavartIntegrator", "evaluateIntegrand", "cannot integrate a face field over a cell"));
  } else if(type == Face and ftype == Cell){
    throw(ErrorHandle("BiotSavartIntegrator", "evaluateIntegrand", "integration of cell fields over faces is not yet supported."));
  }
  int dimSpace = pmesh->getNodeSpaceDimension();
  int dimCurrent = *(pCurrent->getNumValsPerObj());
  const ReferenceElement * refEl = pmesh->getReferenceElement();
  Partitioner * part = pmesh->getPartitioner();
  if(type == Face){
    refEl = refEl->getFaceElement();
  }
  int nNodesEnt = refEl->getNumNodes();
  int nIPs = refEl->getNumIPs();
  const std::vector< std::vector<double> > * phis = refEl->getIPShapeFunctions();
  EMatrix phiMat(nIPs, nNodesEnt);
  phiMat = EMatrix::Zero(nIPs, nNodesEnt);
  for(int ip = 0; ip < nIPs; ip++){
    const std::vector<double> * phiAtIP = &(phis->at(ip));
    for(int iN = 0; iN < nNodesEnt; iN++){
      phiMat(ip, iN) = phiAtIP->at(iN);
    }
  }
  std::vector<double> valsNodes(nNodesEnt*dimCurrent, 0.0);
  std::vector<double> currentIPs(nIPs * dimCurrent, 0.0);
  ipVals->resize(nIPs*dimIntegral, 0.0);
  std::fill(ipVals->begin(), ipVals->end(), 0.0);
  EMatrix nodesEnt = EMatrix::Zero(dimSpace, nNodesEnt);
  EMatrix nodesIPs = EMatrix::Zero(dimSpace, nIPs);
  std::vector<int> cell;
  std::vector<double> nodeBuff(dimSpace, 0.0);
  std::vector<double> dbuff(dimCurrent, 0.0);
  Eigen::Vector3d j = EVector::Zero(3);
  Eigen::Vector3d r = EVector::Zero(3);
  Eigen::Vector3d X = EVector::Zero(3);
  double rnorm = 0.0;
  X.head(pCoords->size()) = EMap<const EVector>(pCoords->data(), pCoords->size());
  if(type == Cell){
    pmesh->getCell(iEnt, &cell);
  } else if(type == Face){
    pmesh->getFace(iEnt, &cell);
  }
  for(int iN = 0; iN < nNodesEnt; iN++){
    int locN = cell[iN];
    if(part != NULL){
      locN = part->global2LocalNode(cell[iN]);
    }
    if(locN != -1){
      pmesh->getPoint(locN, &nodeBuff);
      if(ftype == Node){
        pCurrent->getValues(locN, &dbuff);
      }
    } else{
      pmesh->getGhostPoint(cell[iN], &nodeBuff);
      if(ftype == Node){
        pCurrent->getParValues(cell[iN], &dbuff);
      }
    }
    nodesEnt.col(iN) = EMap<EVector>(nodeBuff.data(), dimSpace);
    if(ftype == Node){
      std::copy(dbuff.begin(), dbuff.end(), valsNodes.begin() + iN*dimCurrent);
    }
  }
  if(ftype == Cell or ftype == Face){
    pCurrent->getValues(iEnt, &valsNodes);
  }
  nodesIPs = nodesEnt*(phiMat.transpose());
  EMap<EMatrix>(currentIPs.data(), dimCurrent, nIPs) = EMap<EMatrix>(valsNodes.data(), dimCurrent, nNodesEnt)*(phiMat.transpose());
  for(int ip = 0; ip < nIPs; ip++){
    j = EVector::Zero(3);
    r = EVector::Zero(3);
    j.head(dimCurrent) = EMap<EVector>(currentIPs.data() + dimCurrent*ip, dimCurrent);
    r.head(dimSpace) = nodesIPs.col(ip);
    r = X - r;
    rnorm = r.norm();
    EMap<EVector>(ipVals->data() + ip*dimIntegral, dimIntegral) = (1.0/std::pow(rnorm, 3)) * j.cross(r);
  }
};//evaluateIntegrand

void BiotSavartIntegrator::integrate(std::vector<double> * integral){
  Integrator::integrate(integral);
  double mu0 = 1.25663706212*(1e-6);
  EMap<EVector>(integral->data(), integral->size()) *= mu0/(4.0*M_PI);
};//integrate

};//hfox
