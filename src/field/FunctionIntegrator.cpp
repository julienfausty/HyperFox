#include "FunctionIntegrator.h"

namespace hfox{

void FunctionIntegrator::setFunction(FieldFunc aFunc){
  setType(Cell);
  std::vector<double> dbuff(pmesh->getNodeSpaceDimension(), 1.0);
  std::vector<double> resbuff;
  aFunc(dbuff, &resbuff);
  setDim(resbuff.size());
  myFunc = aFunc;
};//setFunction

void FunctionIntegrator::evaluateIntegrand(int iEnt, std::vector<double> * ipVals){
  const ReferenceElement * refEl = pmesh->getReferenceElement();
  Partitioner * part = pmesh->getPartitioner();
  if(type == Face){
    refEl = refEl->getFaceElement();
  }
  int dimSpace = pmesh->getNodeSpaceDimension();
  int nNodesEnt = refEl->getNumNodes();
  int nIPs = refEl->getNumIPs();
  const std::vector< std::vector<double> > * phis = refEl->getIPShapeFunctions();
  EMatrix phiMat(nIPs*dimSpace, nNodesEnt*dimSpace);
  phiMat = EMatrix::Zero(nIPs*dimSpace, nNodesEnt*dimSpace);
  for(int ip = 0; ip < nIPs; ip++){
    const std::vector<double> * phiAtIP = &(phis->at(ip));
    for(int iN = 0; iN < nNodesEnt; iN++){
      double val = phiAtIP->at(iN);
      for(int k = 0; k < dimSpace; k++){
        phiMat(ip*dimSpace + k, iN*dimSpace + k) = val;
      }
    }
  }
  ipVals->resize(nIPs*dimIntegral, 0.0);
  std::fill(ipVals->begin(), ipVals->end(), 0.0);
  std::vector<double> elCoords(dimSpace * nNodesEnt, 0.0);
  std::vector<double> ipCoords(dimSpace * nIPs, 0.0);
  std::vector<int> cell;
  std::vector<double> dbuff(dimIntegral, 0.0);
  std::vector<double> valbuff(dimIntegral, 0.0);
  if(type != Face){
    pmesh->getCell(iEnt, &cell);
  } else {
    pmesh->getFace(iEnt, &cell);
  }
  for(int iN = 0; iN < nNodesEnt; iN++){
    int locN = cell[iN];
    if(part != NULL){
      locN = part->global2LocalNode(cell[iN]);
    }
    if(locN != -1){
      pmesh->getPoint(locN, &dbuff);
    } else{
      pmesh->getGhostPoint(cell[iN], &dbuff);
    }
    std::copy(dbuff.begin(), dbuff.end(), elCoords.begin() + iN*dimSpace);
  }
  EMap<EVector>(ipCoords.data(), ipCoords.size()) = phiMat *EMap<EVector>(elCoords.data(), elCoords.size());
  for(int iN = 0; iN < nIPs; iN++){
    dbuff.resize(dimSpace, 0.0);
    std::copy(ipCoords.begin() + iN*dimSpace, ipCoords.begin() + (iN+1)*dimSpace, dbuff.begin());
    myFunc(dbuff, &valbuff);
    std::copy(valbuff.begin(), valbuff.end(), ipVals->begin() + iN*dimIntegral);
  }
};//evaluateIntegrand

};//hfox
