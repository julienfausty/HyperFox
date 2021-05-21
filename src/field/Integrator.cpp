#include "Integrator.h"

namespace hfox{

void Integrator::integrate(std::vector<double> * integral){
  if(type == None){
    throw(ErrorHandle("Integrator", "integrate", "cannot integrate a None type field."));
  }
  if(pmesh == NULL){
    throw(ErrorHandle("Integrator", "integrate", "must set the mesh before trying to integrate"));
  }
  if(dimIntegral <= 0){
    throw(ErrorHandle("Integrator", "integrate", "the dimension of the integral must be superior to zero in order to integrate."));
  }
  integral->resize(dimIntegral);
  std::fill(integral->begin(), integral->end(), 0.0);
  const ReferenceElement * refEl = pmesh->getReferenceElement();
  Partitioner * part = pmesh->getPartitioner();
  if(type == Face){
    refEl = refEl->getFaceElement();
  }
  int dimSpace = pmesh->getNodeSpaceDimension();
  int dimEl = refEl->getDimension();
  int nNodesEnt = refEl->getNumNodes();
  int nIPs = refEl->getNumIPs();
  const std::vector< std::vector<double> > * phis = refEl->getIPShapeFunctions();
  EMatrix phiMat(nIPs*dimIntegral, nNodesEnt*dimIntegral);
  for(int ip = 0; ip < nIPs; ip++){
    const std::vector<double> * phiAtIP = &(phis->at(ip));
    for(int iN = 0; iN < nNodesEnt; iN++){
      double val = phiAtIP->at(iN);
      for(int k = 0; k < dimIntegral; k++){
        phiMat(ip*dimIntegral + k, iN*dimIntegral + k) = val;
      }
    }
  }
  std::vector<int> cell;
  std::vector<double> valsIPs(nIPs*dimIntegral, 0.0);
  std::vector< std::vector<double> > elCoords(nNodesEnt, std::vector<double>(dimSpace, 0.0));
  std::vector<EMatrix> elJacs(nIPs, EMatrix::Zero(dimSpace, dimEl));
  std::vector<double> measure(nIPs, 0.0);
  std::vector<double> coordBuff(dimSpace, 0.0);
  EMap<EVector> integralVec(integral->data(), dimIntegral);
  if(type == Cell or type == Face){
    cell.resize(nNodesEnt);
    int nEnts = pmesh->getNumberCells();
    if(type == Face){
      nEnts = pmesh->getNumberFaces();
    }
    for(int iEl = 0; iEl < nEnts; iEl++){
      if(type == Cell or type == Node){
        pmesh->getCell(iEl, &cell);
      } else{
        pmesh->getFace(iEl, &cell);
      }
      for(int iN = 0; iN < nNodesEnt; iN++){
        int locN = cell[iN];
        if(part != NULL){
          locN = part->global2LocalNode(cell[iN]);
        }
        if(locN != -1){
          pmesh->getPoint(locN, &coordBuff);
        } else{
          pmesh->getGhostPoint(cell[iN], &coordBuff);
        }
        std::copy(coordBuff.begin(), coordBuff.end(), elCoords[iN].begin());
      }
      elJacs = Operator::calcJacobians(elCoords, refEl);
      measure = Operator::calcMeasure(Operator::calcDetJacobians(elJacs), refEl);
      evaluateIntegrand(iEl, &valsIPs);
      for(int ip = 0; ip < nIPs; ip++){
        integralVec += measure[ip]*EMap<EVector>(valsIPs.data() + ip*dimIntegral, dimIntegral);
      }
    }
  } else if(type == Edge or type == Node){
    throw(ErrorHandle("Integrator", "integrate", "edge and node entity types are not supported yet."));
  }
  if(part != NULL){
    MPI_Allreduce(MPI_IN_PLACE, integral->data(), dimIntegral, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
};//integrate

};
