#include "FieldIntegrator.h"

namespace hfox{

void FieldIntegrator::setField(Field * myField){
  pField = myField;
  if(pField == NULL){
    throw(ErrorHandle("FieldIntegrator", "setField", "pointer to field cannot be null."));
  }
  FieldType ftype = *(pField->getFieldType());
  if(ftype == Node or ftype == Cell){
    setType(Cell);
  } else if(ftype == Face){
    setType(Face);
  } else{
    throw(ErrorHandle("FieldIntegrator", "setField", "field types None and Edge are not yet supported."));
  }
  dimIntegral = *(pField->getNumValsPerObj());
  const ReferenceElement * refEl = pmesh->getReferenceElement();
  int nFEnts = *(pField->getNumEntities());
  int nFObj = *(pField->getNumObjPerEnt());
  if(ftype == Cell){
    if(nFEnts != pmesh->getNumberCells() or nFObj != refEl->getNumNodes()){
      throw(ErrorHandle("FieldIntegrator", "setField", "a cell field must have as many entities as there are cells and as many objects as there are nodes per cell to be integrated."));
    }
  } else if(ftype == Node){
    if(nFEnts != pmesh->getNumberPoints() or nFObj != 1){
      throw(ErrorHandle("FieldIntegrator", "setField", "a node field must have as many entities as there are nodes and 1 object per entity to be integrated."));
    }
  } else if(ftype == Face){
    if(nFEnts != pmesh->getNumberFaces() or nFObj != refEl->getFaceElement()->getNumNodes()){
      throw(ErrorHandle("FieldIntegrator", "setField", "a face field must have as many entities as there are faces and as many objects as there are nodes per face to be integrated."));
    }
  }
};//setField

void FieldIntegrator::evaluateIntegrand(int iEnt, std::vector<double> * ipVals){
  FieldType ftype = *(pField->getFieldType());
  if(type == Cell and ftype == Face){
    throw(ErrorHandle("FieldIntegrator", "evaluateIntegrand", "cannot integrate a face field over a cell"));
  } else if(type == Face and ftype == Cell){
    throw(ErrorHandle("FieldIntegrator", "evaluateIntegrand", "integration of cell fields over faces is not yet supported."));
  }
  const ReferenceElement * refEl = pmesh->getReferenceElement();
  Partitioner * part = pmesh->getPartitioner();
  if(type == Face){
    refEl = refEl->getFaceElement();
  }
  int nNodesEnt = refEl->getNumNodes();
  int nIPs = refEl->getNumIPs();
  const std::vector< std::vector<double> > * phis = refEl->getIPShapeFunctions();
  EMatrix phiMat(nIPs*dimIntegral, nNodesEnt*dimIntegral);
  phiMat = EMatrix::Zero(nIPs*dimIntegral, nNodesEnt*dimIntegral);
  for(int ip = 0; ip < nIPs; ip++){
    const std::vector<double> * phiAtIP = &(phis->at(ip));
    for(int iN = 0; iN < nNodesEnt; iN++){
      double val = phiAtIP->at(iN);
      for(int k = 0; k < dimIntegral; k++){
        phiMat(ip*dimIntegral + k, iN*dimIntegral + k) = val;
      }
    }
  }
  std::vector<double> valsNodes(nNodesEnt*dimIntegral, 0.0);
  ipVals->resize(nIPs*dimIntegral, 0.0);
  std::fill(ipVals->begin(), ipVals->end(), 0.0);
  std::vector<int> cell;
  std::vector<double> dbuff(dimIntegral, 0.0);
  if(ftype == Node){
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
        pField->getValues(locN, &dbuff);
      } else{
        pField->getParValues(cell[iN], &dbuff);
      }
      std::copy(dbuff.begin(), dbuff.end(), valsNodes.begin() + iN*dimIntegral);
    }
  } else if(ftype == Cell or ftype == Face){
    pField->getValues(iEnt, &valsNodes);
  }
  EMap<EVector>(ipVals->data(), ipVals->size()) = phiMat * EMap<EVector>(valsNodes.data(), valsNodes.size());
};//evaluateIntegrand

};
