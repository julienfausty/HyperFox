#include "L2InnerProduct.h"

namespace hfox{

void L2InnerProduct::checkField(Field * pField){
  if(pField == NULL){
    throw(ErrorHandle("L2InnerProduct", "checkField", "pointer to field cannot be null."));
  }
  FieldType ftype = *(pField->getFieldType());
  if(ftype == None or ftype == Edge){
    throw(ErrorHandle("L2InnerProduct", "checkField", "field types None and Edge are not yet supported."));
  }
  const ReferenceElement * refEl = pmesh->getReferenceElement();
  int nFEnts = *(pField->getNumEntities());
  int nFObj = *(pField->getNumObjPerEnt());
  if(ftype == Cell){
    if(nFEnts != pmesh->getNumberCells() or nFObj != refEl->getNumNodes()){
      throw(ErrorHandle("L2InnerProduct", "checkField", "a cell field must have as many entities as there are cells and as many objects as there are nodes per cell to be integrated."));
    }
  } else if(ftype == Node){
    if(nFEnts != pmesh->getNumberPoints() or nFObj != 1){
      throw(ErrorHandle("L2InnerProduct", "checkField", "a node field must have as many entities as there are nodes and 1 object per entity to be integrated."));
    }
  } else if(ftype == Face){
    if(nFEnts != pmesh->getNumberFaces() or nFObj != refEl->getFaceElement()->getNumNodes()){
      throw(ErrorHandle("L2InnerProduct", "checkField", "a face field must have as many entities as there are faces and as many objects as there are nodes per face to be integrated."));
    }
  }
};//checkField

void L2InnerProduct::setProduct(Field * lfield, Field * rfield){
  checkField(lfield);
  checkField(rfield);
  FieldType lftype = *(lfield->getFieldType());
  FieldType rftype = *(rfield->getFieldType());
  bool lCond = false;
  bool rCond = false;
  if(lftype != rftype){
    lCond = (lftype == Cell and rftype == Node);
    rCond = (rftype == Cell and lftype == Node);
    if(not (lCond or rCond)){
      throw(ErrorHandle("L2InnerProduct", "setProduct", "FieldField: Fields should either be of same type or a Cell-Node combination."));
    }
  }
  if(lftype == Node and rftype == Node){
    setType(Cell);
  } else if(lCond or rCond){
    setType(Cell);
  } else{
    setType(lftype);
  }
  int ldim = *(lfield->getNumValsPerObj());
  int rdim = *(rfield->getNumValsPerObj());
  if(ldim != rdim){
    throw(ErrorHandle("L2InnerProduct", "setProduct", "FieldField: Fields should have the same number of values per object."));
  }
  setDim(1);
  vectorDim = ldim;
  pmode = FieldField;
  fieldFuncs.resize(0);
  fields.resize(0);
  fields.push_back(lfield);
  fields.push_back(rfield);
};//setProduct(FieldField)

void L2InnerProduct::setProduct(Field * lfield, FieldFunc rfield){
  checkField(lfield);
  FieldType lftype = *(lfield->getFieldType());
  if(lftype == Node){
    setType(Cell);
  } else{
    setType(lftype);
  }
  int ldim = *(lfield->getNumValsPerObj());
  std::vector<double> dbuff(pmesh->getNodeSpaceDimension(), 1.0);
  std::vector<double> resbuff;
  rfield(dbuff, &resbuff);
  int rdim = resbuff.size();
  if(ldim != rdim){
    throw(ErrorHandle("L2InnerProduct", "setProduct", "FieldFunction: objects should have the same number of values per node."));
  }
  setDim(1);
  vectorDim = ldim;
  pmode = FieldFunction;
  fields.resize(0);
  fields.push_back(lfield);
  fieldFuncs.resize(0);
  fieldFuncs.push_back(rfield);
};//setProduct(FieldFunction)

void L2InnerProduct::setProduct(FieldFunc lfield, Field * rfield){
  setProduct(rfield, lfield);
};//setProduct(FieldFunction)

void L2InnerProduct::setProduct(FieldFunc lfield, FieldFunc rfield){
  setType(Cell);
  std::vector<double> dbuff(pmesh->getNodeSpaceDimension(), 1.0);
  std::vector<double> resbuff;
  lfield(dbuff, &resbuff);
  int ldim = resbuff.size();
  resbuff.resize(0);
  rfield(dbuff, &resbuff);
  int rdim = resbuff.size();
  if(ldim != rdim){
    throw(ErrorHandle("L2InnerProduct", "setProduct", "FunctionFunction: objects should have the same number of values per node."));
  }
  setDim(1);
  vectorDim = ldim;
  pmode = FunctionFunction;
  fields.resize(0);
  fieldFuncs.resize(0);
  fieldFuncs.push_back(lfield);
  fieldFuncs.push_back(rfield);
};//setProduct(FunctionFunction)

void L2InnerProduct::evaluateIntegrand(int iEnt, std::vector<double> * ipVals){
  const ReferenceElement * refEl = pmesh->getReferenceElement();
  Partitioner * part = pmesh->getPartitioner();
  if(type == Face){
    refEl = refEl->getFaceElement();
  }
  int dimSpace = pmesh->getNodeSpaceDimension();
  int nNodesEnt = refEl->getNumNodes();
  int nIPs = refEl->getNumIPs();
  const std::vector< std::vector<double> > * phis = refEl->getIPShapeFunctions();
  EMatrix phiMat(nIPs*vectorDim, nNodesEnt*vectorDim);
  phiMat = EMatrix::Zero(nIPs*vectorDim, nNodesEnt*vectorDim);
  for(int ip = 0; ip < nIPs; ip++){
    const std::vector<double> * phiAtIP = &(phis->at(ip));
    for(int iN = 0; iN < nNodesEnt; iN++){
      double val = phiAtIP->at(iN);
      for(int k = 0; k < vectorDim; k++){
        phiMat(ip*vectorDim + k, iN*vectorDim + k) = val;
      }
    }
  }
  ipVals->resize(nIPs*dimIntegral, 0.0);
  std::fill(ipVals->begin(), ipVals->end(), 0.0);
  std::vector<double> valNodes(nNodesEnt * vectorDim, 0.0);
  std::vector<std::vector<double> > valIPs(2, std::vector<double>(nIPs * vectorDim, 0.0));
  std::vector<double> elCoords(dimSpace * nNodesEnt, 0.0);
  std::vector<double> ipCoords(dimSpace * nIPs, 0.0);
  std::vector<int> cell;
  std::vector<double> dbuff(vectorDim, 0.0);
  std::vector<double> valbuff(vectorDim, 0.0);
  if(type != Face){
    pmesh->getCell(iEnt, &cell);
  } else {
    pmesh->getFace(iEnt, &cell);
  }
  for(int i = 0; i < fields.size(); i++){
    std::fill(valNodes.begin(), valNodes.end(), 0.0);
    FieldType ftype = *(fields[i]->getFieldType());
    if(ftype == Node){
      for(int iN = 0; iN < nNodesEnt; iN++){
        int locN = cell[iN];
        if(part != NULL){
          locN = part->global2LocalNode(cell[iN]);
        }
        if(locN != -1){
          fields[i]->getValues(locN, &dbuff);
        } else{
          fields[i]->getParValues(cell[iN], &dbuff);
        }
        std::copy(dbuff.begin(), dbuff.end(), valNodes.begin() + iN*vectorDim);
      }
    } else if(ftype == Cell or ftype == Face){
      fields[i]->getValues(iEnt, &(valNodes));
    }
    EMap<EVector>(valIPs[i].data(), valIPs[i].size()) = phiMat*EMap<EVector>(valNodes.data(), valNodes.size());
  }
  if(pmode != FieldField){
    phiMat = EMatrix::Zero(nIPs * dimSpace, nNodesEnt * dimSpace);
    for(int ip = 0; ip < nIPs; ip++){
      const std::vector<double> * phiAtIP = &(phis->at(ip));
      for(int iN = 0; iN < nNodesEnt; iN++){
        double val = phiAtIP->at(iN);
        for(int k = 0; k < dimSpace; k++){
          phiMat(ip*dimSpace + k, iN*dimSpace + k) = val;
        }
      }
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
  }
  for(int i = 0; i < fieldFuncs.size(); i++){
    for(int iN = 0; iN < nIPs; iN++){
      dbuff.resize(dimSpace, 0.0);
      std::copy(ipCoords.begin() + iN*dimSpace, ipCoords.begin() + (iN+1)*dimSpace, dbuff.begin());
      fieldFuncs[i](dbuff, &valbuff);
      std::copy(valbuff.begin(), valbuff.end(), valIPs[i + fields.size()].begin() + iN*vectorDim);
    }
  }
  for(int ip = 0; ip < nIPs; ip++){
    ipVals->at(ip) = EMap<EVector>(valIPs[0].data() + ip*vectorDim, vectorDim).dot(EMap<EVector>(valIPs[1].data() + ip*vectorDim, vectorDim));
  }
};//evaluateIntegrand

};
