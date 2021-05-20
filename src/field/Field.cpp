#include "Field.h"

namespace hfox{

Field::Field(){
  pmesh = NULL;
  type = None;
  numObjPerEnt = 0;
  numValsPerObj = 0;
};//empty constructor

Field::Field(Mesh * mesh){
  pmesh = mesh;
  type = None;
  numObjPerEnt = 0;
  numValsPerObj = 0;
};//empty constructor

Field::Field(Mesh * mesh, FieldType t, int nObjPerEnt, int nValsPerObj){
  pmesh = mesh;
  type = t;
  computeNumEntities();
  numObjPerEnt = nObjPerEnt;
  numValsPerObj = nValsPerObj;
  allocate();
};//allocation constructor

Field::~Field(){
};//Destructor

void Field::computeNumEntities(){
  switch(type){
    case None: {numEntities = 0; break;}
    case Node: {numEntities = pmesh->getNumberPoints(); break;}
    case Edge: {throw(ErrorHandle("Field", "allocate", "edge fields are not supported yet.")); break;}
    case Face: {numEntities = pmesh->getNumberFaces(); break;}
    case Cell: {numEntities = pmesh->getNumberCells(); break;}
  }
};//computeNumEntities

void Field::allocate(){
  if(pmesh != NULL){
    int length = numEntities*numObjPerEnt*numValsPerObj;
    values.resize(length, 0.0);
  }
};//allocate

void Field::getValues(int i, std::vector<double> * vals){
  int numValsPerEnt = numObjPerEnt*numValsPerObj;
  vals->resize(numValsPerEnt);
  std::copy(values.begin() + i*numValsPerEnt, values.begin() + (i+1)*numValsPerEnt, vals->begin());
}; //getValue

void Field::getSliceValues(std::vector<int> & is, std::vector<double> * vals){
  int numValsPerEnt = numObjPerEnt*numValsPerObj;
  vals->resize(is.size()*numValsPerEnt);
  for(int i = 0; i < is.size(); i++){
    std::copy(values.begin() + is[i]*numValsPerEnt, values.begin() + (is[i]+1)*numValsPerEnt, 
        vals->begin() + i*numValsPerEnt);
  }
};//getSliceValues

void Field::getParValues(int i, std::vector<double> * vals){
  int numValsPerEnt = numObjPerEnt*numValsPerObj;
  vals->resize(numValsPerEnt);
  std::vector<int>::iterator it = std::find(parIds.begin(), parIds.end(), i);
  if(it == parIds.end()){
    throw(ErrorHandle("Field", "getParValues", "could not find the global ID " + std::to_string(i) + " in the parallel ID list."));
  }
  int locInd = std::distance(parIds.begin(), it);
  std::copy(parValues.begin() + locInd*numValsPerEnt, parValues.begin() + (locInd+1)*numValsPerEnt, vals->begin());
}; //getParValue

void Field::integrate(std::vector<double> * integral){
  if(type == None){
    throw(ErrorHandle("Field", "integrate", "cannot integrate a None type field."));
  }
  integral->resize(numValsPerObj);
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
  EMatrix phiMat(nIPs*numValsPerObj, nNodesEnt*numValsPerObj);
  for(int ip = 0; ip < nIPs; ip++){
    const std::vector<double> * phiAtIP = &(phis->at(ip));
    for(int iN = 0; iN < nNodesEnt; iN++){
      double val = phiAtIP->at(iN);
      for(int k = 0; k < numValsPerObj; k++){
        phiMat(ip*numValsPerObj + k, iN*numValsPerObj + k) = val;
      }
    }
  }
  std::vector<int> cell;
  std::vector<double> valsNodes(nNodesEnt*numValsPerObj, 0.0);
  std::vector<double> valsIPs(nIPs*numValsPerObj, 0.0);
  std::vector<double> elCoords(nNodesEnt*dimSpace, 0.0);
  std::vector<EMatrix> elJacs(nIPs, EMatrix::Zero(dimSpace, dimEl));
  std::vector<double> measure(nIPs, 0.0);
  std::vector<double> coordBuff(dimSpace, 0.0);
  std::vector<double> dbuff(numValsPerObj, 0.0);
  EMap<EVector> integralVec(integral->data(), numValsPerObj);
  if(type == Cell or type == Node or type == Face){
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
          if(type == Node){
            this->getValues(locN, &dbuff);
          }
        } else{
          pmesh->getGhostPoint(cell[iN], &coordBuff);
          if(type == Node){
            this->getParValues(cell[iN], &dbuff);
          }
        }
        std::copy(coordBuff.begin(), coordBuff.end(), elCoords.begin() + iN*dimSpace);
        if(type == Node){
          std::copy(dbuff.begin(), dbuff.end(), valsNodes.begin() + iN*numValsPerObj);
        }
      }
      if(type == Cell or type == Face){
        this->getValues(iEl, &valsNodes);
      }
      elJacs = Operator::calcJacobians(elCoords, refEl);
      measure = Operator::calcMeasure(Operator::calcDetJacobians(elJacs), refEl);
      EMap<EVector>(valsIPs.data(), valsIPs.size()) = phiMat * EMap<EVector>(valsNodes.data(), valsNodes.size());
      for(int ip = 0; ip < nIPs; ip++){
        integralVec += measure[ip]*EMap<EVector>(valsIPs.data() + ip*numValsPerObj, numValsPerObj);
      }
    }
  } else if(type == Edge){
    throw(ErrorHandle("Field", "integrate", "edge fields are not supported yet."));
  }
  if(part != NULL){
    MPI_Allreduce(MPI_IN_PLACE, integral->data(), numValsPerObj, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
};//integrate

} //hfox
