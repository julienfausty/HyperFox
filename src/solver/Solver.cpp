#include "Solver.h"

namespace hfox{

void Solver::initialize(){
  if(linSystem != NULL){
    linSystem->destroySystem();
    linSystem->initialize();
    linSystem->configure();
  };
  initialized = 1;
};//initialize

void Solver::constructLocalFields(std::vector<int> & entities, std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::iterator it;
  for(it = fm->begin(); it != fm->end(); it++){
    (*fieldMap)[it->first]->getSliceValues(entities, &(it->second));
  }
};//constructLocalFields

void Solver::constructLocalFields(std::vector<int> & entities, std::vector<int> locEntities, std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::iterator it;
  std::vector<double> buffer;
  int buffSize;
  for(it = fm->begin(); it != fm->end(); it++){
    Field * F = fieldMap->at(it->first);
    buffSize = *(F->getNumObjPerEnt()) * (*(F->getNumValsPerObj()));
    buffer.resize(buffSize, 0);
    it->second.resize(buffSize*entities.size(),0);
    for(int i = 0; i < entities.size(); i++){
      if(locEntities[i] != -1){
        F->getValues(locEntities[i], &buffer);
      } else { 
        F->getParValues(entities[i], &buffer);
      }
      std::copy(buffer.begin(), buffer.end(), it->second.begin() + buffSize*i);
    }
  }
};//constructLocalFields

std::map<std::string, std::vector<double> > Solver::prepareLocalFieldMap(FieldType ft){
  const ReferenceElement * refEl = myMesh->getReferenceElement();
  std::map<std::string, std::vector<double> > localFieldMap;
  int elSize = 0;
  switch(ft){
    case Node:{elSize = refEl->getNumNodes() * nDOFsPerNode; break;}
    case Face:{elSize = refEl->getFaceElement()->getNumNodes() * nDOFsPerNode; break;}
    case Cell:{elSize = refEl->getNumNodes() * nDOFsPerNode; break;}
  }
  std::map<std::string, Field * >::iterator it;
  for(it = fieldMap->begin(); it != fieldMap->end(); it++){
    if(*(it->second->getFieldType()) == ft){
      localFieldMap[it->first] = std::vector<double>(elSize, 0.0);
    }
  }
  return localFieldMap;
};//prepareLocalFieldMap

void Solver::reorder(std::vector<double> * solutionVec, std::vector<int> * thisRange, std::vector<int> * dofsWeNeed){
  Partitioner * part = myMesh->getPartitioner();
  std::vector<int> currentMap, sendIds;
  std::vector<double> sendVals;
  std::vector<int>::iterator it;
  for(int i = 0; i < thisRange->size(); i++){
    it = std::find(dofsWeNeed->begin(), dofsWeNeed->end(), thisRange->at(i));
    if(it != dofsWeNeed->end()){
      currentMap.push_back(i);
      currentMap.push_back(std::distance(dofsWeNeed->begin(), it));
    } else {
      sendIds.push_back(thisRange->at(i));
      sendVals.push_back(solutionVec->at(i));
    }
  }
  //do the reordering we can currently do on this partition
  std::vector<double> buffer(solutionVec->size(), 0.0);
  for(int i = 0; i < currentMap.size()/2; i++){
    buffer[currentMap[2*i+1]] = solutionVec->at(currentMap[2*i]);
  }
  std::copy(buffer.begin(), buffer.end(), solutionVec->begin());
  buffer.clear();
  std::vector<int> recvIds;
  for(int i = 0; i < dofsWeNeed->size(); i++){
    it = std::find(thisRange->begin(), thisRange->end(), dofsWeNeed->at(i));
    if(it == thisRange->end()){
      recvIds.push_back(dofsWeNeed->at(i));
      recvIds.push_back(i);
    }
  }
  //gather the mismatched dofs to all nodes
  int nSendDofs = sendIds.size();
  std::vector<int> nDofs(part->getNumPartitions());
  std::vector<int> offsets(part->getNumPartitions());
  MPI_Allgather(&nSendDofs, 1, MPI_INT, nDofs.data(), 1, MPI_INT, MPI_COMM_WORLD);
  for(int i = 0; i < nDofs.size(); i++){
    offsets[i] = std::accumulate(nDofs.begin(), nDofs.begin() + i, 0);
  }
  std::vector<int> allDofIds(std::accumulate(nDofs.begin(), nDofs.end(), 0));
  std::vector<double> allDofVals(allDofIds.size());
  MPI_Allgatherv(sendIds.data(), sendIds.size(), MPI_INT, allDofIds.data(), nDofs.data(), offsets.data(), MPI_INT, MPI_COMM_WORLD);
  MPI_Allgatherv(sendVals.data(), sendVals.size(), MPI_DOUBLE, allDofVals.data(), nDofs.data(), offsets.data(), MPI_DOUBLE, MPI_COMM_WORLD);
  for(int i = 0; i < recvIds.size()/2; i++){
    it = std::find(allDofIds.begin(), allDofIds.end(), recvIds[2*i]);
    if(it != allDofIds.end()){
      (*solutionVec)[recvIds[2*i + 1]] = allDofVals[std::distance(allDofIds.begin(), it)];
    } else {
      throw(ErrorHandle("CGSolver", "solve", "one of the dof indexes of the partition could not be found in the communication structure."));
    }
  }

};

}//hfox
