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
};

}//hfox
