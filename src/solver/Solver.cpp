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

}//hfox
