#include "FEModel.h"

namespace hfox{

FEModel::FEModel(const ReferenceElement * re){
  refEl = re;
};//constructor

FEModel::~FEModel(){
  for(std::map<std::string, Operator*>::iterator it = operatorMap.begin(); it != operatorMap.end(); it++){
    delete it->second;
  }
};//destructor

void FEModel::setTimeScheme(TimeScheme * ts){
  if(allocated or fieldSet){
    throw(ErrorHandle("FEModel", "setTimeScheme", "the time scheme must be set before allocation or field setting"));
  }
  timeScheme = ts;
};//setTimeScheme

void FEModel::compute(){
  if(allocated){
    computeLocalMatrix();
    computeLocalRHS();
    if(timeScheme != NULL){
      timeScheme->apply(&localMatrix, &localRHS);
    }
  } else {
    throw(ErrorHandle("FEModel", "compute", "the FEModel must be allocated before computing."));
  }
};//compute

}//hfox
