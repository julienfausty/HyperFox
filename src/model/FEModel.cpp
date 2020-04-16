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

void FEModel::compute(){
  if(allocated){
    computeLocalMatrix();
    computeLocalRHS();
  } else {
    throw(ErrorHandle("FEModel", "compute", "the FEModel must be allocated before computing."));
  }
};//compute

}//hfox
