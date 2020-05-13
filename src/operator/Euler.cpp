#include "Euler.h"

namespace hfox{

Euler::Euler(const ReferenceElement * re, bool isExplicitUser) : TimeScheme(re){
  isExplicit = isExplicitUser;
};//constructor

void Euler::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator it;
  it = fm->find("Solution");
  if(it == fm->end()){
    throw(ErrorHandle("Euler", "setFieldMap", "the field map must have a Solution field"));
  }
  fieldMap["Solution"] = &(it->second);
};//setFieldMap

void Euler::apply(EMatrix * stiffness, EVector * rhs){
  if(!allocated){
    throw(ErrorHandle("Euler", "apply", "the Euler time scheme should at least be allocated (and hopefully assembled) before being applied"));
  }
  if(fieldMap.empty()){
    throw(ErrorHandle("Euler", "apply", "the fieldMap must be set before attempting to apply the time scheme"));
  }
  if(deltat == 0.0){
    throw(ErrorHandle("Euler", "apply", "the time step needs to be set before applying and it should not be 0"));
  }
  //op *= (1/deltat);
  *stiffness *= deltat;
  *rhs *= deltat;
  //if(!isExplicit){
    //*stiffness += op;
    //*rhs += op*EMap<const EVector>(fieldMap["Solution"]->data(), op.cols());
  //} else {
    //*rhs += (op - (*stiffness))*EMap<const EVector>(fieldMap["Solution"]->data(), op.cols());
    //*stiffness = op;
  //}
  if(!isExplicit){
    stiffness->block(0, 0, op.rows(), op.cols()) += op;
    rhs->segment(0, op.rows()) += op*EMap<const EVector>(fieldMap["Solution"]->data(), op.cols());
  } else {
    rhs->segment(0, op.rows()) += (op - (*stiffness))*EMap<const EVector>(fieldMap["Solution"]->data(), op.cols());
    stiffness->block(0, 0, op.rows(), op.cols()) = op;
  }
};//apply

}//hfox
