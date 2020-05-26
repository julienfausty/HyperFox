#include "RungeKutta.h"

namespace hfox{

RungeKutta::RungeKutta(const ReferenceElement * re, RKType type, std::vector<std::string> fields) : TimeScheme(re){
  setButcherTable(type);
  setAuxiliaryFields(fields);
  stageCounter = 0;
};//constructor


void RungeKutta::setButcherTable(RKType type){
  if(butcherDB.empty()){
    setUpDB();
  };
  bTable = butcherDB[type];
};//setButcherTable type

void RungeKutta::setButcherTable(EMatrix butchTable){
  if(butchTable.cols() != butchTable.rows()){
    throw(ErrorHandle("RungeKutta", "setButcherTable", "the Butcher table should be square"));
  }
  for(int i = 0; i < butchTable.rows(); i++){
    for(int j = i+1; j < butchTable.cols(); j++){
      if(butchTable(i, j) != 0){
        throw(ErrorHandle("RungeKutta", "setButcherTable", "the upper triangular part of the Butcher" 
             " table should be null (no fully implicit implementation as of yet)"));
      }
    }
  }
  bTable = butchTable;
};//setButcherTable

void RungeKutta::setAuxiliaryFields(std::vector<std::string> fields){
  if(!fieldMap.empty()){
    throw(ErrorHandle("RungeKutta", "setAuxiliaryField", "the auxiliary fields must be specified before the field map"));
  }
  auxiliaryFields = fields;
};//setButcherTable type

int RungeKutta::getNumStages() const{
  return (bTable.cols() - 1);
};//getNumStages

void RungeKutta::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator it;
  it = fm->find("Solution");
  if(it == fm->end()){
    throw(ErrorHandle("RungeKutta", "setFieldMap", "the field map must have a Solution field"));
  }
  fieldMap["Solution"] = &(it->second);
  for(int i = 0; i < auxiliaryFields.size(); i++){
    it = fm->find(auxiliaryFields[i]);
    if(it == fm->end()){
      throw(ErrorHandle("RungeKutta", "setFieldMap", "the field map must provide the auxiliary fields. Failure to find: " + auxiliaryFields[i]));
    }
    fieldMap[auxiliaryFields[i]] = &(it->second);
  }
  it = fm->find("OldSolution");
  if(it == fm->end()){
    throw(ErrorHandle("RungeKutta", "setFieldMap", "the field map must have a OldSolution field"));
  }
  fieldMap["OldSolution"] = &(it->second);
  for(int i = 0; i < auxiliaryFields.size(); i++){
    it = fm->find("Old" + auxiliaryFields[i]);
    if(it == fm->end()){
      throw(ErrorHandle("RungeKutta", "setFieldMap", "the field map must provide the 'Old' auxiliary fields. Failure to find: Old" + auxiliaryFields[i]));
    }
    fieldMap["Old" + auxiliaryFields[i]] = &(it->second);
  }
  int nStages = getNumStages();
  for(int k = 0; k < nStages; k++){
    std::string nameField = "RKStage_" + std::to_string(k);
    it = fm->find(nameField);
    if(it == fm->end()){
      throw(ErrorHandle("RungeKutta", "setFieldMap", "the field map must have the same number of RKStage fields as stages"));
    }
    fieldMap[nameField] = &(it->second);
    for(int i = 0; i < auxiliaryFields.size(); i++){
      nameField = "RKStage_" + auxiliaryFields[i] + "_" + std::to_string(k);
      it = fm->find(nameField);
      if(it == fm->end()){
        throw(ErrorHandle("RungeKutta", "setFieldMap", "the field map must provide all the 'RKStage' auxiliary fields. Failure to find: " + nameField));
      }
      fieldMap[nameField] = &(it->second);
    }
  }
};//setFieldMap

void RungeKutta::apply(EMatrix * stiffness, EVector * s){
  if(!allocated){
    throw(ErrorHandle("RungeKutta", "apply", "the time scheme should at least be allocated (and hopefully assembled) before being applied"));
  }
  if(fieldMap.empty()){
    throw(ErrorHandle("RungeKutta", "apply", "the fieldMap must be set before attempting to apply the time scheme"));
  }
  if(deltat == 0.0){
    throw(ErrorHandle("RungeKutta", "apply", "the time step needs to be set before applying and it should not be 0"));
  }
  int lenRows = fieldMap["Solution"]->size();
  int lenCols = lenRows;
  for(int i = 0; i < auxiliaryFields.size(); i++){
    lenCols += fieldMap[auxiliaryFields[i]]->size();
  }
  if((stiffness->rows() != lenRows) or (stiffness->cols() != lenCols)){
    throw(ErrorHandle("RungeKutta", "apply", "the stiffness matrix does not have the correct dimensions, it should be (" + std::to_string(lenRows) + ", " + std::to_string(lenCols) + ")"));
  }
  if(s->size() != lenRows){
    throw(ErrorHandle("RungeKutta", "apply", "the right hand side matrix does not have the correct length, it should be " + std::to_string(lenRows)));
  }
  EMatrix thisStage = bTable.block(stageCounter, 1, 1, getNumStages());
  EVector uj = EVector::Zero(lenCols);
  std::string nameField;
  int start = 0;
  for(int i = 0; i < stageCounter; i++){
    nameField = "RKStage_" + std::to_string(i);
    uj.segment(start, lenRows) += thisStage(0, i)*EMap<const EVector>(fieldMap[nameField]->data(), lenRows);
    start += lenRows;
    for(int j = 0; j < auxiliaryFields.size(); j++){
      nameField = "RKStage_" + auxiliaryFields[j] + "_" + std::to_string(i);
      uj.segment(start, fieldMap[nameField]->size()) += thisStage(0, i)*EMap<const EVector>(fieldMap[nameField]->data(), fieldMap[nameField]->size());
      start += fieldMap[nameField]->size();
    }
  }
  uj *= deltat;
  EVector ut = EVector::Zero(lenCols);
  start = 0;
  ut.segment(start, lenRows) = EMap<const EVector>(fieldMap["OldSolution"]->data(), fieldMap["OldSolution"]->size());
  start += lenRows;
  for(int j = 0; j < auxiliaryFields.size(); j++){
    nameField = "Old" + auxiliaryFields[j];
    ut.segment(start, fieldMap[nameField]->size()) = EMap<const EVector>(fieldMap[nameField]->data(), fieldMap[nameField]->size());
    start += fieldMap[nameField]->size();
  }
  uj += ut;
  (*stiffness) *= deltat;
  (*s) *= deltat;
  //s->segment(0, op.rows()) -= stiffness->block(0, 0, op.rows(), op.cols())*uj;//explicit part
  (*s) -= (*stiffness)*uj;//explicit part
  //stiffness->block(0, 0, op.rows(), op.cols()) *= thisStage(0, stageCounter);//implicit part
  (*stiffness) *= thisStage(0, stageCounter);//implicit part
  stiffness->block(0, 0, lenRows, lenRows) += op;
  (*s) += (*stiffness)*ut;
};//apply

void RungeKutta::computeStage(std::map<std::string, Field*> * fm){
  if(stageCounter >= getNumStages()){
    throw(ErrorHandle("RungeKutta", "computeStage", "cannot compute more stages than the method allows, think about computing the solution"));
  }
  EMatrix thisStage = bTable.block(stageCounter, 1, 1, getNumStages());
  double invdt = 1.0/(deltat);
  double buffer;
  Field * rkF = fm->at("RKStage_" + std::to_string(stageCounter));
  Field * solF = fm->at("Solution");
  Field * oldSolF = fm->at("OldSolution");
  for(int i = 0; i < solF->getLength(); i++){
    rkF->getValues()->at(i) = invdt*(solF->getValues()->at(i) - oldSolF->getValues()->at(i));
    buffer = 0.0;
    for(int j = 0; j < stageCounter+1; j++){
      buffer += thisStage(0, j) * (fm->at("RKStage_" + std::to_string(j))->getValues()->at(i));
    }
    solF->getValues()->at(i) = oldSolF->getValues()->at(i) + deltat * buffer;
  }
  for(int k = 0; k < auxiliaryFields.size(); k++){
    rkF = fm->at("RKStage_"+ auxiliaryFields[k] + "_" + std::to_string(stageCounter));
    solF = fm->at(auxiliaryFields[k]);
    oldSolF = fm->at("Old"+auxiliaryFields[k]);
    for(int i = 0; i < solF->getLength(); i++){
      rkF->getValues()->at(i) = invdt*(solF->getValues()->at(i) - oldSolF->getValues()->at(i));
      buffer = 0.0;
      for(int j = 0; j < stageCounter+1; j++){
        buffer += thisStage(0, j) * (fm->at("RKStage_" + auxiliaryFields[k] + "_" + std::to_string(j))->getValues()->at(i));
      }
      solF->getValues()->at(i) = oldSolF->getValues()->at(i) + deltat * buffer;
    }
  }
  stageCounter += 1;
};//computeStage

void RungeKutta::computeSolution(std::map<std::string, Field*> * fm){
  if(stageCounter != getNumStages()){
    throw(ErrorHandle("RungeKutta", "computeSolution", "all stages must be computed before computing the solution"));
  }
  EMatrix bs = bTable.block(stageCounter, 1, 1, getNumStages());
  Field * solF = fm->at("Solution");
  Field * oldSolF = fm->at("OldSolution");
  double bkdt = 0.0;
  Field * rkF;
  for(int i = 0; i < solF->getLength(); i++){
    solF->getValues()->at(i) = oldSolF->getValues()->at(i);
  }
  for(int k = 0; k < getNumStages(); k++){
    bkdt = bs(0, k)*deltat;
    rkF = fm->at("RKStage_" + std::to_string(k));
    for(int i = 0; i < rkF->getLength(); i++){
      solF->getValues()->at(i) += bkdt*(rkF->getValues()->at(i));
    }
  }
  for(int j = 0; j < auxiliaryFields.size(); j++){
    solF = fm->at(auxiliaryFields[j]);
    oldSolF = fm->at("Old"+auxiliaryFields[j]);
    for(int i = 0; i < solF->getLength(); i++){
      solF->getValues()->at(i) = oldSolF->getValues()->at(i);
    }
    for(int k = 0; k < getNumStages(); k++){
      bkdt = bs(0, k)*deltat;
      rkF = fm->at("RKStage_" + auxiliaryFields[j] + "_" + std::to_string(k));
      for(int i = 0; i < rkF->getLength(); i++){
        solF->getValues()->at(i) += bkdt*(rkF->getValues()->at(i));
      }
    }
  }
  stageCounter = 0;
};//computeSolution

void RungeKutta::setUpDB(){
  butcherDB[FEuler].resize(2,2);
  butcherDB[FEuler] << 0, 0, 0, 1;
  butcherDB[EMidpoint].resize(3,3);
  butcherDB[EMidpoint] << 
    0, 0, 0, 
    0.5, 0.5, 0.0,
    0, 0, 1.0;
  butcherDB[Heun].resize(3,3);
  butcherDB[Heun] << 
    0, 0, 0, 
    1.0, 1.0, 0.0,
    0, 0.5, 0.5;
  butcherDB[Kutta3].resize(4,4);
  butcherDB[Kutta3] << 
    0, 0, 0, 0, 
    0.5, 0.5, 0.0, 0.0,
    1.0, -1.0, 2.0, 0.0,
    0.0, 1.0/6.0, 2.0/3.0, 1.0/6.0;
  butcherDB[Heun3].resize(4,4);
  butcherDB[Heun3] << 
    0, 0, 0, 0, 
    1.0/3.0, 1.0/3.0, 0.0, 0.0,
    2.0/3.0, 0.0, 2.0/3.0, 0.0,
    0.0, 1.0/4.0, 0.0, 3.0/4.0;
  butcherDB[SSPRK3].resize(4,4);
  butcherDB[SSPRK3] << 
    0, 0, 0, 0, 
    1.0, 1.0, 0.0, 0.0,
    1.0/2.0, 1.0/4.0, 1.0/4.0, 0.0,
    0.0, 1.0/6.0, 1.0/6.0, 2.0/3.0;
  butcherDB[RK4].resize(5,5);
  butcherDB[RK4] << 
    0, 0, 0, 0, 0, 
    1.0/2.0, 1.0/2.0, 0.0, 0.0, 0.0,
    1.0/2.0, 0.0, 1.0/2.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0;
  butcherDB[BEuler].resize(2,2);
  butcherDB[BEuler] << 1, 1, 0, 1;
  butcherDB[IMidpoint].resize(2,2);
  butcherDB[IMidpoint] << 1.0/2.0, 1.0/2.0, 0, 1;
  butcherDB[CrankNicolson].resize(3,3);
  butcherDB[CrankNicolson] << 
    0.0, 0.0, 0.0,
    1.0, 1.0/2.0, 1.0/2.0,
    0.0, 1.0/2.0, 1.0/2.0;
  butcherDB[KS2].resize(3,3);
  butcherDB[KS2] << 
    1.0/2.0, 1.0/2.0, 0.0,
    3.0/2.0, -1.0/2.0, 2.0,
    0.0, -1.0/2.0, 3.0/2.0;
  butcherDB[QZ2].resize(3,3);
  butcherDB[QZ2] << 
    1.0/4.0, 1.0/4.0, 0.0,
    3.0/4.0, 1.0/2.0, 1.0/4.0,
    0.0, 1.0/2.0, 1.0/2.0;
  butcherDB[ALX2].resize(3,3);
  double gamma = 1.0 - std::sqrt(2.0)/2.0;
  butcherDB[ALX2] << 
    gamma, gamma, 0.0,
    1.0, 1.0-gamma, gamma,
    0.0, 1.0-gamma, gamma;
  butcherDB[RK43].resize(5,5);
  butcherDB[RK43] << 
    1.0/2.0, 1.0/2.0, 0.0, 0.0, 0.0,
    2.0/3.0, 1.0/6.0, 1.0/2.0, 0.0, 0.0,
    1.0/2.0, -1.0/2.0, 1.0/2.0, 1.0/2.0, 0.0,
    1.0, 3.0/2.0, -3.0/2.0, 1.0/2.0, 1.0/2.0,
    0.0, 3.0/2.0, -3.0/2.0, 1.0/2.0, 1.0/2.0;
};//setUpDB

std::map<RKType, EMatrix> RungeKutta::butcherDB;

}//hfox
