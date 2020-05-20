#include "RungeKutta.h"

namespace hfox{

RungeKutta::RungeKutta(const ReferenceElement * re, RKType type) : TimeScheme(re){
  setButcherTable(type);
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
  it = fm->find("OldSolution");
  if(it == fm->end()){
    throw(ErrorHandle("RungeKutta", "setFieldMap", "the field map must have a OldSolution field"));
  }
  fieldMap["OldSolution"] = &(it->second);
  int nStages = getNumStages();
  for(int k = 0; k < nStages; k++){
    std::string nameField = "RKStage_" + std::to_string(k);
    it = fm->find(nameField);
    if(it == fm->end()){
      throw(ErrorHandle("RungeKutta", "setFieldMap", "the field map must have the same number of RKStage fields as stages"));
    }
    fieldMap[nameField] = &(it->second);
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
  EMatrix thisStage = bTable.block(stageCounter, 1, 1, getNumStages());
  EVector uj = EVector::Zero(fieldMap["OldSolution"]->size());
  std::string nameField;
  for(int i = 0; i < stageCounter; i++){
    nameField = "RKStage_" + std::to_string(i);
    uj += thisStage(0, i)*EMap<const EVector>(fieldMap[nameField]->data(), fieldMap[nameField]->size());
  }
  uj *= deltat;
  EMap<const EVector> ut(fieldMap["OldSolution"]->data(), fieldMap["OldSolution"]->size());
  uj += ut;
  (*stiffness) *= deltat;
  s->segment(0, op.rows()) *= deltat;
  s->segment(0, op.rows()) += op*ut - stiffness->block(0, 0, op.rows(), op.cols())*uj;//explicit part
  stiffness->block(0, 0, op.rows(), op.cols()) *= thisStage(0, stageCounter);//implicit part
  stiffness->block(0, 0, op.rows(), op.cols()) += op;
};//apply

void RungeKutta::computeStage(std::map<std::string, Field*> * fm){
  if(stageCounter >= getNumStages()){
    throw(ErrorHandle("RungeKutta", "computeStage", "cannot compute more stages than the method allows, think about computing the solution"));
  }
  double invdt = 1.0/(deltat);
  Field * rkF = fm->at("RKStage_" + std::to_string(stageCounter));
  Field * solF = fm->at("Solution");
  Field * oldSolF = fm->at("OldSolution");
  for(int i = 0; i < solF->getLength(); i++){
    rkF->getValues()->at(i) = invdt*(solF->getValues()->at(i) - oldSolF->getValues()->at(i));
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
  for(int i = 0; i < solF->getLength(); i++){
    solF->getValues()->at(i) = oldSolF->getValues()->at(i);
  }
  double bkdt = 0.0;
  Field * rkF;
  for(int k = 0; k < getNumStages(); k++){
    bkdt = bs(0, k)*deltat;
    rkF = fm->at("RKStage_" + std::to_string(k));
    for(int i = 0; i < rkF->getLength(); i++){
      solF->getValues()->at(i) += bkdt*(rkF->getValues()->at(i));
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
  butcherDB[BEuler].resize(2,2);
  butcherDB[BEuler] << 1.0/2.0, 1.0/2.0, 0, 1;
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
  butcherDB[QZ2].resize(5,5);
  butcherDB[QZ2] << 
    1.0/2.0, 1.0/2.0, 0.0, 0.0, 0.0,
    2.0/3.0, 1.0/6.0, 1.0/2.0, 0.0, 0.0,
    1.0/2.0, -1.0/2.0, 1.0/2.0, 1.0/2.0, 0.0,
    1.0, 3.0/2.0, -3.0/2.0, 1.0/2.0, 1.0/2.0,
    0.0, 3.0/2.0, -3.0/2.0, 1.0/2.0, 1.0/2.0;
};//setUpDB

std::map<RKType, EMatrix> RungeKutta::butcherDB;

}//hfox
