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

} //hfox
