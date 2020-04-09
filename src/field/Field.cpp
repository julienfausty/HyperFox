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
    values.resize(length);
  }
};//allocate

} //hfox
