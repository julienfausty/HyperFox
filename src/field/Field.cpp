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
  numObjPerEnt = nObjPerEnt;
  numValsPerObj = nValsPerObj;
  allocate();
};//allocation constructor

void Field::allocate(){
  if(pmesh != NULL){
    int length;
    switch(type){
      case None: {length = 0; break;}
      case Node: {length = pmesh->getNumberPoints(); break;}
      case Edge: {throw(ErrorHandle("Field", "allocate", "edge fields are not supported yet.")); break;}
      case Face: {length = pmesh->getNumberFaces(); break;}
      case Cell: {length = pmesh->getNumberCells(); break;}
    }
    length *= numObjPerEnt*numValsPerObj;
    values.resize(length);
  }
};//allocate

} //hfox
