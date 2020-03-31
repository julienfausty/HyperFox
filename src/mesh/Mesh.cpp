#include "Mesh.h"

namespace hfox{

Mesh::Mesh(): refElement(NULL){
  initializeMBInterface();
};//Empty constructor

Mesh::Mesh(int dim, int order, std::string geom, std::vector< std::vector<double> > & point_candidate, 
    std::vector< std::vector<int> > & connectivity_candidate): refElement(NULL){
  setReferenceElement(dim, order, geom);
  initializeMBInterface();
  setMesh(point_candidate, connectivity_candidate);
  computeFaces();
};//point/connectivity constructor

Mesh::~Mesh(){
  if(refElement != NULL){
    delete refElement;
  }
  delete mbInterface;
};//Destructor

void Mesh::setReferenceElement(int dim, int order, std::string geom){
  if(refElement != NULL){
    delete refElement;
  }
  refElement = new ReferenceElement(dim, order, geom);
};//setReferenceElement

void Mesh::initializeMBInterface(){
  moab::ErrorCode mbErr;
  mbInterface = new (std::nothrow) moab::Core;
  if(mbInterface == NULL){
    throw(ErrorHandle("Mesh", "Constructor", "could not construct new moab mesh."));
  }
  mbErr = mbInterface->create_meshset(moab::MESHSET_SET, meshset);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "Constructor", "could not initialize meshset."));
  }
};//initializeMBMesh

void Mesh::setMesh(std::vector< std::vector<double> > & points_candidate, 
    std::vector< std::vector<int> > & connectivity_candidate){
  moab::ErrorCode mbErr;
  if(refElement == NULL){
    throw("Mesh", "setMesh", "cannot set mesh before reference element is set.");
  }
  std::vector<double> point(3, 0.0);
  std::vector<moab::EntityHandle> cell;
  std::vector<moab::EntityHandle> orderedVertList;
  std::vector<moab::EntityHandle> orderedElList;
  for(int i = 0; i < points_candidate.size(); i++){
    moab::EntityHandle ent;
    for(int j = 0; j < points_candidate[i].size(); j++){
      point[j] = points_candidate[i][j];
    }
    mbErr = mbInterface->create_vertex(point.data(), ent);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("Mesh", "setMesh", "could not create vertex (moab)."));
    }
    orderedVertList.push_back(ent);
  }
  mbErr = mbInterface->add_entities(meshset, orderedVertList.data(), orderedVertList.size());
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "setMesh", "could not set vertexes in mesh (moab)."));
  }
  moab::EntityType type = determineMOABType();
  for(int i = 0; i < connectivity_candidate.size(); i++){
    moab::EntityHandle ent;
    cell.resize(connectivity_candidate[i].size());
    for(int j = 0; j < connectivity_candidate[i].size(); j++){
      cell[j] = orderedVertList[connectivity_candidate[i][j]];
    }
    mbErr = mbInterface->create_element(type, cell.data(), cell.size(), ent);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("Mesh", "setMesh", "could not create element (moab)."));
    }
    orderedElList.push_back(ent);
  }
  mbErr = mbInterface->add_entities(meshset, orderedElList.data(), orderedElList.size());
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "setMesh", "could not set elements in mesh (moab)."));
  }
};//setPoints

moab::EntityType Mesh::determineMOABType() const{
  moab::EntityType type;
  int dim = refElement->getDimension();
  elementGeometry elemGeom = refElement->getGeometry();
  switch(dim){
    case 1:
      {
      type = moab::MBEDGE;
      break;
      }
    case 2:
      switch(elemGeom){
        case simplex:
          {
            type = moab::MBTRI;
            break;
          }
        case orthotope:
          {
            type = moab::MBQUAD;
            break;
          }
      }
      break;
    case 3:
      switch(elemGeom){
        case simplex:
          {
            type = moab::MBTET;
            break;
          }
        case orthotope:
          {
            type = moab::MBHEX;
            break;
          }
      }
      break;
  }
  return type;
};//determineMOABType

void Mesh::getPoints(std::vector< std::vector<double> > * points) const{
  moab::ErrorCode mbErr;
  moab::Range vertexes;
  mbErr = mbInterface->get_entities_by_dimension(meshset, 0, vertexes);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getPoints", "could not get vertexes from meshset"));
  }
  points->resize(vertexes.size());
  int dim = refElement->getDimension();
  double pV[3];
  for(moab::Range::iterator it = vertexes.begin(); it != vertexes.end(); it++){
    mbErr = mbInterface->get_coords(&(*it), 1, pV);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("Mesh", "getPoints", "could not get vertex coords from list"));
    }
    int index = mbInterface->id_from_handle(*it) -1;
    (*points)[index] = std::vector<double>(pV, pV+dim);
  }
};//getPoints

void Mesh::getConnectivity(std::vector< std::vector<int> > * cells) const{
  moab::ErrorCode mbErr;
  moab::Range elems;
  int dim = refElement->getDimension();
  mbErr = mbInterface->get_entities_by_dimension(meshset, dim, elems);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getConnectivity", "could not get elements from meshset"));
  }
  cells->resize(elems.size());
  for(moab::Range::iterator it = elems.begin(); it != elems.end(); it++){
    std::vector<moab::EntityHandle> verts;
    mbErr = mbInterface->get_connectivity(&(*it), 1, verts);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("Mesh", "getConnectivity", "could not get vertices from element"));
    }
    int index = mbInterface->id_from_handle(*it) - 1;
    (*cells)[index].resize(verts.size());
    for(int i = 0; i < verts.size(); i++){
      (*cells)[index][i] = mbInterface->id_from_handle(verts[i])-1;
    }
  }
};//getConnectivity

void Mesh::getPoint(int i, std::vector<double> * point) const{
  moab::ErrorCode mbErr;
  moab::EntityHandle ent;
  mbErr = mbInterface->handle_from_id(moab::MBVERTEX, i+1, ent);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getPoint", "could not get vertex from meshset."));
  }
  double pV[3];
  mbErr = mbInterface->get_coords(&ent, 1, pV);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getPoint", "could not get coords from vertex."));
  }
  int dim = refElement->getDimension();
  *point = std::vector<double>(pV, pV+dim);
};//getPoint

void Mesh::getCell(int i, std::vector<int> * cell) const{
  moab::ErrorCode mbErr;
  moab::EntityHandle ent;
  mbErr = mbInterface->handle_from_id(determineMOABType(), i+1, ent);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getCell", "could not get element from meshset."));
  }
  std::vector<moab::EntityHandle> verts(refElement->getNumNodes());
  mbErr = mbInterface->get_connectivity(&ent, 1, verts);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getCell", "could not get vertices from element"));
  }
  cell->resize(verts.size());
  for(int i = 0; i < verts.size(); i++){
    (*cell)[i] = mbInterface->id_from_handle(verts[i]) - 1;
  }
};//getCell

const ReferenceElement * Mesh::getReferenceElement() const{
  return refElement;
};//getReferenceElement

int Mesh::getNumberPoints() const{
  moab::ErrorCode mbErr;
  moab::Range vertexes;
  mbErr = mbInterface->get_entities_by_dimension(meshset, 0, vertexes);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getNumberPoints", "could not get vertexes from meshset"));
  }
  return vertexes.size();
};//getNumberPoints

int Mesh::getNumberCells() const{
  moab::ErrorCode mbErr;
  moab::Range elems;
  int dim = refElement->getDimension();
  mbErr = mbInterface->get_entities_by_dimension(meshset, dim, elems);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getNumberCells", "could not get cells from meshset"));
  }
  return elems.size();
};//getNumberCells

int Mesh::getDimension() const{
  return refElement->getDimension();
};//getDimension

void Mesh::getInnerFaces(std::vector< std::array<int, 4> > * innerFaces) const{

};//getInnerFaces

void Mesh::getInnerFace(int i, std::array<int, 4> * innerFace) const{

};//getInnerFace

void Mesh::getOuterFaces(std::vector< std::array<int, 2> > * outerFaces) const{

};//getOuterFaces

void Mesh::getOuterFace(int i, std::array<int, 2> * outerFace) const{

};//getOuterFace

void Mesh::computeFaces(){
};//computeFaces

}//hfox
