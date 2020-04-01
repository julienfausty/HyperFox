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
};//point/connectivity constructor

Mesh::Mesh(int dim, int order, std::string geom, moab::EntityHandle meshset_candidate): refElement(NULL){
  setReferenceElement(dim, order, geom);
  initializeMBInterface();
  setMeshSet(meshset_candidate);
};//meshset constructor

Mesh::~Mesh(){
  if(refElement != NULL){
    delete refElement;
  }
  moab::ErrorCode mbErr = mbInterface->clear_meshset(&meshset, 1);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "Destructor", "could not clear meshset."));
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
  computeFaces();
};//setPoints

void Mesh::setMeshSet(moab::EntityHandle meshset_candidate){
  moab::ErrorCode mbErr;
  mbErr = mbInterface->clear_meshset(&meshset, 1);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "setMeshSet", "could not clear current meshset (moab)."));
  }
  meshset = meshset_candidate;
  computeFaces();
};//setMeshSet

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
      {
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
      }
    case 3:
      {
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
  }
  return type;
};//determineMOABType

moab::EntityType Mesh::determineMOABTypeFace() const{
  moab::EntityType type;
  int dim = refElement->getDimension();
  elementGeometry elemGeom = refElement->getGeometry();
  switch(dim){
    case 1:
      {
      type = moab::MBVERTEX;
      break;
      }
    case 2:
      {
      type = moab::MBEDGE;
      break;
      }
    case 3:
      {
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
      }
  }
  return type;
};//determineMOABTypeFace

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

void Mesh::getCells(std::vector< std::vector<int> > * cells) const{
  moab::ErrorCode mbErr;
  moab::Range elems;
  int dim = refElement->getDimension();
  mbErr = mbInterface->get_entities_by_dimension(meshset, dim, elems);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getCells", "could not get elements from meshset"));
  }
  cells->resize(elems.size());
  for(moab::Range::iterator it = elems.begin(); it != elems.end(); it++){
    std::vector<moab::EntityHandle> verts;
    mbErr = mbInterface->get_connectivity(&(*it), 1, verts);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("Mesh", "getCells", "could not get vertices from element"));
    }
    int index = mbInterface->id_from_handle(*it) - 1;
    (*cells)[index].resize(verts.size());
    for(int i = 0; i < verts.size(); i++){
      (*cells)[index][i] = mbInterface->id_from_handle(verts[i])-1;
    }
  }
};//getCells


void Mesh::getSlicePoints(const std::vector<int> & slice, std::vector< std::vector<double> > * points) const{
  moab::ErrorCode mbErr;
  std::vector<moab::EntityHandle> vertexes(slice.size());
  for(int i = 0; i < slice.size(); i++){
    mbErr = mbInterface->handle_from_id(moab::MBVERTEX, slice[i]+1, vertexes[i]);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("Mesh", "getSlicePoints", "could not get vertexes from meshset"));
    }
  }
  points->resize(slice.size());
  int dim = refElement->getDimension();
  std::vector<double> vs(3*slice.size());
  mbErr = mbInterface->get_coords(vertexes.data(), vertexes.size(), vs.data());
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getSlicePoints", "could not get vertex coords from list"));
  }
  for(int i = 0; i < slice.size(); i++){
    (*points)[i] = std::vector<double>(vs.begin() + 3*i, vs.begin()+ 3*i + dim);
  }
};//getPoints

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


void Mesh::getSliceCells(const std::vector<int> & slice, std::vector< std::vector<int> > * cells) const{
  moab::ErrorCode mbErr;
  std::vector<moab::EntityHandle> elements(slice.size());
  for(int i = 0; i < slice.size(); i++){
    mbErr = mbInterface->handle_from_id(determineMOABType(), slice[i]+1, elements[i]);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("Mesh", "getSliceCell", "could not get element from meshset."));
    }
  }
  int nNodesEl = refElement->getNumNodes();
  std::vector<moab::EntityHandle> verts(elements.size()*nNodesEl);
  mbErr = mbInterface->get_connectivity(elements.data(), elements.size(), verts);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getSliceCell", "could not get vertices from element"));
  }
  cells->resize(elements.size());
  for(int i = 0; i < elements.size(); i++){
    for(int j = 0; j < nNodesEl; j++){
      (*cells)[i][j] = mbInterface->id_from_handle(verts[i*nNodesEl + j]) - 1;
    }
  }
};//getSliceCell

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

void Mesh::computeFaces(){
  moab::ErrorCode mbErr;
  moab::Range cells;
  mbErr = mbInterface->get_entities_by_dimension(meshset, refElement->getDimension(), cells);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "computeFaces", "could not get cells from meshset (moab)"));
  }
  moab::Range faces;
  mbErr = mbInterface->get_adjacencies(cells, refElement->getDimension() - 1, 1, faces, moab::Interface::UNION);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "computeFaces", "could not develop face adjacencies from cells (moab)"));
  }
  mbErr = mbInterface->add_entities(meshset, faces);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "computeFaces", "could not add faces to meshset (moab)"));
  }
  moab::Range adj;
  mbErr = mbInterface->get_adjacencies(faces, refElement->getDimension(), 1, adj, moab::Interface::UNION);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "computeFaces", "could not get adjacencies from faces to cells (moab)"));
  }
  computeFace2CellMap();
  computeCell2FaceMap();
};//computeFaces

int Mesh::getNumberFaces() const{
  moab::ErrorCode mbErr;
  moab::Range faces;
  int dim = refElement->getDimension();
  mbErr = mbInterface->get_entities_by_dimension(meshset, dim-1, faces);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getNumberFaces", "could not get faces from meshset"));
  }
  return faces.size();
};//getNumberFaces

void Mesh::getFaces(std::vector< std::vector<int> > * faces) const{
  moab::ErrorCode mbErr;
  moab::Range fcs;
  int dim = refElement->getDimension();
  mbErr = mbInterface->get_entities_by_dimension(meshset, dim-1, fcs);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getFaces", "could not get faces from meshset"));
  }
  faces->resize(fcs.size());
  for(moab::Range::iterator it = fcs.begin(); it != fcs.end(); it++){
    std::vector<moab::EntityHandle> verts;
    mbErr = mbInterface->get_connectivity(&(*it), 1, verts);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("Mesh", "getFaces", "could not get vertices from face"));
    }
    int index = mbInterface->id_from_handle(*it) - 1;
    (*faces)[index].resize(verts.size());
    for(int i = 0; i < verts.size(); i++){
      (*faces)[index][i] = mbInterface->id_from_handle(verts[i])-1;
    }
  }
};//getFaces

void Mesh::getSliceFaces(const std::vector<int> & slice, std::vector< std::vector<int> > * faces) const{
  moab::ErrorCode mbErr;
  std::vector<moab::EntityHandle> fcs(slice.size());
  for(int i = 0; i < slice.size(); i++){
    mbErr = mbInterface->handle_from_id(determineMOABTypeFace(), slice[i]+1, fcs[i]);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("Mesh", "getSliceFaces", "could not get face from meshset."));
    }
  }
  int nNodesFace = refElement->getFaceElement()->getNumNodes();
  std::vector<moab::EntityHandle> verts(fcs.size()*nNodesFace);
  mbErr = mbInterface->get_connectivity(fcs.data(), fcs.size(), verts);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getSliceFaces", "could not get vertices from faces"));
  }
  faces->resize(fcs.size());
  for(int i = 0; i < fcs.size(); i++){
    for(int j = 0; j < nNodesFace; j++){
      (*faces)[i][j] = mbInterface->id_from_handle(verts[i*nNodesFace + j]) - 1;
    }
  }
};//getSliceFaces


void Mesh::getFace(int i, std::vector<int> * face) const{
  moab::ErrorCode mbErr;
  moab::EntityHandle ent;
  mbErr = mbInterface->handle_from_id(determineMOABTypeFace(), i+1, ent);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getFace", "could not get face from meshset."));
  }
  std::vector<moab::EntityHandle> verts(refElement->getFaceElement()->getNumNodes());
  mbErr = mbInterface->get_connectivity(&ent, 1, verts);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getCell", "could not get vertices from element"));
  }
  face->resize(verts.size());
  for(int i = 0; i < verts.size(); i++){
    (*face)[i] = mbInterface->id_from_handle(verts[i]) - 1;
  }
};//getFace

void Mesh::computeFace2CellMap(){
  moab::ErrorCode mbErr;
  moab::Range faces;
  mbErr = mbInterface->get_entities_by_dimension(meshset, refElement->getDimension() - 1, faces);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getFace2CellMap", "could not get faces from meshset"));
  }
  face2CellMap.resize(faces.size());
  for(moab::Range::iterator it = faces.begin(); it != faces.end(); it++){
    moab::Range adj;
    mbErr = mbInterface->get_adjacencies(faces, refElement->getDimension(), 0, adj);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("Mesh", "getFace2CellMap", "could not develop cell adjacencies from face (moab)"));
    }
    int index = mbInterface->id_from_handle(*it) - 1;
    face2CellMap[index].resize(0);
    for(moab::Range::iterator itadj = adj.begin(); itadj != adj.end(); itadj++){
      face2CellMap[index].push_back(mbInterface->id_from_handle(*itadj) - 1);
    }
  }
};//computeFace2CellMap

void Mesh::computeCell2FaceMap(){
  moab::ErrorCode mbErr;
  moab::Range cells;
  mbErr = mbInterface->get_entities_by_dimension(meshset, refElement->getDimension(), cells);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getCell2FaceMap", "could not get cells from meshset"));
  }
  cell2FaceMap.resize(cells.size());
  std::cout << "Number of Faces: " << getNumberFaces() << std::endl;
  std::vector< std::vector<int> > faces;
  getFaces(&faces);
  for(int i = 0; i < faces.size(); i++){
    std::cout << "face " << i << std::endl;
    for(int j = 0; j < faces[i].size(); j++){
      std::cout << faces[i][j] << std::endl;
    }
  }
  for(moab::Range::iterator itcell = cells.begin(); itcell != cells.end(); itcell++){
    moab::Range adj;
    mbErr = mbInterface->get_adjacencies(&(*itcell), 1, refElement->getDimension() - 1, 0, adj);
    //std::cout << "adjacencies found for cell " << mbInterface->id_from_handle(*itcell)-1 << std::endl;
    //adj.print();
    //std::vector<double> coords(3);
    //for(moab::Range::iterator itadj = adj.begin(); itadj != adj.end(); itadj++){
      //std::vector<moab::EntityHandle> verts;
      //mbErr = mbInterface->get_connectivity(&(*itadj), 1, verts);
      //for(int k = 0; k < verts.size(); k++){
        //mbErr = mbInterface->get_coords(&(verts[k]), 1, coords.data());
        //std::cout << "face coords " << k << std::endl; 
        //for(int i = 0; i < refElement->getDimension(); i++){
          //std::cout << coords[i] << std::endl;
        //}
      //}
    //}
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("Mesh", "getCell2FaceMap", "could not get face adjacencies from cell (moab)"));
    }
    int index = mbInterface->id_from_handle(*itcell) - 1;
    cell2FaceMap[index].resize(adj.size());
    for(moab::Range::iterator itadj = adj.begin(); itadj != adj.end(); itadj++){
      int indexFace = getFaceIndex(*itcell, *itadj);
      cell2FaceMap[index][indexFace] = mbInterface->id_from_handle(*itadj) - 1;
    }
  }
};//computeCell2FaceMap

const std::vector< std::vector<int> > * Mesh::getFace2CellMap() const{
  return &face2CellMap;
};//getFace2CellMap

const std::vector< std::vector<int> > * Mesh::getCell2FaceMap() const{
  return &cell2FaceMap;
};//getFace2CellMap

const std::vector<int> * Mesh::getFace2Cell(int i) const{
  return &(face2CellMap[i]);
};//getCell2Face

const std::vector<int> * Mesh::getCell2Face(int i) const{
  return &(cell2FaceMap[i]);
};//getCell2Face

int Mesh::getFaceIndex(const moab::EntityHandle & cell, const moab::EntityHandle & face) const{
  moab::ErrorCode mbErr;
  std::vector<moab::EntityHandle> cellVerts(refElement->getNumNodes());
  mbErr = mbInterface->get_connectivity(&cell, 1, cellVerts);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getFaceIndex", "could not get vertices from cell (moab)"));
  }
  std::vector<moab::EntityHandle> faceVerts(refElement->getFaceElement()->getNumNodes());
  mbErr = mbInterface->get_connectivity(&face, 1, faceVerts);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getFaceIndex", "could not get vertices from face (moab)"));
  }
  std::vector<int> indexes(faceVerts.size());
  std::vector<moab::EntityHandle>::iterator it;
  for(int i = 0; i < faceVerts.size(); i++){
    it = std::find(cellVerts.begin(), cellVerts.end(), faceVerts[i]);
    if(it == cellVerts.end()){
      throw(ErrorHandle("Mesh", "getFaceIndex", "could not find one of the vertex faces in the cell."));
    }
    indexes[i] = std::distance(cellVerts.begin(), it);
  }
  const std::vector< std::vector<int> > * faceNodes = refElement->getFaceNodes();
  bool isFace = 0;
  int index = -1;
  for(int i = 0; i < faceNodes->size(); i++){
    for(int j = 0; j < indexes.size(); j++){
      isFace = (std::find((*faceNodes)[i].begin(), (*faceNodes)[i].end(), indexes[j]) != (*faceNodes)[i].end());
      if(!isFace){
        break;
      }
    }
    if(isFace){
      index = i;
      break;
    }
  }
  if(index == -1){
    throw(ErrorHandle("Mesh", "getFaceIndex", "could not find corresponding face"));
  }
  return index;
};//getFaceIndex

}//hfox
