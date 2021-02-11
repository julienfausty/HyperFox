#include "Mesh.h"

namespace hfox{

Mesh::Mesh(): refElement(NULL), part(NULL){
};//Empty constructor

Mesh::Mesh(int dim, int order, std::string geom): refElement(NULL), part(NULL){
  setReferenceElement(dim, order, geom);
};//refElement constructor

Mesh::Mesh(int dim, int order, std::string geom, int dimPointSpace, std::vector<double> & point_candidate, 
    std::vector<int> & connectivity_candidate): refElement(NULL){
  setReferenceElement(dim, order, geom);
  setMesh(dimPointSpace, point_candidate, connectivity_candidate);
};//point/connectivity constructor

Mesh::~Mesh(){
  if(refElement != NULL){
    delete refElement;
  }
};//Destructor

void Mesh::setReferenceElement(int dim, int order, std::string geom){
  if(refElement != NULL){
    delete refElement;
  }
  refElement = new ReferenceElement(dim, order, geom);
};//setReferenceElement

void Mesh::setMesh(int dimPointSpace, std::vector<double> & points_candidate, 
    std::vector<int> & connectivity_candidate){
  moab::ErrorCode mbErr;
  if(refElement == NULL){
    throw("Mesh", "setMesh", "cannot set mesh before reference element is set.");
  }
  dimNodeSpace = dimPointSpace;
  nNodesPerCell = refElement->getNumNodes();
  nNodes = points_candidate.size()/dimNodeSpace;
  nCells = connectivity_candidate.size()/nNodesPerCell;
  nodes = points_candidate;
  nodes.resize(nodes.size());
  cells = connectivity_candidate;
  cells.resize(cells.size());
  computeFaces();
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

const std::vector<double> * Mesh::getPoints() const{
  return &nodes;
};//getPoints

const std::vector<int> * Mesh::getCells() const{
  return &cells;
};//getCells


void Mesh::getSlicePoints(const std::vector<int> & slice, std::vector< std::vector<double> > * points) const{
  points->resize(slice.size());
  for(int i = 0; i < slice.size(); i++){
    getPoint(slice[i], &((*points)[i]));
  }
};//getSlicePoints

void Mesh::getPoint(int i, std::vector<double> * point) const{
  point->resize(dimNodeSpace);
  point->assign(nodes.begin() + i*dimNodeSpace, nodes.begin() + (i+1)*dimNodeSpace);
};//getPoint


void Mesh::getSliceCells(const std::vector<int> & slice, std::vector< std::vector<int> > * userCells) const{
  userCells->resize(slice.size());
  for(int i = 0; i < slice.size(); i++){
    getCell(slice[i], &((*userCells)[i]));
  }
};//getSliceCell

void Mesh::getCell(int i, std::vector<int> * cell) const{
  cell->resize(nNodesPerCell);
  cell->assign(cells.begin() + i*nNodesPerCell, cells.begin() + (i+1)*nNodesPerCell);
};//getCell

const ReferenceElement * Mesh::getReferenceElement() const{
  return refElement;
};//getReferenceElement

int Mesh::getNumberPoints() const{
  return nNodes;
};//getNumberPoints

int Mesh::getNumberCells() const{
  return nCells;
};//getNumberCells

int Mesh::getDimension() const{
  return refElement->getDimension();
};//getDimension

int Mesh::getNodeSpaceDimension() const{
  return dimNodeSpace;
};//getDimension

void Mesh::computeFaces(){
  //creation of moab instance
  moab::ErrorCode mbErr;
  moab::Interface * mbInterface = new (std::nothrow) moab::Core;
  if(mbInterface == NULL){
    throw(ErrorHandle("Mesh", "computeFaces", "could not construct new moab mesh."));
  }
  moab::EntityHandle meshset;
  mbErr = mbInterface->create_meshset(moab::MESHSET_SET, meshset);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "computeFaces", "could not initialize meshset."));
  }
  //skeleton mesh creation
  std::vector<double> point(3, 0.0);
  std::vector<moab::EntityHandle> cell;
  std::vector<moab::EntityHandle> orderedVertList;
  std::vector<moab::EntityHandle> orderedElList;
  for(int i = 0; i < nNodes; i++){
    moab::EntityHandle ent;
    for(int j = 0; j < dimNodeSpace; j++){
      point[j] = nodes[i*dimNodeSpace + j];
    }
    mbErr = mbInterface->create_vertex(point.data(), ent);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("Mesh", "computeFaces", "could not create vertex (moab)."));
    }
    orderedVertList.push_back(ent);
  }
  mbErr = mbInterface->add_entities(meshset, orderedVertList.data(), orderedVertList.size());
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "computeFaces", "could not set vertexes in mesh (moab)."));
  }
  std::string strGeom;
  switch(refElement->getGeometry()){
    case(simplex): {strGeom = "simplex"; break;}
    case(orthotope): {strGeom = "orthotope"; break;}
  }
  int skelOrder = (refElement->getOrder() != 0) ? 1 : 0;
  ReferenceElement skelRefEl(refElement->getDimension(), skelOrder, strGeom);
  int nSkelNodes = skelRefEl.getNumNodes();
  moab::EntityType type = determineMOABType();
  for(int i = 0; i < nCells; i++){
    moab::EntityHandle ent;
    cell.resize(nSkelNodes);
    for(int j = 0; j < nSkelNodes; j++){
      cell[j] = orderedVertList[cells[i*nNodesPerCell + j]];
    }
    mbErr = mbInterface->create_element(type, cell.data(), cell.size(), ent);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("Mesh", "computeFaces", "could not create element (moab)."));
    }
    orderedElList.push_back(ent);
  }
  mbErr = mbInterface->add_entities(meshset, orderedElList.data(), orderedElList.size());
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "computeFaces", "could not set elements in mesh (moab)."));
  }
  //face computing
  moab::Range moabCells;
  mbErr = mbInterface->get_entities_by_dimension(meshset, refElement->getDimension(), moabCells);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "computeFaces", "could not get cells from meshset (moab)"));
  }
  moab::Range moabFaces;
  mbErr = mbInterface->get_adjacencies(moabCells, refElement->getDimension() - 1, 1, moabFaces, moab::Interface::UNION);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "computeFaces", "could not develop face adjacencies from cells (moab)"));
  }
  mbErr = mbInterface->add_entities(meshset, moabFaces);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "computeFaces", "could not add faces to meshset (moab)"));
  }
  computeCell2FaceMap(mbInterface, meshset, skelRefEl);
  computeFace2CellMap(mbInterface, meshset);
  computeBoundary(mbInterface, meshset);
  computeFaceConnectivity(mbInterface, meshset);
  //destruction of moab instance
  moab::Range allEnts;
  mbErr = mbInterface->get_entities_by_handle(meshset, allEnts);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "computeFaces", "could not get all entities"));
  }
  mbErr = mbInterface->delete_entities(allEnts);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "computeFaces", "could not delete all entities"));
  }
  mbErr = mbInterface->delete_entities(&meshset, 1);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "computeFaces", "could not clear meshset."));
  }
  delete mbInterface;
};//computeFaces

int Mesh::getNumberFaces() const{
  return nFaces;
};//getNumberFaces

const std::vector<int> * Mesh::getFaces() const{
  return &faces;
};//getFaces

void Mesh::getSliceFaces(const std::vector<int> & slice, std::vector< std::vector<int> > * userFaces) const{
  userFaces->resize(slice.size());
  for(int i = 0; i < slice.size(); i++){
    getFace(slice[i], &((*userFaces)[i]));
  }
};//getSliceFaces


void Mesh::getFace(int i, std::vector<int> * face) const{
  face->resize(nNodesPerFace);
  face->assign(faces.begin() + i*nNodesPerFace, faces.begin() + (i+1)*nNodesPerFace);
};//getFace

void Mesh::getGhostPoint(int i, std::vector<double> * point) const{
  point->resize(dimNodeSpace, 0.0);
  bool found = 0;
  int nDataPNode = 2 + dimNodeSpace;
  for(int k = 0; k < ghostNodes.size()/nDataPNode; k++){
    if(ghostNodes[k*nDataPNode + 1] == i){
      point->assign(ghostNodes.begin() + k*nDataPNode + 2, ghostNodes.begin() + k*nDataPNode + 2 + dimNodeSpace);
      found = 1;
      break;
    }
  }
  if(!found){
    throw(ErrorHandle("Mesh", "getGhostPoint", "node index " + std::to_string(i) + " not found in ghostNodes"));
  }
};//getGhostPoint

void Mesh::getGhostCell(int i, std::vector<int> * cell) const{
  cell->resize(nNodesPerCell, 0);
  bool found = 0;
  int nDataPCell = 2 + nNodesPerCell;
  for(int k = 0; k < ghostCells.size()/nDataPCell; k++){
    if(ghostCells[k*nDataPCell + 1] == i){
      cell->assign(ghostCells.begin() + k*nDataPCell + 2, ghostCells.begin() + k*nDataPCell + 2 + nNodesPerCell);
      found = 1;
      break;
    }
  }
  if(!found){
    throw(ErrorHandle("Mesh", "getGhostCell", "cell index " + std::to_string(i) + " not found in ghostCells"));
  }
};//getGhostCell

void Mesh::getGhostFace(int i, std::vector<int> * cell) const{
  cell->resize(nNodesPerFace, 0);
  bool found = 0;
  int nDataPCell = 2 + nNodesPerFace;
  for(int k = 0; k < ghostFaces.size()/nDataPCell; k++){
    if(ghostFaces[k*nDataPCell + 1] == i){
      cell->assign(ghostFaces.begin() + k*nDataPCell + 2, ghostFaces.begin() + k*nDataPCell + 2 + nNodesPerFace);
      found = 1;
      break;
    }
  }
  if(!found){
    throw(ErrorHandle("Mesh", "getGhostFace", "face index " + std::to_string(i) + " not found in ghostFaces"));
  }
};//getGhostCell

void Mesh::getGhostFace2Cell(int i, std::vector<int> * face2Cell) const{
  face2Cell->resize(2, 0);
  bool found = 0;
  int nDataP = 4;
  for(int k = 0; k < ghostFace2CellMap.size()/nDataP; k++){
    if(ghostFace2CellMap[k*nDataP + 1] == i){
      face2Cell->assign(ghostFace2CellMap.begin() + k*nDataP + 2, ghostFace2CellMap.begin() + k*nDataP + 4);
      found = 1;
      break;
    }
  }
  if(!found){
    throw(ErrorHandle("Mesh", "getGhostFace2Cell", "face index " + std::to_string(i) + " not found in ghostFace2CellMap"));
  }
};//getGhostFace2Cell

void Mesh::getGhostCell2Face(int i, std::vector<int> * cell2Face) const{
  cell2Face->resize(nFacesPerCell, 0);
  bool found = 0;
  int nDataP = 2 + nFacesPerCell;
  for(int k = 0; k < ghostCell2FaceMap.size()/nDataP; k++){
    if(ghostCell2FaceMap[k*nDataP + 1] == i){
      cell2Face->assign(ghostCell2FaceMap.begin() + k*nDataP + 2, ghostCell2FaceMap.begin() + (k+1)*nDataP);
      found = 1;
      break;
    }
  }
  if(!found){
    throw(ErrorHandle("Mesh", "getGhostCell2Face", "cell index " + std::to_string(i) + " not found in ghostCell2FaceMap"));
  }
};//getGhostFace2Cell

void Mesh::computeFace2CellMap(moab::Interface * mbInterface, moab::EntityHandle & meshset){
  moab::ErrorCode mbErr;
  moab::Range moabFaces;
  mbErr = mbInterface->get_entities_by_dimension(meshset, refElement->getDimension() - 1, moabFaces);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "computeFace2CellMap", "could not get faces from meshset"));
  }
  face2CellMap.resize(moabFaces.size() * 2);
  for(moab::Range::iterator it = moabFaces.begin(); it != moabFaces.end(); it++){
    moab::Range adj;
    mbErr = mbInterface->get_adjacencies(&(*it), 1, refElement->getDimension(), 0, adj);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("Mesh", "computeFace2CellMap", "could not develop cell adjacencies from face (moab)"));
    }
    int index = mbInterface->id_from_handle(*it) - 1;
    int locInd = 0;
    for(moab::Range::iterator itadj = adj.begin(); itadj != adj.end(); itadj++){
      face2CellMap[index *2 + locInd] = mbInterface->id_from_handle(*itadj) - 1;
      if(adj.size() == 2){
        locInd += 1;
      } else if(adj.size() == 1){
        face2CellMap[index * 2 + 1] = -1;
      }
    }
  }
};//computeFace2CellMap

void Mesh::computeCell2FaceMap(moab::Interface * mbInterface, moab::EntityHandle & meshset, ReferenceElement & skelRefEl){
  moab::ErrorCode mbErr;
  moab::Range moabCells;
  mbErr = mbInterface->get_entities_by_dimension(meshset, refElement->getDimension(), moabCells);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "computeCell2FaceMap", "could not get cells from meshset"));
  }
  nFacesPerCell = refElement->getFaceNodes()->size();
  cell2FaceMap.resize(moabCells.size() * nFacesPerCell);
  for(moab::Range::iterator itcell = moabCells.begin(); itcell != moabCells.end(); itcell++){
    moab::Range adj;
    mbErr = mbInterface->get_adjacencies(&(*itcell), 1, refElement->getDimension() - 1, 0, adj);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("Mesh", "computeCell2FaceMap", "could not get face adjacencies from cell (moab)"));
    }
    int index = mbInterface->id_from_handle(*itcell) - 1;
    for(moab::Range::iterator itadj = adj.begin(); itadj != adj.end(); itadj++){
      int indexFace = getFaceIndex(mbInterface, *itcell, *itadj, skelRefEl);
      cell2FaceMap[index*nFacesPerCell + indexFace] = mbInterface->id_from_handle(*itadj) - 1;
    }
  }
};//computeCell2FaceMap

const std::vector<int> * Mesh::getFace2CellMap() const{
  return &face2CellMap;
};//getFace2CellMap

const std::vector<int> * Mesh::getCell2FaceMap() const{
  return &cell2FaceMap;
};//getFace2CellMap

void Mesh::getFace2Cell(int i, std::vector<int> * face2Cell) const{
  if(face2CellMap[i*2 + 1] != -1){
    face2Cell->resize(2);
    face2Cell->assign(face2CellMap.begin() + i*2, face2CellMap.begin() + (i+1)*2);
  } else {
    face2Cell->resize(1);
    (*face2Cell)[0] = face2CellMap[i*2];
  }
};//getCell2Face

void Mesh::getCell2Face(int i, std::vector<int> * cell2Face) const{
  cell2Face->resize(nFacesPerCell);
  cell2Face->assign(cell2FaceMap.begin() + i*nFacesPerCell, cell2FaceMap.begin() + (i+1)*nFacesPerCell);
};//getCell2Face

int Mesh::getNumFacesPerCell() const{
  return nFacesPerCell;
};//getNumFacesPerCell

int Mesh::getFaceIndex(moab::Interface * mbInterface, const moab::EntityHandle & cell, const moab::EntityHandle & face, ReferenceElement & skelRefEl) const{
  moab::ErrorCode mbErr;
  std::vector<moab::EntityHandle> cellVerts(skelRefEl.getNumNodes());
  mbErr = mbInterface->get_connectivity(&cell, 1, cellVerts);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "getFaceIndex", "could not get vertices from cell (moab)"));
  }
  std::vector<moab::EntityHandle> faceVerts(skelRefEl.getFaceElement()->getNumNodes());
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
  const std::vector< std::vector<int> > * faceNodes = skelRefEl.getFaceNodes();
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

void Mesh::computeBoundary(moab::Interface * mbInterface, moab::EntityHandle & meshset){
  moab::ErrorCode mbErr;
  moab::Range moabCells, boundary;
  mbErr = mbInterface->get_entities_by_dimension(meshset, refElement->getDimension(), moabCells);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "computeBoundary", "could not get cells from meshset"));
  }
  moab::Skinner mbSkinner(mbInterface);
  mbErr = mbSkinner.find_skin(meshset, moabCells, refElement->getDimension()-1, boundary);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "computeBoundary", "could not find the boundary"));
  }
  for(moab::Range::iterator itbound = boundary.begin(); itbound != boundary.end(); itbound++){
    boundaryFaces.insert(mbInterface->id_from_handle(*itbound) - 1);
  }
};//computeBoundary

std::set<int> * Mesh::getBoundaryFaces(){
  return &boundaryFaces;
};//getBoundaryFaces

void Mesh::computeFaceConnectivity(moab::Interface * mbInterface, moab::EntityHandle & meshset){
  moab::ErrorCode mbErr;
  moab::Range moabFaces;
  mbErr = mbInterface->get_entities_by_dimension(meshset, refElement->getDimension() - 1, moabFaces);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("Mesh", "computeFace2CellMap", "could not get faces from meshset"));
  }
  const std::vector< std::vector<int> > * faceIndexes = refElement->getFaceNodes();
  nFaces = moabFaces.size();
  nNodesPerFace = refElement->getFaceElement()->getNumNodes();
  faces.resize(nFaces * nNodesPerFace);
  int indEl, locFaceIndex;
  for(int i = 0; i < nFaces; i++){
    indEl = face2CellMap[i*2];
    locFaceIndex = std::distance(cell2FaceMap.begin() + indEl*nFacesPerCell, 
        std::find(cell2FaceMap.begin() + indEl*nFacesPerCell, cell2FaceMap.end() + (indEl+1)*nFacesPerCell, i));
    for(int j = 0; j < ((*faceIndexes)[locFaceIndex]).size(); j++){
      faces[i*nNodesPerFace + j] = cells[indEl*nNodesPerCell + ((*faceIndexes)[locFaceIndex])[j]];
    }
  }
};//computeFaceConnectivity

}//hfox
