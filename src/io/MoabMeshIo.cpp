#include "MoabMeshIo.h"

namespace hfox{

MoabMeshIo::MoabMeshIo(){
  myMesh = NULL;
};//empty constructor

MoabMeshIo::MoabMeshIo(Mesh * mesh){
  setMesh(mesh);
};//mesh constructor

void MoabMeshIo::load(std::string filename){
  if(myMesh == NULL){
    throw(ErrorHandle("MoabMeshIo", "load", "must enter a mesh into the io before loading"));
  }
  moab::ErrorCode mbErr;
  moab::Interface * mbIFace = new (std::nothrow) moab::Core;
  if(mbIFace == NULL){
    throw(ErrorHandle("MoabMeshIo", "load", "could not construct new moab mesh."));
  }
  moab::EntityHandle fs;
  mbErr = mbIFace->create_meshset(moab::MESHSET_SET, fs);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("MoabMeshIo", "load", "could not create meshset"));
  }
  mbErr = mbIFace->load_file(filename.c_str(), &fs);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("MoabMeshIo", "load", "could not load mesh file: " + filename));
  }
  setMeshSet(mbIFace, fs);
  delete mbIFace;
};//load

void MoabMeshIo::setMeshSet(moab::Interface * mbIFace, moab::EntityHandle & meshset){
  std::vector<double> points;
  std::vector<int> cells;
  moab::ErrorCode mbErr;
  moab::Range vertexes;
  mbErr = mbIFace->get_entities_by_dimension(meshset, 0, vertexes);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("MoabMeshIo", "setMeshSet", "could not get vertexes from meshset"));
  }
  //points
  int dim = myMesh->getReferenceElement()->getDimension();
  points.resize(vertexes.size()*dim);
  double pV[3];
  for(moab::Range::iterator it = vertexes.begin(); it != vertexes.end(); it++){
    mbErr = mbIFace->get_coords(&(*it), 1, pV);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("MoabMeshIo", "setMeshSet", "could not get vertex coords from list"));
    }
    int index = mbIFace->id_from_handle(*it) -1;
    for(int i = 0; i < dim ; i++){
      points[index*dim + i] = *(pV+i);
    }
  }
  //cells
  moab::Range elems;
  mbErr = mbIFace->get_entities_by_dimension(meshset, dim, elems);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("MoabMeshIo", "setMeshSet", "could not get elements from meshset"));
  }
  int nNodesPerCell = myMesh->getReferenceElement()->getNumNodes();
  cells.resize(elems.size()*nNodesPerCell);
  for(moab::Range::iterator it = elems.begin(); it != elems.end(); it++){
    std::vector<moab::EntityHandle> verts;
    mbErr = mbIFace->get_connectivity(&(*it), 1, verts);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("MoabMeshIo", "setMeshSet", "could not get vertices from element"));
    }
    int index = mbIFace->id_from_handle(*it) - 1;
    for(int i = 0; i < verts.size(); i++){
      cells[index*nNodesPerCell + i] = mbIFace->id_from_handle(verts[i])-1;
    }
  }
  myMesh->setMesh(dim, points, cells);
};//setMeshSet

void MoabMeshIo::write(std::string filename){
  if(myMesh == NULL){
    throw(ErrorHandle("MoabMeshIo", "write", "must enter a mesh into the io before writing"));
  }
  moab::ErrorCode mbErr;
  moab::Interface * mbIFace = new (std::nothrow) moab::Core;
  if(mbIFace == NULL){
    throw(ErrorHandle("MoabMeshIo", "write", "could not construct new moab mesh."));
  }
  moab::EntityHandle meshset;
  mbErr = mbIFace->create_meshset(moab::MESHSET_SET, meshset);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("MoabMeshIo", "write", "could not create meshset"));
  }
  const std::vector<double> * points = myMesh->getPoints();
  int nPoints = myMesh->getNumberPoints();
  int dim = myMesh->getNodeSpaceDimension();
  const std::vector<int> * cells = myMesh->getCells();
  int nCells = myMesh->getNumberCells();
  int nNodesPerCell = myMesh->getReferenceElement()->getNumNodes();
  std::vector<double> point(3, 0.0);
  std::vector<moab::EntityHandle> cell(nNodesPerCell);
  std::vector<moab::EntityHandle> orderedVertList(nPoints);
  std::vector<moab::EntityHandle> orderedElList(nCells);
  for(int i = 0; i < nPoints; i++){
    moab::EntityHandle ent;
    for(int j = 0; j < dim; j++){
      point[j] = (*points)[i*dim + j];
    }
    mbErr = mbIFace->create_vertex(point.data(), ent);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("MoabMeshIo", "write", "could not create vertex (moab)."));
    }
    orderedVertList[i] = ent;
  }
  mbErr = mbIFace->add_entities(meshset, orderedVertList.data(), orderedVertList.size());
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("MoabMeshIo", "write", "could not set vertexes in mesh (moab)."));
  }
  moab::EntityType type = determineMOABType();
  for(int i = 0; i < nCells; i++){
    moab::EntityHandle ent;
    for(int j = 0; j < nNodesPerCell; j++){
      cell[j] = orderedVertList[(*cells)[i*nNodesPerCell + j]];
    }
    mbErr = mbIFace->create_element(type, cell.data(), cell.size(), ent);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("MoabMeshIo", "write", "could not create element (moab)."));
    }
    orderedElList[i] = ent;
  }
  mbErr = mbIFace->add_entities(meshset, orderedElList.data(), orderedElList.size());
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("MoabMeshIo", "write", "could not set elements in mesh (moab)."));
  }

  mbErr = mbIFace->write_file(filename.c_str(), 0, 0, &meshset, 1);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("MoabMeshIo", "write", "could not write mesh file: " + filename));
  }
  delete mbIFace;
};//write

moab::EntityType MoabMeshIo::determineMOABType() const{
  moab::EntityType type;
  const ReferenceElement * refElement = myMesh->getReferenceElement();
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

}//hfox
