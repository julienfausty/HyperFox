#include "MoabMeshIo.h"

namespace hfox{

MoabMeshIo::MoabMeshIo() : myMesh(NULL), mbIFace(NULL){
};//empty constructor

MoabMeshIo::MoabMeshIo(Mesh * mesh) : myMesh(NULL), mbIFace(NULL){
  setMesh(mesh);
};//mesh constructor

void MoabMeshIo::setMesh(Mesh * mesh){
  if(mesh->getReferenceElement() == NULL){
    throw(ErrorHandle("MoabMeshIo", "setMesh", "reference element of mesh is not defined."));
  }
  myMesh = mesh;
  mbIFace = mesh->getMoabInterface();
};//setMesh

void MoabMeshIo::load(std::string filename){
  if(myMesh == NULL){
    throw(ErrorHandle("MoabMeshIo", "load", "must enter a mesh into the io before loading"));
  }
  moab::ErrorCode mbErr;
  moab::EntityHandle fs;
  mbErr = mbIFace->create_meshset(moab::MESHSET_SET, fs);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("MoabMeshIo", "load", "could not create meshset"));
  }
  myMesh->setMeshSet(fs);
  mbErr = mbIFace->load_file(filename.c_str(), &fs);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("MoabMeshIo", "load", "could not load mesh file: " + filename));
  }
  myMesh->setMeshSet(fs);
};//load

void MoabMeshIo::write(std::string filename){
  if(myMesh == NULL){
    throw(ErrorHandle("MoabMeshIo", "write", "must enter a mesh into the io before writing"));
  }
  moab::ErrorCode mbErr;
  mbErr = mbIFace->write_file(filename.c_str(), 0, 0, myMesh->getMeshSet(), 1);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("MoabMeshIo", "write", "could not write mesh file: " + filename));
  }
};//write

}//hfox
