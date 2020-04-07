#include "HDF5Io.h"

namespace hfox{

HDF5Io::HDF5Io(){
  myMesh = NULL;
};//empty constructor

HDF5Io::HDF5Io(Mesh * mesh){
  setMesh(mesh);
};//mesh constructor

void HDF5Io::load(std::string filename){
  if(myMesh == NULL){
    throw(ErrorHandle("HDF5Io", "load", "the mesh must be set before loading a file."));
  }
  if(!H5Fis_hdf5(filename.c_str())){
    throw(ErrorHandle("HDF5Io", "load", "the file " + filename + " is not an hdf5 file."));
  }
  hid_t file;
  hid_t tempid;
  file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  bool meshExists, fieldExists;
  meshExists = H5Lexists(file, "Mesh", tempid);
  fieldExists = H5Lexists(file, "FieldData", tempid);
  if(meshExists){
    hid_t meshGrp = H5Gopen(file, "Mesh", tempid);
    loadMesh(meshGrp);
    H5Gclose(meshGrp);
  }
  if(fieldExists){
    hid_t fieldGrp = H5Gopen(file, "FieldData", tempid);
    loadFields(fieldGrp);
    H5Gclose(fieldGrp);
  }
  if((!meshExists) and (!fieldExists)){
    throw(ErrorHandle("HDFIo", "load", "could not find Mesh or FieldData groups in file"));
  }
  H5Fclose(file);
};//load

void HDF5Io::write(std::string filename){
  hid_t file;
  hid_t lcplid, gcplid, gaplid;
  file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  bool meshExists, fieldExists;
  meshExists = (myMesh != NULL);
  fieldExists = (!(fieldMap.empty()));
  if(meshExists){
    hid_t meshGrp = H5Gcreate(file, "Mesh", lcplid, gcplid, gaplid);
    writeMesh(meshGrp);
    H5Gclose(meshGrp);
  }
  if(fieldExists){
    hid_t fieldGrp = H5Gcreate(file, "FieldData", lcplid, gcplid, gaplid);
    writeFields(fieldGrp);
    H5Gclose(fieldGrp);
  }
  if((!meshExists) and (!fieldExists)){
    throw(ErrorHandle("HDFIo", "write", "could not find anything to write"));
  }
  H5Fclose(file);
};//write

void HDF5Io::loadMesh(hid_t meshGrp){
  //std::vector< std::vector<double> > nodes;
  //std::vector< std::vector<int> > cells;
  //hid_t dataspace;
  //std::vector<hsize_t> shape(2);
  //int status;
  ////nodes
  //hid_t nodeSet = H5Dopen(meshGrp, "Nodes", H5P_DEFAULT);
  //dataspace = H5Dget_space(nodeSet);
  //status = H5Dget_simple_extent_dims(dataspace, shape.data(), NULL);
  //if(status < 0){
    //throw(ErrorHandle("HDF5Io", "loadMesh", "could not get shape of dataspace"));
  //}
};//loadMesh

void HDF5Io::loadFields(hid_t fieldGrp){
};//loadField

void HDF5Io::writeMesh(hid_t meshGrp){
};//writeMesh

void HDF5Io::writeFields(hid_t fieldGrp){
};//writeField

}
