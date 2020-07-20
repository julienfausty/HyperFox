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
    throw(ErrorHandle("HDF5Io", "load", "must enter a mesh into the io before loading a file."));
  }
  bool is_master = 1;
  int mpiInit = 0;
  int mpiRank;
  int dimNodeSpace = 0;
  MPI_Initialized(&mpiInit);
  if(mpiInit){
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    is_master = (mpiRank == 0);
  }
  if(is_master){
    if(!(H5Fis_hdf5(filename.c_str()) > 0)){
      throw(ErrorHandle("HDF5Io", "load", "the file " + filename + " is not an hdf5 file."));
    }
    hid_t file;
    hid_t tempid;
    file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    bool meshExists, fieldExists;
    herr_t status;
    status = H5Lexists(file, "Mesh", H5P_DEFAULT);
    meshExists = (status > 0);
    status = H5Lexists(file, "FieldData", H5P_DEFAULT);
    fieldExists = (status > 0);
    if(meshExists){
      hid_t meshGrp = H5Gopen(file, "Mesh", H5P_DEFAULT);
      loadMesh(meshGrp);
      H5Gclose(meshGrp);
    }
    if(fieldExists){
      hid_t fieldGrp = H5Gopen(file, "FieldData", H5P_DEFAULT);
      loadFields(fieldGrp);
      H5Gclose(fieldGrp);
    }
    if(!(meshExists || fieldExists)){
      throw(ErrorHandle("HDFIo", "load", "could not find Mesh or FieldData groups in file"));
    }
    H5Fclose(file);
    dimNodeSpace = myMesh->getNodeSpaceDimension();
  }
  if(mpiInit){
    MPI_Bcast(&dimNodeSpace, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(!is_master){
      std::vector<double> emptynodes;
      std::vector<int> emptycells;
      myMesh->setMesh(dimNodeSpace, emptynodes, emptycells);
      myMesh->update();
    }
  }
};//load

void HDF5Io::write(std::string filename){
  hid_t file;
  hid_t meshGrp;
  hid_t fieldGrp;
  bool is_master = 1;
  int mpiInit = 0;
  int mpiRank;
  int dimNodeSpace = 0;
  MPI_Initialized(&mpiInit);
  if(mpiInit){
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    is_master = (mpiRank == 0);
  }
  if(is_master){
    file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }
  bool meshExists, fieldExists;
  meshExists = (myMesh != NULL);
  fieldExists = (!(fieldMap.empty()));
  if(meshExists){
    if(is_master){
      meshGrp = H5Gcreate(file, "Mesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    writeMesh(meshGrp);
    if(is_master){
      H5Gclose(meshGrp);
    }
  }
  if(fieldExists){
    if(is_master){
      fieldGrp = H5Gcreate(file, "FieldData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    writeFields(fieldGrp);
    if(is_master){
      H5Gclose(fieldGrp);
    }
  }
  if((!meshExists) and (!fieldExists)){
    throw(ErrorHandle("HDFIo", "write", "could not find anything to write"));
  }
  if(is_master){
    H5Fclose(file);
  }
};//write

void HDF5Io::loadMesh(hid_t meshGrp){
  int dimNodeSpace;
  std::vector<double> nodes;
  std::vector<int> cells;
  hid_t dataspace, datatype;
  std::vector<hsize_t> shape(2);
  int status;
  //nodes
  hid_t nodeSet = H5Dopen(meshGrp, "Nodes", H5P_DEFAULT);
  datatype = H5Dget_type(nodeSet);
  dataspace = H5Dget_space(nodeSet); 
  status = H5Sget_simple_extent_dims(dataspace, shape.data(), NULL);
  if(status < 0){
    throw(ErrorHandle("HDF5Io", "loadMesh", "could not get shape of node dataspace"));
  }
  dimNodeSpace = shape[1];
  nodes.resize(shape[0]*shape[1]);
  status = H5Dread(nodeSet, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, nodes.data());
  if(status < 0){
    throw(ErrorHandle("HDF5Io", "loadMesh", "could not load nodes into buffer"));
  }
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Dclose(nodeSet);
  //cells
  hid_t cellSet = H5Dopen(meshGrp, "Cells", H5P_DEFAULT);
  datatype = H5Dget_type(cellSet);
  dataspace = H5Dget_space(cellSet);
  status = H5Sget_simple_extent_dims(dataspace, shape.data(), NULL);
  if(status < 0){
    throw(ErrorHandle("HDF5Io", "loadMesh", "could not get shape of cell dataspace"));
  }
  cells.resize(shape[0]*shape[1]);
  status = H5Dread(cellSet, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cells.data());
  if(status < 0){
    throw(ErrorHandle("HDF5Io", "loadMesh", "could not load cells into buffer"));
  }
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Dclose(cellSet);
  myMesh->setMesh(dimNodeSpace, nodes, cells);
};//loadMesh

void HDF5Io::loadFields(hid_t fieldGrp){
  herr_t status;
  hid_t datatype, dataspace, fieldSet, fType;
  std::vector<hsize_t> shape(3);
  for(std::map<std::string, Field*>::iterator itfMap = fieldMap.begin(); itfMap != fieldMap.end(); itfMap++){
    status = H5Lexists(fieldGrp, itfMap->first.c_str(), H5P_DEFAULT);
    if(status < 0){
      throw(ErrorHandle("HDF5Io", "loadFields", "field with name " + itfMap->first + " was not found in FieldData"));
    }
    fieldSet = H5Dopen(fieldGrp, itfMap->first.c_str(), H5P_DEFAULT);
    datatype = H5Dget_type(fieldSet);
    dataspace = H5Dget_space(fieldSet); 
    status = H5Sget_simple_extent_dims(dataspace, shape.data(), NULL);
    if(status < 0){
      throw(ErrorHandle("HDF5Io", "loadFields", "could not get shape of field dataspace"));
    }
    status = H5Aexists(fieldSet, "ftype");
    if(!(status > 0)){
      throw(ErrorHandle("HDFIo", "loadFields", "attribute ftype could not be found for field: " + itfMap->first));
    }
    fType = H5Aopen_by_name(fieldSet, ".", "ftype", H5P_DEFAULT, H5P_DEFAULT);
    int ibuff;
    status = H5Aread(fType, H5T_NATIVE_INT, &ibuff);
    *(itfMap->second->getFieldType()) = (FieldType)ibuff;
    *(itfMap->second->getNumEntities()) = shape[0];
    *(itfMap->second->getNumObjPerEnt()) = shape[1];
    *(itfMap->second->getNumValsPerObj()) = shape[2];
    itfMap->second->allocate();
    status = H5Dread(fieldSet, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, itfMap->second->getValues()->data());
    if(status < 0){
      throw(ErrorHandle("HDF5Io", "loadMesh", "could not load field values into field"));
    }
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Aclose(fType);
    H5Dclose(fieldSet);
  }
};//loadField

void HDF5Io::writeMesh(hid_t meshGrp){
  Partitioner * part = myMesh->getPartitioner();
  std::vector<int> cells;
  std::vector<double> nodes;
  const std::vector<int> * meshCells = myMesh->getCells();
  const std::vector<double> * meshNodes = myMesh->getPoints();
  int nNodesPCell = myMesh->getReferenceElement()->getNumNodes();
  int dimSpace = myMesh->getNodeSpaceDimension();
  int nNodes;
  int nCells;
  if(part == NULL){
    nNodes = myMesh->getNumberPoints();
    nCells = myMesh->getNumberCells();
    cells.resize(meshCells->size());
    std::copy(meshCells->begin(), meshCells->end(), cells.begin());
    nodes.resize(meshNodes->size());
    std::copy(meshNodes->begin(), meshNodes->end(), nodes.begin());
  } else {
    std::vector<int> nodeIds;
    std::vector<int> cellIds;
    nNodes = part->getTotalNumberNodes();
    nCells = part->getTotalNumberEls();
    if(part->getRank() == 0){
      nodeIds.resize(nNodes, 0);
      nodes.resize(nNodes * (myMesh->getNodeSpaceDimension()));
      cellIds.resize(nCells, 0);
      cells.resize(nCells * nNodesPCell);
    }
    std::vector<int> nObjPPart(part->getNumPartitions());
    std::vector<int> offsets(part->getNumPartitions());
    int nLocObj = myMesh->getNumberPoints();
    MPI_Allgather(&nLocObj, 1, MPI_INT, nObjPPart.data(), 1, MPI_INT, MPI_COMM_WORLD);
    for(int i = 0; i < part->getNumPartitions(); i++){
      offsets[i] = std::accumulate(nObjPPart.begin(), nObjPPart.begin() + i, 0);
    }
    MPI_Gatherv(part->getNodeIds()->data(), nLocObj, MPI_INT, nodeIds.data(), nObjPPart.data(), offsets.data(), MPI_INT, 0, MPI_COMM_WORLD);
    for(int i = 0; i < part->getNumPartitions(); i++){
      offsets[i] *= myMesh->getNodeSpaceDimension();
      nObjPPart[i] *= myMesh->getNodeSpaceDimension();
    }
    nLocObj *= myMesh->getNodeSpaceDimension();
    MPI_Gatherv(myMesh->getPoints()->data(), nLocObj, MPI_DOUBLE, nodes.data(), nObjPPart.data(), offsets.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    nLocObj = myMesh->getNumberCells();
    MPI_Allgather(&nLocObj, 1, MPI_INT, nObjPPart.data(), 1, MPI_INT, MPI_COMM_WORLD);
    for(int i = 0; i < part->getNumPartitions(); i++){
      offsets[i] = std::accumulate(nObjPPart.begin(), nObjPPart.begin() + i, 0);
    }
    MPI_Gatherv(part->getCellIds()->data(), nLocObj, MPI_INT, cellIds.data(), nObjPPart.data(), offsets.data(), MPI_INT, 0, MPI_COMM_WORLD);
    for(int i = 0; i < part->getNumPartitions(); i++){
      offsets[i] *= nNodesPCell;
      nObjPPart[i] *= nNodesPCell;
    }
    nLocObj *= nNodesPCell;
    MPI_Gatherv(myMesh->getCells()->data(), nLocObj, MPI_INT, cells.data(), nObjPPart.data(), offsets.data(), MPI_INT, 0, MPI_COMM_WORLD);
    if(part->getRank() == 0){
      std::vector<double> dbuffer(nodes.size(), 0.0);
      for(int i = 0; i < nNodes; i++){
        for(int k = 0; k < dimSpace; k++){
          dbuffer[nodeIds[i]*dimSpace + k] = nodes[i*dimSpace + k];
        }
      }
      std::copy(dbuffer.begin(), dbuffer.end(), nodes.begin());
      dbuffer.clear();
      std::vector<int> ibuffer(cells.size(), 0);
      for(int i = 0; i < nCells; i++){
        for(int k = 0; k < nNodesPCell; k++){
          ibuffer[cellIds[i]*nNodesPCell + k] = cells[i*nNodesPCell + k];
        }
      }
      std::copy(ibuffer.begin(), ibuffer.end(), cells.begin());
      ibuffer.clear();
    }
  }
  bool is_master = 1;
  int mpiInit = 0;
  int mpiRank;
  int dimNodeSpace = 0;
  MPI_Initialized(&mpiInit);
  if(mpiInit){
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    is_master = (mpiRank == 0);
  }
  if(is_master){
    hid_t dataspace, nodeSet, cellSet;
    herr_t status;
    std::vector<hsize_t> shape(2);
    //nodes
    shape[0] = nNodes;
    shape[1] = myMesh->getNodeSpaceDimension();
    dataspace = H5Screate_simple(2, shape.data(), NULL);
    nodeSet = H5Dcreate(meshGrp, "Nodes", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(nodeSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nodes.data());
    if(status < 0){
      throw(ErrorHandle("HDF5Io", "writeMesh", "problem writing nodes of mesh"));
    }
    H5Sclose(dataspace);
    H5Dclose(nodeSet);
    //cells
    shape[0] = nCells;
    shape[1] = nNodesPCell;
    dataspace = H5Screate_simple(2, shape.data(), NULL);
    cellSet = H5Dcreate(meshGrp, "Cells", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(cellSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cells.data());
    if(status < 0){
      throw(ErrorHandle("HDF5Io", "writeMesh", "problem writing cells of mesh"));
    }
    H5Sclose(dataspace);
    H5Dclose(cellSet);
  }
};//writeMesh

void HDF5Io::writeFields(hid_t fieldGrp){
  herr_t status;
  hid_t dataspace, fieldSet, fType, attrDataspace;
  std::vector<hsize_t> shape(3);
  int nEntities;
  std::vector<double> vals;
  int nFVals;
  for(std::map<std::string, Field*>::iterator it = fieldMap.begin(); it != fieldMap.end(); it++){
    Mesh * fmesh = it->second->getMesh();
    Partitioner * part = fmesh->getPartitioner();
    nFVals = (*(it->second->getNumObjPerEnt())) * (*(it->second->getNumValsPerObj()));
    if(part == NULL){
      nEntities = *(it->second->getNumEntities());
      vals.resize(it->second->getValues()->size());
      std::copy(it->second->getValues()->begin(), it->second->getValues()->end(), vals.begin());
    } else {
      std::vector<int> fIds;
      const std::vector<int> * ids;
      switch(*(it->second->getFieldType())){
        case None: {nEntities = 0; break;}
        case Node: {nEntities = part->getTotalNumberNodes(); ids = part->getNodeIds(); break;}
        case Edge: {throw(ErrorHandle("HDF5Io", "writeFields", "edge fields are not supported yet.")); break;}
        case Face: {nEntities = part->getTotalNumberFaces(); ids = part->getFaceIds(); break;}
        case Cell: {nEntities = part->getTotalNumberEls(); ids = part->getCellIds(); break;}
      }
      if(part->getRank() == 0){
        fIds.resize(nEntities, 0);
        vals.resize(nEntities * nFVals);
      }
      std::vector<int> nObjPPart(part->getNumPartitions());
      std::vector<int> offsets(part->getNumPartitions());
      int nLocObj = *(it->second->getNumEntities());
      MPI_Allgather(&nLocObj, 1, MPI_INT, nObjPPart.data(), 1, MPI_INT, MPI_COMM_WORLD);
      for(int i = 0; i < part->getNumPartitions(); i++){
        offsets[i] = std::accumulate(nObjPPart.begin(), nObjPPart.begin() + i, 0);
      }
      MPI_Gatherv(ids->data(), nLocObj, MPI_INT, fIds.data(), nObjPPart.data(), offsets.data(), MPI_INT, 0, MPI_COMM_WORLD);
      for(int i = 0; i < part->getNumPartitions(); i++){
        offsets[i] *= nFVals;
        nObjPPart[i] *= nFVals;
      }
      nLocObj *= nFVals;
      MPI_Gatherv(it->second->getValues()->data(), nLocObj, MPI_DOUBLE, vals.data(), nObjPPart.data(), offsets.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if(part->getRank() == 0){
        std::vector<double> dbuffer(vals.size(), 0.0);
        for(int i = 0; i < nEntities; i++){
          for(int k = 0; k < nFVals; k++){
            dbuffer[fIds[i]*nFVals + k] = vals[i*nFVals + k];
          }
        }
        std::copy(dbuffer.begin(), dbuffer.end(), vals.begin());
        dbuffer.clear();
      }
    }
    bool is_master = 1;
    int mpiInit = 0;
    int mpiRank;
    int dimNodeSpace = 0;
    MPI_Initialized(&mpiInit);
    if(mpiInit){
      MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
      is_master = (mpiRank == 0);
    }
    if(is_master){
      shape[0] = nEntities;
      shape[1] = *(it->second->getNumObjPerEnt());
      shape[2] = *(it->second->getNumValsPerObj());
      dataspace = H5Screate_simple(3, shape.data(), NULL);
      hsize_t attrDim = 1;
      attrDataspace = H5Screate_simple(1, &attrDim, NULL);
      fieldSet = H5Dcreate(fieldGrp, it->first.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      fType = H5Acreate(fieldSet, "ftype", H5T_NATIVE_INT, attrDataspace, H5P_DEFAULT, H5P_DEFAULT);
      int ibuff = (int)*(it->second->getFieldType());
      status = H5Awrite(fType, H5T_NATIVE_INT, &ibuff);
      if(status < 0){
        throw(ErrorHandle("HDF5Io", "writeFields", "problem writing ftype attribute for field: " + it->first));
      }
      status = H5Dwrite(fieldSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vals.data());
      if(status < 0){
        throw(ErrorHandle("HDF5Io", "writeFields", "problem writing values of field: " + it->first));
      }
      H5Sclose(attrDataspace);
      H5Aclose(fType);
      H5Sclose(dataspace);
      H5Dclose(fieldSet);
    }
  }
};//writeField

}
