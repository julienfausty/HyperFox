#include "ZoltanPartitioner.h"

namespace hfox{

ZoltanPartitioner::~ZoltanPartitioner(){
  if(!(zChanges.isfree)){
    freeZChanges();
  }
  if(zObj != NULL){
    delete zObj;
  }
};//destructor

void ZoltanPartitioner::initialize(){
  Partitioner::initialize();
  float version;
  Zoltan_Initialize(myOpts.argc, myOpts.argv, &version);
  if(zObj != NULL){
    delete zObj;
  }
  zObj = new Zoltan(MPI_COMM_WORLD);
  zObj->Set_Param("LB_METHOD", myOpts.method);
  zObj->Set_Param("GRAPH_PACKAGE", myOpts.graphPackage);
  zObj->Set_Param("LB_APPROACH", myOpts.approach);
  zObj->Set_Param("NUM_GID_ENTRIES", "1");
  zObj->Set_Param("NUM_LID_ENTRIES", "1");
  zObj->Set_Param("DEBUG_LEVEL", myOpts.debugLevel);
  zObj->Set_Param("OBJ_WEIGHT_DIM", "0");
  zObj->Set_Param("EDGE_WEIGHT_DIM", "0");
  zObj->Set_Param("TFLOPS_SPECIAL", "1");
  zChanges.isfree = 1;
};//initialize


void ZoltanPartitioner::computePartition(){
  if(!initialized){
    throw(ErrorHandle("ZoltanPartitioner", "computePartition", "must initialize the Zoltan Partitioner before computing the partition"));
  }
  if((myOpts.debugLevel != "0") and (rank == 0)){
    std::cout << "Zoltan Configuration: " << std::endl;
  }
  zObj->Set_Num_Obj_Fn(ZoltanPartitioner::getNumberCells, myMesh);
  zObj->Set_Obj_List_Fn(ZoltanPartitioner::getListCells, &elementIDs);
  std::pair<Mesh*, Partitioner*> thePair(myMesh, this);
  zObj->Set_Num_Edges_Multi_Fn(ZoltanPartitioner::getNumNeighborsCellSlice, &thePair);
  zObj->Set_Edge_List_Multi_Fn(ZoltanPartitioner::getListNeighborsCellSlice, &thePair);
  //run partition
  if(myOpts.debugLevel != "0" and rank == 0){
    std::cout << "Running Zoltan... " << std::endl;
  }
  if(!(zChanges.isfree)){
    freeZChanges();
  }
  int ierr = zObj->LB_Partition(zChanges.changes, zChanges.num_gid_entries, zChanges.num_lid_entries, zChanges.num_import,
      zChanges.import_global_ids, zChanges.import_local_ids, zChanges.import_procs, zChanges.import_to_part, 
      zChanges.num_export, zChanges.export_global_ids, zChanges.export_local_ids, zChanges.export_procs, zChanges.export_to_part);
  zChanges.isfree = 0;
  if(ierr != ZOLTAN_OK){
    throw(ErrorHandle("ZoltanPartitioner", "computePartition", "Zoltan failed to compute partition, might want to rerun with higher debug level in options to see what's wrong."));
  }
  //format the cell partition
  emptyImportExportMaps();
  for(int i = 0; i < zChanges.num_export; i++){
    if(exportMap.find(zChanges.export_procs[i]) == exportMap.end()){
      exportMap[zChanges.export_procs[i]] = {{Cell, std::vector<int>()}};
    }
    exportMap[zChanges.export_procs[i]][Cell].push_back(zChanges.export_global_ids[i]);
  }
  for(int i = 0; i < zChanges.num_import; i++){
    if(importMap.find(zChanges.import_procs[i]) == importMap.end()){
      importMap[zChanges.import_procs[i]] = {{Cell, std::vector<int>()}};
    }
    importMap[zChanges.import_procs[i]][Cell].push_back(zChanges.import_global_ids[i]);
  }
  freeZChanges();
  std::map<int, std::map<FieldType, std::vector<int> > >::iterator itMap;
  std::vector<int> cell;
  std::vector<int> cell2Face;
  int localCellIndex;
  if(myOpts.debugLevel != "0" and rank == 0){
    std::cout << "Partitioning faces and nodes based on cell partition" << std::endl;
  }
  //decide how to partition the nodes and faces based on the cell partition
  for(itMap = exportMap.begin(); itMap != exportMap.end(); itMap++){
    std::set<int> exportCandFaces;
    std::set<int> exportCandNodes;
    for(int i = 0; i < (itMap->second[Cell]).size(); i++){
      //cells
      localCellIndex = global2LocalElement((itMap->second)[Cell][i]);
      myMesh->getCell(localCellIndex, &cell);
      myMesh->getCell2Face(localCellIndex, &cell2Face);
      for(int k = 0; k < cell2Face.size(); k++){
        exportCandFaces.insert(cell2Face[k]);
      }
      for(int k = 0; k < cell.size(); k++){
        exportCandNodes.insert(cell[k]);
      }
    }
    std::set<int>::iterator itset;
    itMap->second[Face] = std::vector<int>();
    itMap->second[Node] = std::vector<int>();
    std::map<int, std::map<FieldType, std::vector<int> > >::iterator itBuff;
    for(itset = exportCandFaces.begin(); itset != exportCandFaces.end(); itset++){
      localCellIndex = global2LocalFace(*itset);
      if(localCellIndex != -1){
        bool exported = 0;
        for(itBuff = exportMap.begin(); itBuff != itMap; itBuff++){
          exported = (std::find(itBuff->second[Face].begin(), itBuff->second[Face].end(), *itset) != itBuff->second[Face].end());
          if(exported){
            break;
          }
        }
        if(!exported){
          itMap->second[Face].push_back(*itset);
        }
      }
    }
    for(itset = exportCandNodes.begin(); itset != exportCandNodes.end(); itset++){
      localCellIndex = global2LocalNode(*itset);
      if(localCellIndex != -1){
        bool exported = 0;
        for(itBuff = exportMap.begin(); itBuff != itMap; itBuff++){
          exported = (std::find(itBuff->second[Node].begin(), itBuff->second[Node].end(), *itset) != itBuff->second[Node].end());
          if(exported){
            break;
          }
        }
        if(!exported){
          itMap->second[Node].push_back(*itset);
        }
      }
    }
  }
  if(myOpts.debugLevel != "0" and rank == 0){
    std::cout << "Communicating face and node partitions" << std::endl;
  }
  int tag;
  //communicate the node and face partition to relevant processes
  for(itMap = exportMap.begin(); itMap != exportMap.end(); itMap++){
    tag = 0;
    int s = itMap->second.at(Face).size();
    MPI_Send(&s, 1, MPI_INT, itMap->first, tag, MPI_COMM_WORLD);
    tag += 1;
    s = itMap->second.at(Node).size();
    MPI_Send(&s, 1, MPI_INT, itMap->first, tag, MPI_COMM_WORLD);
    tag += 1;
    MPI_Send(itMap->second.at(Face).data(), itMap->second.at(Face).size(), MPI_INT, itMap->first, tag, MPI_COMM_WORLD);
    tag += 1;
    MPI_Send(itMap->second.at(Node).data(), itMap->second.at(Node).size(), MPI_INT, itMap->first, tag, MPI_COMM_WORLD);
  }
  for(itMap = importMap.begin(); itMap != importMap.end(); itMap++){
    tag = 0;
    int s;
    MPI_Recv(&s, 1, MPI_INT, itMap->first, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    itMap->second[Face] = std::vector<int>(s, 0);
    tag += 1;
    MPI_Recv(&s, 1, MPI_INT, itMap->first, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    itMap->second[Node] = std::vector<int>(s, 0);
    tag += 1;
    MPI_Recv(itMap->second.at(Face).data(), itMap->second.at(Face).size(), MPI_INT, itMap->first, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    tag += 1;
    MPI_Recv(itMap->second.at(Node).data(), itMap->second.at(Node).size(), MPI_INT, itMap->first, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  if(myOpts.debugLevel != "0" and rank == 0){
    std::cout << "finished partitioning with Zoltan" << std::endl;
  }
};//computePartition


int ZoltanPartitioner::getNumberCells(void * data, int * ierr){
  Mesh * thisMesh = (Mesh*) data;
  *ierr = ZOLTAN_OK;
  return thisMesh->getNumberCells();
};//getNumberCells

void ZoltanPartitioner::getListCells(void * data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim, float * obj_weights, int * ierr){
  std::vector<unsigned int> * elIDs = (std::vector<unsigned int>*) data;
  *ierr = ZOLTAN_OK;
  std::copy(elIDs->begin(), elIDs->end(), global_ids);
  std::vector<unsigned int> locIds(elIDs->size());
  std::iota(locIds.begin(), locIds.end(), 0);
  std::copy(locIds.begin(), locIds.end(), local_ids);
};//getListCells

void ZoltanPartitioner::getNumNeighborsCellSlice(void * data, int num_gid_entries, int num_lid_entries, int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int * num_edges, int * ierr){
  std::pair<Mesh*, Partitioner*> * myPair = (std::pair<Mesh*, Partitioner*> *) data;
  *ierr = ZOLTAN_OK;
  std::vector<int> cell2Face(myPair->first->getReferenceElement()->getNumFaces());
  std::vector<int> face2Cell;
  int numNeighbors;
  int locFace;
  for(int i = 0; i < num_obj; i++){
    numNeighbors = 0;
    myPair->first->getCell2Face(local_ids[i], &cell2Face);
    for(int j = 0; j < cell2Face.size(); j++){
      locFace = myPair->second->global2LocalFace(cell2Face[j]);
      if(locFace != -1){
        myPair->first->getFace2Cell(locFace, &face2Cell);
        numNeighbors += face2Cell.size() - 1;
      } else{
        numNeighbors += 1;
      }
    }
    num_edges[i] = numNeighbors;
  }
};//getNumNeighborsCellSlice

void ZoltanPartitioner::getListNeighborsCellSlice(void * data, int num_gid_entries, int num_lid_entries, int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int * num_edges, ZOLTAN_ID_PTR nbor_global_ids, int * nbor_procs, int wgt_dim, float * ewgts, int * ierr){
  std::pair<Mesh*, Partitioner*> * myPair = (std::pair<Mesh*, Partitioner*> *) data;
  *ierr = ZOLTAN_OK;
  std::vector<int> cell2Face(myPair->first->getReferenceElement()->getNumFaces());
  std::vector<int> face2Cell;
  int locFace;
  int posNbors = 0;
  for(int i = 0; i < num_obj; i++){// for the cells
    myPair->first->getCell2Face(local_ids[i], &cell2Face);
    for(int j = 0; j < cell2Face.size(); j++){ //for the faces on the cells
      locFace = myPair->second->global2LocalFace(cell2Face[j]);// get the local index of the face
      if(locFace != -1){//if the face is held by this partition
        myPair->first->getFace2Cell(locFace, &face2Cell);
        if(face2Cell.size() != 1){//if the face is not part of the boundary of the mesh
          nbor_procs[posNbors] = myPair->second->getRank();
          if(face2Cell[0] == global_ids[i]){//if the current cell is the first adjacency
            nbor_global_ids[posNbors] = face2Cell[1];//take the second adjacency
          } else if(face2Cell[1] == global_ids[i]){//if the current cell is the second adjacency
            nbor_global_ids[posNbors] = face2Cell[0];//take the first
          } else {
            throw(ErrorHandle("ZoltanPartitioner", "getListNeighborsCellSlice", "big issue: one of the faces of the element does not have the element as an adjacency in the face2Cell map"));
          }
          posNbors += 1;
        }
      }
      else{//if the face is in fact shared and on another partition
        bool foundEl = 0;
        for(int k = 0; k < myPair->second->getSharedFaceList()->size()/3; k++){//for all the shared faces
          if(myPair->second->getSharedFaceList()->at(k*3) == cell2Face[j]){//if the shared face is the face we are looking for
            foundEl = 1;
            nbor_procs[posNbors] = myPair->second->getSharedFaceList()->at(k*3 + 1);
            nbor_global_ids[posNbors] = myPair->second->getSharedFaceList()->at(k*3 + 2);
            posNbors += 1;
            break;
          }
        }
        if(!foundEl){
          throw(ErrorHandle("ZoltanPartitioner", "getListNeighborsCellSlice", "could not find shared face in the sharedFaceList"));
        }
      }
    }
  }
};//getListNeighborsCellSlice

void ZoltanPartitioner::freeZChanges(){
  zObj->LB_Free_Part(&(zChanges.import_global_ids), &(zChanges.import_local_ids), &(zChanges.import_procs), &(zChanges.import_to_part));
  zObj->LB_Free_Part(&(zChanges.export_global_ids), &(zChanges.export_local_ids), &(zChanges.export_procs), &(zChanges.export_to_part));
  zChanges.isfree = 1;
};//freeZChanges

};//hfox
