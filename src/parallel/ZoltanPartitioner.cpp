#include "ZoltanPartitioner.h"

namespace hfox{

ZoltanPartitioner::~ZoltanPartitioner(){
  MPI_Barrier(MPI_COMM_WORLD);
  //std::cout << "Destroying" << std::endl;
  if(!(zChanges.isfree)){
    //std::cout << "trying to free zChanges" << std::endl;
    freeZChanges();
  }
  if(zObj != NULL){
    //std::cout << "trying to free zObj" <<std::endl;
    delete zObj;
  }
  MPI_Barrier(MPI_COMM_WORLD);
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
  zChanges.isfree = 1;
  MPI_Barrier(MPI_COMM_WORLD);
};//initialize


void ZoltanPartitioner::computePartition(){
  if(!initialized){
    throw(ErrorHandle("ZoltanPartitioner", "computePartition", "must initialize the Zoltan Partitioner before computing the partition"));
  }
  if((myOpts.debugLevel != "0") and (rank == 0)){
    std::cout << "Zoltan Configuration: " << std::endl;
  }
  zObj->Set_Param("OBJ_WEIGHT_DIM", "0");
  zObj->Set_Param("EDGE_WEIGHT_DIM", "0");
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
  if(myOpts.debugLevel != "0" and rank == 0){
    std::cout << "finished partitioning with Zoltan" << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
};//computePartition

void ZoltanPartitioner::update(){
  if(!initialized){
    throw(ErrorHandle("ZoltanPartitioner", "update", "must compute initialize the Partitioner and compute the Partition before updating the fields"));
  }
  if(zChanges.isfree){
    throw(ErrorHandle("ZoltanPartitioner", "update", "must compute the partition before updating the mesh and fields"));
  }
  std::map<int, std::map< FieldType, std::vector<int> > > exportMap, importMap;
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
  int nNodesPCell = myMesh->getReferenceElement()->getNumNodes();
  int nNodesPFace = myMesh->getReferenceElement()->getFaceElement()->getNumNodes();
  int localCellIndex;
  for(itMap = exportMap.begin(); itMap != exportMap.end(); itMap++){
    std::set<int> exportCandFaces;
    std::set<int> exportCandNodes;
    for(int i = 0; i < (itMap->second).size(); i++){
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
    for(itset = exportCandFaces.begin(); itset != exportCandFaces.end(); itset++){
      localCellIndex = global2LocalFace(*itset);
      if(localCellIndex != -1){
        itMap->second[Face].push_back(*itset);
      }
    }
    for(itset = exportCandNodes.begin(); itset != exportCandNodes.end(); itset++){
      localCellIndex = global2LocalNode(*itset);
      if(localCellIndex != -1){
        itMap->second[Node].push_back(*itset);
      }
    }
  }
  //TO DO:
  //make fieldMap (FieldType, *Field)
  //generate buffers to receive the sent items
  //send in order: the cells and the cell2faces the cell fields, the faces and the faces2cells and the boundary and the face fields, the nodes and the node fields
  //recieve the things sent from the procs in the importMap and buffers
  //make modifications to mesh and fields
  //compute the shared faces list
  //update the mesh/fields with shared information
  MPI_Barrier(MPI_COMM_WORLD);
};//update

int ZoltanPartitioner::getNumberCells(void * data, int * ierr){
  //std::cout << "Getting number of cells... " << std::endl;
  Mesh * thisMesh = (Mesh*) data;
  *ierr = ZOLTAN_OK;
  //std::cout << "nCells = " << thisMesh->getNumberCells() << std::endl;
  return thisMesh->getNumberCells();
};//getNumberCells

void ZoltanPartitioner::getListCells(void * data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim, float * obj_weights, int * ierr){
  //std::cout << "Getting list of cells..." << std::endl;
  std::vector<unsigned int> * elIDs = (std::vector<unsigned int>*) data;
  *ierr = ZOLTAN_OK;
  std::copy(elIDs->begin(), elIDs->end(), global_ids);
  std::vector<unsigned int> locIds(elIDs->size());
  std::iota(locIds.begin(), locIds.end(), 0);
  std::copy(locIds.begin(), locIds.end(), local_ids);
  //std::cout << "... finished list of cells" << std::endl;
};//getListCells

void ZoltanPartitioner::getNumNeighborsCellSlice(void * data, int num_gid_entries, int num_lid_entries, int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int * num_edges, int * ierr){
  //std::cout << "Getting num of neighbors for cell slice..." << std::endl;
  std::pair<Mesh*, Partitioner*> * myPair = (std::pair<Mesh*, Partitioner*> *) data;
  *ierr = ZOLTAN_OK;
  std::vector<int> cell2Face(myPair->first->getReferenceElement()->getNumFaces());
  std::vector<int> face2Cell;
  int numNeighbors;
  int locFace;
  //std::cout << "Cell Number/Number of neighbors:\n";
  //for(int i = 0; i < num_obj; i++){std::cout << global_ids[i] << " ";}
  //std::cout << "\n";
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
    //std::cout << num_edges[i] << " ";
  }
  //std::cout << std::endl;
};//getNumNeighborsCellSlice

void ZoltanPartitioner::getListNeighborsCellSlice(void * data, int num_gid_entries, int num_lid_entries, int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int * num_edges, ZOLTAN_ID_PTR nbor_global_ids, int * nbor_procs, int wgt_dim, float * ewgts, int * ierr){
  //std::cout << "Getting neighbors for cell slice..." << std::endl;
  std::pair<Mesh*, Partitioner*> * myPair = (std::pair<Mesh*, Partitioner*> *) data;
  *ierr = ZOLTAN_OK;
  std::vector<int> cell2Face(myPair->first->getReferenceElement()->getNumFaces());
  std::vector<int> face2Cell;
  int locFace;
  int posNbors = 0;
  //std::cout << "Cell - Neigbors : " << std::endl;
  for(int i = 0; i < num_obj; i++){// for the cells
    //std::cout << global_ids[i] << " - ";
    myPair->first->getCell2Face(local_ids[i], &cell2Face);
    for(int j = 0; j < cell2Face.size(); j++){ //for the faces on the cells
      locFace = myPair->second->global2LocalFace(cell2Face[j]);// get the local index of the face
      if(locFace != -1){//if the face is held by this partition
        myPair->first->getFace2Cell(locFace, &face2Cell);
        if(face2Cell.size() != 1){//if the face is not part of the boundary of the mesh
          nbor_procs[posNbors] = myPair->second->getRank();
          if(face2Cell[0] == global_ids[i]){//if the current cell is the first adjacency
            nbor_global_ids[posNbors] = face2Cell[1];//take the second adjacency
            //std::cout << face2Cell[1] << " ";
          } else if(face2Cell[1] == global_ids[i]){//if the current cell is the second adjacency
            nbor_global_ids[posNbors] = face2Cell[0];//take the first
            //std::cout << face2Cell[0] << " ";
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
            //std::cout << nbor_global_ids[posNbors] << " ";
            posNbors += 1;
            break;
          }
        }
        if(!foundEl){
          throw(ErrorHandle("ZoltanPartitioner", "getListNeighborsCellSlice", "could not find shared face in the sharedFaceList"));
        }
      }
    }
    //std::cout << std::endl;
  }
  //std::cout << "finished getting neighbors for slice" << std::endl;
};//getListNeighborsCellSlice

void ZoltanPartitioner::freeZChanges(){
  zObj->LB_Free_Part(&(zChanges.import_global_ids), &(zChanges.import_local_ids), &(zChanges.import_procs), &(zChanges.import_to_part));
  zObj->LB_Free_Part(&(zChanges.export_global_ids), &(zChanges.export_local_ids), &(zChanges.export_procs), &(zChanges.export_to_part));
  zChanges.isfree = 1;
};//freeZChanges

};//hfox
