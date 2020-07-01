#include "ZoltanPartitioner.h"

namespace hfox{

ZoltanPartitioner::~ZoltanPartitioner(){
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
};//initialize


void ZoltanPartitioner::computePartition(){
  zObj->Set_Param("OBJ_WEIGHT_DIM", "0");
  zObj->Set_Num_Obj_Fn(ZoltanPartitioner::getNumberCells, myMesh);
  zObj->Set_Obj_List_Fn(ZoltanPartitioner::getListCells, &elementIDs);
  std::pair<Mesh*, Partitioner*> thePair(myMesh, this);
  zObj->Set_Num_Edges_Multi_Fn(ZoltanPartitioner::getNumNeighborsCellSlice, &thePair);
  zObj->Set_Edge_List_Multi_Fn(ZoltanPartitioner::getListNeighborsCellSlice, &thePair);
};//computePartition

void ZoltanPartitioner::update(){
};//update

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

};//hfox
