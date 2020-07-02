#include "Partitioner.h"

namespace hfox{

void Partitioner::initialize(){
  int mpiInit = 0;
  initializedMPI = 0;
  MPI_Initialized(&mpiInit);
  if(!mpiInit){
    MPI_Init(NULL, NULL);
    initializedMPI = 1;
  }
  MPI_Comm_size(MPI_COMM_WORLD, &nPartitions);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //here we need to implement the initializing of the global ids
  int locNNodes = myMesh->getNumberPoints();
  int locNFaces = myMesh->getNumberFaces();
  int locNCells = myMesh->getNumberCells();
  std::vector<int> nLocAll(nPartitions, 0);
  int globalStartingIndex = 0;
  MPI_Allgather(&locNNodes, 1, MPI_INT, nLocAll.data(), nLocAll.size(), MPI_INT, MPI_COMM_WORLD);
  if(rank != 0){
    globalStartingIndex = std::accumulate(nLocAll.begin(), nLocAll.begin() + rank - 1, 0);
  }
  nodeIDs.resize(locNNodes, 0);
  std::iota(nodeIDs.begin(), nodeIDs.end(), globalStartingIndex);
  MPI_Allgather(&locNFaces, 1, MPI_INT, nLocAll.data(), nLocAll.size(), MPI_INT, MPI_COMM_WORLD);
  if(rank != 0){
    globalStartingIndex = std::accumulate(nLocAll.begin(), nLocAll.begin() + rank - 1, 0);
  }
  faceIDs.resize(locNFaces, 0);
  std::iota(faceIDs.begin(), faceIDs.end(), globalStartingIndex);
  MPI_Allgather(&locNCells, 1, MPI_INT, nLocAll.data(), nLocAll.size(), MPI_INT, MPI_COMM_WORLD);
  if(rank != 0){
    globalStartingIndex = std::accumulate(nLocAll.begin(), nLocAll.begin() + rank - 1, 0);
  }
  elementIDs.resize(locNCells, 0);
  std::iota(elementIDs.begin(), elementIDs.end(), globalStartingIndex);
  computeSharedFaces();
  initialized = 1;
};//initialize

void Partitioner::computeSharedFaces(){
  bool isShared;
  int localCellInd;
  int globalFaceInd;
  std::vector<int> localSharedFaces;
  std::vector<int> cell;
  for(int i = 0; i < myMesh->getNumberFaces(); i++){
    myMesh->getFace2Cell(i, &cell);
    if(cell.size() == 2){//if not boundary face
      for(int j = 0; j < cell.size(); j++){//for adjacent elements
        localCellInd = global2LocalElement(cell[j]);//find local index of element
        if(localCellInd == -1){//if adjacent element not on this partition
          globalFaceInd = local2GlobalFace(i);
          localSharedFaces.push_back(globalFaceInd);
          localSharedFaces.push_back(rank);
          localSharedFaces.push_back(cell[j]);
          break;
        }
      }
    }
  }
  int nSharedFaces = localSharedFaces.size()/3;
  int totNSharedFaces = 0;
  MPI_Allreduce(&nSharedFaces, &totNSharedFaces, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  std::vector<int> allSharedFaces(totNSharedFaces*3, 0);
  MPI_Allgather(localSharedFaces.data(), localSharedFaces.size(), MPI_INT, allSharedFaces.data(), allSharedFaces.size(), MPI_INT, MPI_COMM_WORLD);
  int sharedFaceInd, procRank, adjEl;
  bool found;
  sharedFaceList.resize(nSharedFaces*3, 0);
  for(int i = 0; i < nSharedFaces; i++){
    found = 0;
    sharedFaceInd = localSharedFaces[i*3];
    for(int j = 0; j < totNSharedFaces; j++){
      if(sharedFaceInd == allSharedFaces[j*3]){
        if(allSharedFaces[j*3 + 1] != rank){
          found = 1;
          std::copy(allSharedFaces.begin() + j*3, allSharedFaces.end() + (j+1)*3, sharedFaceList.begin() + i*3);
          break;
        }
      }
    }
    if(!found){
      throw(ErrorHandle("Partitioner", "computeSharedFaces", "could not find the entry for shared face " + std::to_string(sharedFaceInd)));
    }
  }
};//computeSharedFaces

Partitioner::~Partitioner(){
  if(initializedMPI){
    MPI_Finalize();
  }
};//destructor

void Partitioner::setMesh(Mesh * pMesh){
  myMesh = pMesh;
};//setMesh

int Partitioner::getTotalNumberNodes() const{
  int nNodes = myMesh->getNumberPoints();
  int totnNodes = 0;
  MPI_Allreduce(&nNodes, &totnNodes, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return totnNodes;
};//getTotalNumberNodes

int Partitioner::getTotalNumberFaces() const{
  int nFaces = myMesh->getNumberFaces();
  int totnFaces = 0;
  MPI_Allreduce(&nFaces, &totnFaces, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return totnFaces;
};//getTotalNumberFaces

int Partitioner::getTotalNumberEls() const{
  int nCells = myMesh->getNumberCells();
  int totnCells = 0;
  MPI_Allreduce(&nCells, &totnCells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return totnCells;
};//getTotalNumberEls

void Partitioner::local2GlobalNodeSlice(const std::vector<int> & loc, std::vector<int> * glob) const{
  Utils::slice(loc, &nodeIDs, glob);
};//local2GlobalNodeSlice

void Partitioner::local2GlobalFaceSlice(const std::vector<int> & loc, std::vector<int> * glob) const{
  Utils::slice(loc, &faceIDs, glob);
};//local2GlobalFaceSlice

void Partitioner::local2GlobalElementSlice(const std::vector<int> & loc, std::vector<int> * glob) const{
  Utils::slice(loc, &elementIDs, glob);
};//local2GlobalElementSlice

int Partitioner::global2LocalNode(int glob) const{
  std::vector<int>::const_iterator it = std::find(nodeIDs.begin(), nodeIDs.end(), glob);
  int loc = -1;
  if(it != nodeIDs.end()){
    loc = std::distance(nodeIDs.begin(), it);
  }
  return loc;
};//global2LocalNode

int Partitioner::global2LocalFace(int glob) const{
  std::vector<int>::const_iterator it = std::find(faceIDs.begin(), faceIDs.end(), glob);
  int loc = -1;
  if(it != nodeIDs.end()){
    loc = std::distance(faceIDs.begin(), it);
  }
  return loc;
};//global2LocalFace

int Partitioner::global2LocalElement(int glob) const{
  std::vector<int>::const_iterator it = std::find(elementIDs.begin(), elementIDs.end(), glob);
  int loc = -1;
  if(it != nodeIDs.end()){
    loc = std::distance(elementIDs.begin(), it);
  }
  return loc;
};//global2LocalElement

void Partitioner::global2LocalNodeSlice(const std::vector<int> & glob, std::vector<int> * loc) const{
  loc->resize(glob.size());
  for(int i = 0; i < glob.size(); i++){
    loc->at(i) = global2LocalNode(glob[i]);
  }
};//global2LocalNodeSlice

void Partitioner::global2LocalFaceSlice(const std::vector<int> & glob, std::vector<int> * loc) const{
  loc->resize(glob.size());
  for(int i = 0; i < glob.size(); i++){
    loc->at(i) = global2LocalFace(glob[i]);
  }
};//global2LocalFaceSlice

void Partitioner::global2LocalElementSlice(const std::vector<int> & glob, std::vector<int> * loc) const{
  loc->resize(glob.size());
  for(int i = 0; i < glob.size(); i++){
    loc->at(i) = global2LocalElement(glob[i]);
  }
};//global2LocalElementSlice

void Partitioner::multiplyIndexes(int dim, const std::vector<int> * indexes, std::vector<int> * multiIndexes){
  multiIndexes->resize(dim*(indexes->size()));
  for(int i = 0; i < indexes->size(); i++){
    for(int j = 0; j < dim; j++){
      multiIndexes->at(i*dim + j) = indexes->at(i)*dim + j;
    }
  }
};

}
