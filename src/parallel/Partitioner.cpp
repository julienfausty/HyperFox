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
  initialized = 1;
};//initialize

Partitioner::~Partitioner(){
  if(initializedMPI){
    MPI_Finalize();
  }
}

int Partitioner::getTotalNumberNodes() const{
  int nNodes = myMesh->getNumberPoints();
  int totnNodes = 0;
  MPI_Allreduce(&nNodes, &totnNodes, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return totnNodes;
}

int Partitioner::getTotalNumberFaces() const{
  int nFaces = myMesh->getNumberFaces();
  int totnFaces = 0;
  MPI_Allreduce(&nFaces, &totnFaces, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return totnFaces;
}

int Partitioner::getTotalNumberEls() const{
  int nCells = myMesh->getNumberCells();
  int totnCells = 0;
  MPI_Allreduce(&nCells, &totnCells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return totnCells;
}

void Partitioner::local2GlobalNodeSlice(const std::vector<int> & loc, std::vector<int> * glob) const{
  Utils::slice(loc, &nodeIDs, glob);
};

void Partitioner::local2GlobalFaceSlice(const std::vector<int> & loc, std::vector<int> * glob) const{
  Utils::slice(loc, &faceIDs, glob);
};

void Partitioner::local2GlobalElementSlice(const std::vector<int> & loc, std::vector<int> * glob) const{
  Utils::slice(loc, &elementIDs, glob);
};

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
