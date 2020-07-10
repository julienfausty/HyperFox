#include "Partitioner.h"

namespace hfox{

void Partitioner::initialize(){
  int mpiInit = 0;
  MPI_Initialized(&mpiInit);
  if(!mpiInit){
    throw(ErrorHandle("Partitioner", "inititalize", "MPI should be initialized before initializing the partitioner"));
  }
  MPI_Comm_size(MPI_COMM_WORLD, &nPartitions);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int locNNodes = myMesh->getNumberPoints();
  int locNFaces = myMesh->getNumberFaces();
  int locNCells = myMesh->getNumberCells();
  std::vector<int> nLocAll;
  nLocAll.resize(nPartitions);
  int globalStartingIndex = 0;
  MPI_Allgather(&locNNodes, 1, MPI_INT, nLocAll.data(), 1, MPI_INT, MPI_COMM_WORLD);
  if(rank != 0){
    globalStartingIndex = std::accumulate(nLocAll.begin(), nLocAll.begin() + rank, 0);
  }
  nodeIDs.resize(locNNodes, 0);
  std::iota(nodeIDs.begin(), nodeIDs.end(), globalStartingIndex);
  MPI_Allgather(&locNFaces, 1, MPI_INT, nLocAll.data(), 1, MPI_INT, MPI_COMM_WORLD);
  if(rank != 0){
    globalStartingIndex = std::accumulate(nLocAll.begin(), nLocAll.begin() + rank, 0);
  }
  faceIDs.resize(locNFaces, 0);
  std::iota(faceIDs.begin(), faceIDs.end(), globalStartingIndex);
  MPI_Allgather(&locNCells, 1, MPI_INT, nLocAll.data(), 1, MPI_INT, MPI_COMM_WORLD);
  if(rank != 0){
    globalStartingIndex = std::accumulate(nLocAll.begin(), nLocAll.begin() + rank, 0);
  }
  elementIDs.resize(locNCells, 0);
  std::iota(elementIDs.begin(), elementIDs.end(), globalStartingIndex);
  computeSharedFaces();
  initialized = 1;
};//initialize

void Partitioner::computeSharedFaces(){
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
  std::vector<int> nSharedFacesPerProc(nPartitions, 0);
  MPI_Allgather(&nSharedFaces, 1, MPI_INT, nSharedFacesPerProc.data(), 1, MPI_INT, MPI_COMM_WORLD);
  int totNSharedFaces = std::accumulate(nSharedFacesPerProc.begin(), nSharedFacesPerProc.end(), 0);
  std::vector<int> allSharedFaces(totNSharedFaces*3, 0);
  std::vector<int> faceOffsets(nPartitions, 0);
  for(int i = 0; i < nPartitions; i++){
    faceOffsets[i] = std::accumulate(nSharedFacesPerProc.begin(), nSharedFacesPerProc.begin()+i, 0) * 3;
  }
  for(int i = 0; i < nPartitions; i++){
    nSharedFacesPerProc[i] *= 3;
  }
  MPI_Allgatherv(localSharedFaces.data(), localSharedFaces.size(), MPI_INT, allSharedFaces.data(), nSharedFacesPerProc.data(), faceOffsets.data(), MPI_INT, MPI_COMM_WORLD);
  int adjEl;
  sharedFaceList.resize(0);
  for(int i = 0; i < totNSharedFaces; i++){
    adjEl = global2LocalElement(allSharedFaces[i*3+2]);
    if(adjEl != -1){
      for(int k = 0; k < 3; k++){
        sharedFaceList.push_back(allSharedFaces[i*3 + k]);
      }
    }
  }
};//computeSharedFaces

Partitioner::~Partitioner(){
};//destructor

void Partitioner::setMesh(Mesh * pMesh){
  myMesh = pMesh;
};//setMesh

void Partitioner::emptyImportExportMaps(){
  exportMap.clear();
  importMap.clear();
}

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
  //std::cout << "rank: " << rank << "\n";
  //std::cout << "loc Faces: " << nFaces << "\n";
  //std::cout << "glob Faces: " << totnFaces << std::endl;
  return totnFaces;
};//getTotalNumberFaces

int Partitioner::getTotalNumberEls() const{
  int nCells = myMesh->getNumberCells();
  int totnCells = 0;
  MPI_Allreduce(&nCells, &totnCells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  //std::cout << "rank: " << rank << "\n";
  //std::cout << "loc Cells: " << nCells << "\n";
  //std::cout << "glob Cells: " << totnCells << std::endl;
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
  if(it != faceIDs.end()){
    loc = std::distance(faceIDs.begin(), it);
  }
  return loc;
};//global2LocalFace

int Partitioner::global2LocalElement(int glob) const{
  std::vector<int>::const_iterator it = std::find(elementIDs.begin(), elementIDs.end(), glob);
  int loc = -1;
  if(it != elementIDs.end()){
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

void Partitioner::update(){
  if(!initialized){
    throw(ErrorHandle("Partitioner", "update", "must compute initialize the Partitioner and compute the Partition before updating the fields"));
  }
  if(rank == 0){
    std::cout << "Updating mesh and fields with new partition..." << std::endl;
  }
  std::map<int, std::map<FieldType, std::vector<int> > >::iterator itMap;
  std::vector<int> cell;
  std::vector<int> cell2Face;
  int dimSpace = myMesh->getNodeSpaceDimension();
  int nNodesPCell = myMesh->getReferenceElement()->getNumNodes();
  int nNodesPFace = myMesh->getReferenceElement()->getFaceElement()->getNumNodes();
  int nFacesPCell = myMesh->getNumFacesPerCell();
  int localCellIndex;
  int tag;
  //make fieldMap (FieldType, *Field)
  std::map<FieldType, std::vector<Field*> > fieldMap;
  FieldType * fT;
  for(int i = 0; i < fields.size(); i++){
    fT = fields[i]->getFieldType();
    if(fieldMap.find(*fT) == fieldMap.end()){
      fieldMap[*fT] = std::vector<Field*>();
    }
    fieldMap[*fT].push_back(fields[i]);
  }
  if(rank == 0){
    std::cout << "Generating send and recieve buffers for data" << std::endl;
  }
  //generate send/recieve buffers
  std::map<int, std::vector<int> > iSendBuffer, iRecvBuffer;
  std::map<int, std::vector<double> > dSendBuffer, dRecvBuffer;
  int inExportData, inImportData, dnExportData, dnImportData;
  const std::set<int> * bSet = myMesh->getBoundaryFaces();
  std::set<int>::const_iterator bIt;
  std::vector<int> boundaryBuffer;
  for(itMap = exportMap.begin(); itMap != exportMap.end(); itMap++){
    int nCells = itMap->second[Cell].size();
    int nFaces = itMap->second[Face].size();
    int nNodes = itMap->second[Node].size();
    int nExpCells = nCells * nNodesPCell;//the cells
    int nExpCell2Face = nCells * nFacesPCell;//the cell2Faces
    int nExpFaces = nFaces * nNodesPFace;//the faces
    int nExpFace2Cell = nFaces * 2;//the face2Cell
    int nExpNodes = nNodes*dimSpace;
    inExportData = nExpCells + nExpCell2Face + nExpFaces + nExpFace2Cell;
    for(int i = 0; i < nFaces; i++){
      bIt = bSet->find(itMap->second[Face][i]);
      if(bIt != bSet->end()){
        boundaryBuffer.push_back(*bIt);
      }
    }
    inExportData += boundaryBuffer.size();
    iSendBuffer[itMap->first] = std::vector<int>();
    iSendBuffer[itMap->first].resize(inExportData, 0);
    dnExportData = nExpNodes;
    int nExpFCell = 0, nExpFFace = 0, nExpFNode = 0;
    for(int i = 0; i < fieldMap[Cell].size(); i++){
      nExpFCell += nCells * (*(fieldMap[Cell][i]->getNumObjPerEnt())) * (*(fieldMap[Cell][i]->getNumValsPerObj()));
    }
    for(int i = 0; i < fieldMap[Face].size(); i++){
      nExpFFace += nFaces * (*(fieldMap[Face][i]->getNumObjPerEnt())) * (*(fieldMap[Face][i]->getNumValsPerObj()));
    }
    for(int i = 0; i < fieldMap[Node].size(); i++){
      nExpFNode += nNodes * (*(fieldMap[Node][i]->getNumObjPerEnt())) * (*(fieldMap[Node][i]->getNumValsPerObj()));
    }
    dnExportData += nExpFCell + nExpFFace + nExpFNode;
    dSendBuffer[itMap->first] = std::vector<double>();
    dSendBuffer[itMap->first].resize(dnExportData, 0.0);
    int foffset, nFVals;
    for(int i = 0; i < nCells; i++){
      localCellIndex = global2LocalElement(itMap->second[Cell][i]);
      myMesh->getCell(localCellIndex, &cell);
      std::copy(cell.begin(), cell.end(), iSendBuffer[itMap->first].begin() + i*nNodesPCell);
      myMesh->getCell2Face(localCellIndex, &cell2Face);
      std::copy(cell2Face.begin(), cell2Face.end(), iSendBuffer[itMap->first].begin() + nExpCells + i*nFacesPCell);
      foffset = nExpNodes;
      for(int k = 0; k < fieldMap[Cell].size(); k++){
        Field * F = fieldMap[Cell][k];
        nFVals = (*(F->getNumObjPerEnt())) * (*(F->getNumValsPerObj()));
        std::copy(F->getValues()->begin() + localCellIndex*nFVals, F->getValues()->begin() + (localCellIndex+1)*nFVals, dSendBuffer[itMap->first].begin() + foffset + i*nFVals);
        foffset += nCells*nFVals;
      }
    }
    int offset = nExpCells + nExpCell2Face;
    for(int i = 0; i < nFaces; i++){
      localCellIndex = global2LocalFace(itMap->second[Face][i]);
      myMesh->getFace(localCellIndex, &cell);
      std::copy(cell.begin(), cell.end(), iSendBuffer[itMap->first].begin() + offset + i*nNodesPFace);
      std::copy(myMesh->getFace2CellMap()->begin() + localCellIndex*2, myMesh->getFace2CellMap()->begin() + (localCellIndex + 1)*2, iSendBuffer[itMap->first].begin() + offset + nExpFaces + i*2);
      foffset = nExpNodes + nExpFCell;
      for(int k = 0; k < fieldMap[Face].size(); k++){
        Field * F = fieldMap[Face][k];
        nFVals = (*(F->getNumObjPerEnt())) * (*(F->getNumValsPerObj()));
        std::copy(F->getValues()->begin() + localCellIndex*nFVals, F->getValues()->begin() + (localCellIndex+1)*nFVals, dSendBuffer[itMap->first].begin() + foffset + i*nFVals);
        foffset += nFaces*nFVals;
      }
    }
    offset += nExpFaces + nExpFace2Cell;
    std::copy(boundaryBuffer.begin(), boundaryBuffer.end(), iSendBuffer[itMap->first].begin() + offset);
    const std::vector<double> * pNodes = myMesh->getPoints();
    for(int i = 0; i < nNodes; i++){
      localCellIndex = global2LocalNode(itMap->second[Node][i]);
      std::copy(pNodes->begin() + localCellIndex*dimSpace, pNodes->begin() + (localCellIndex+1)*dimSpace, dSendBuffer[itMap->first].begin() + i*dimSpace);
      foffset = nExpNodes + nExpFCell + nExpFFace;
      for(int k = 0; k < fieldMap[Node].size(); k++){
        Field * F = fieldMap[Node][k];
        nFVals = (*(F->getNumObjPerEnt())) * (*(F->getNumValsPerObj()));
        std::copy(F->getValues()->begin() + localCellIndex*nFVals, F->getValues()->begin() + (localCellIndex+1)*nFVals, dSendBuffer[itMap->first].begin() + foffset + i*nFVals);
        foffset += nNodes*nFVals;
      }
    }
  }
  if(rank == 0){
    std::cout << "Communicating send buffers of data" << std::endl;
  }
  //send/recieve the buffers
  for(itMap = exportMap.begin(); itMap != exportMap.end(); itMap++){
    tag = 0;
    int s = iSendBuffer[itMap->first].size();
    MPI_Send(&s, 1, MPI_INT, itMap->first, tag, MPI_COMM_WORLD);
    tag += 1;
    s = dSendBuffer[itMap->first].size();
    MPI_Send(&s, 1, MPI_INT, itMap->first, tag, MPI_COMM_WORLD);
    tag += 1;
    MPI_Send(iSendBuffer[itMap->first].data(), iSendBuffer[itMap->first].size(), MPI_INT, itMap->first, tag, MPI_COMM_WORLD);
    tag += 1;
    MPI_Send(dSendBuffer[itMap->first].data(), dSendBuffer[itMap->first].size(), MPI_DOUBLE, itMap->first, tag, MPI_COMM_WORLD);
  }
  for(itMap = importMap.begin(); itMap != importMap.end(); itMap++){
    iRecvBuffer[itMap->first] = std::vector<int>();
    dRecvBuffer[itMap->first] = std::vector<double>();
    tag = 0;
    int s;
    MPI_Recv(&s, 1, MPI_INT, itMap->first, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    iRecvBuffer[itMap->first].resize(s, 0);
    tag += 1;
    MPI_Recv(&s, 1, MPI_INT, itMap->first, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    dRecvBuffer[itMap->first].resize(s, 0.0);
    tag += 1;
    MPI_Recv(iRecvBuffer[itMap->first].data(), iRecvBuffer[itMap->first].size(), MPI_INT, itMap->first, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    tag += 1;
    MPI_Recv(dRecvBuffer[itMap->first].data(), dRecvBuffer[itMap->first].size(), MPI_DOUBLE, itMap->first, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  if(rank == 0){
    std::cout << "Making modifications to mesh and fields" << std::endl;
  }
  //make modifications to mesh and fields
  Modifier< std::vector<int> > remover;
  remover.setType(REMOVE);
  Modifier< std::vector<int> > iAppender;
  iAppender.setType(APPEND);
  Modifier< std::vector<double> > dAppender;
  dAppender.setType(APPEND);
  //remove sent items
  std::vector<int> locIndexes;
  std::vector<int> indexes2Del;
  int nFVals;
  for(itMap = exportMap.begin(); itMap != exportMap.end(); itMap++){
    //cell associated removal
    locIndexes.resize(itMap->second[Cell].size());
    global2LocalElementSlice(itMap->second[Cell], &locIndexes);
    multiplyIndexes(nNodesPCell, &locIndexes, &indexes2Del);
    remover.setValues(&indexes2Del);
    myMesh->modifyCells(&remover);
    multiplyIndexes(nFacesPCell, &locIndexes, &indexes2Del);
    remover.setValues(&indexes2Del);
    myMesh->modifyCell2FaceMap(&remover);
    for(int k = 0; k < fieldMap[Cell].size(); k++){
      Field * F = fieldMap[Cell][k];
      nFVals = (*(F->getNumObjPerEnt())) * (*(F->getNumValsPerObj()));
      multiplyIndexes(nFVals, &locIndexes, &indexes2Del);
      remover.setValues(&indexes2Del);
      remover.modify(F->getValues());
    }
    remover.setValues(&locIndexes);
    remover.modify(&elementIDs);
    //face associated removal
    locIndexes.resize(itMap->second[Face].size());
    global2LocalFaceSlice(itMap->second[Face], &locIndexes);
    multiplyIndexes(nNodesPFace, &locIndexes, &indexes2Del);
    remover.setValues(&indexes2Del);
    myMesh->modifyFaces(&remover);
    multiplyIndexes(2, &locIndexes, &indexes2Del);
    remover.setValues(&indexes2Del);
    myMesh->modifyFace2CellMap(&remover);
    for(int k = 0; k < fieldMap[Face].size(); k++){
      Field * F = fieldMap[Face][k];
      nFVals = (*(F->getNumObjPerEnt())) * (*(F->getNumValsPerObj()));
      multiplyIndexes(nFVals, &locIndexes, &indexes2Del);
      remover.setValues(&indexes2Del);
      remover.modify(F->getValues());
    }
    std::set<int> * bFaces = myMesh->getBoundaryFaces();
    std::set<int>::iterator bIt;
    for(int i = 0; i < itMap->second[Face].size(); i++){
      bIt = bFaces->find(itMap->second[Face][i]);
      if(bIt != bFaces->end()){
        bFaces->erase(bIt);
      }
    }
    remover.setValues(&locIndexes);
    remover.modify(&faceIDs);
    //node associated removal
    locIndexes.resize(itMap->second[Node].size());
    global2LocalNodeSlice(itMap->second[Node], &locIndexes);
    multiplyIndexes(dimSpace, &locIndexes, &indexes2Del);
    remover.setValues(&indexes2Del);
    myMesh->modifyNodes(&remover);
    for(int k = 0; k < fieldMap[Node].size(); k++){
      Field * F = fieldMap[Node][k];
      nFVals = (*(F->getNumObjPerEnt())) * (*(F->getNumValsPerObj()));
      multiplyIndexes(nFVals, &locIndexes, &indexes2Del);
      remover.setValues(&indexes2Del);
      remover.modify(F->getValues());
    }
    remover.setValues(&locIndexes);
    remover.modify(&nodeIDs);
    myMesh->update();
  }
  std::vector<int> iBuffData;
  std::vector<double> dBuffData;
  int foffset;
  for(itMap = importMap.begin(); itMap != importMap.end(); itMap++){
    //populating cell data
    int nCells = itMap->second[Cell].size();
    iAppender.setValues(&(itMap->second[Cell]));
    iAppender.modify(&elementIDs);
    int nImpCells = nCells * nNodesPCell;
    iBuffData.resize(nImpCells);
    std::copy(iRecvBuffer[itMap->first].begin(), iRecvBuffer[itMap->first].begin() + nImpCells, iBuffData.begin());
    iAppender.setValues(&iBuffData);
    myMesh->modifyCells(&iAppender);
    int nImpCell2Faces = nCells * nFacesPCell;
    iBuffData.resize(nImpCell2Faces);
    std::copy(iRecvBuffer[itMap->first].begin() + nImpCells, iRecvBuffer[itMap->first].begin() + nImpCells + nImpCell2Faces, iBuffData.begin());
    iAppender.setValues(&iBuffData);
    myMesh->modifyCell2FaceMap(&iAppender);
    int nNodes = itMap->second[Node].size();
    foffset = nNodes * dimSpace;
    for(int k = 0; k < fieldMap[Cell].size(); k++){
      Field * F = fieldMap[Cell][k];
      nFVals = (*(F->getNumObjPerEnt())) * (*(F->getNumValsPerObj()));
      int nField = nFVals * nCells;
      dBuffData.resize(nField);
      std::copy(dRecvBuffer[itMap->first].begin() + foffset, dRecvBuffer[itMap->first].begin() + foffset + nField, dBuffData.begin());
      dAppender.setValues(&dBuffData);
      dAppender.modify(F->getValues());
      foffset += nField;
    }
    //populating face data
    int nFaces = itMap->second[Face].size();
    iAppender.setValues(&(itMap->second[Face]));
    iAppender.modify(&faceIDs);
    int nImpFaces = nFaces * nNodesPFace;
    iBuffData.resize(nImpFaces);
    int offset = nImpCells + nImpCell2Faces;
    std::copy(iRecvBuffer[itMap->first].begin() + offset, iRecvBuffer[itMap->first].begin() + offset + nImpFaces, iBuffData.begin());
    iAppender.setValues(&iBuffData);
    myMesh->modifyFaces(&iAppender);
    int nImpFace2Cells = nFaces * 2;
    iBuffData.resize(nImpFace2Cells);
    offset += nImpFaces;
    std::copy(iRecvBuffer[itMap->first].begin() + offset, iRecvBuffer[itMap->first].begin() + offset + nImpFace2Cells, iBuffData.begin());
    iAppender.setValues(&iBuffData);
    myMesh->modifyFace2CellMap(&iAppender);
    offset += nImpFace2Cells;
    for(int k = 0; k < fieldMap[Face].size(); k++){
      Field * F = fieldMap[Face][k];
      nFVals = (*(F->getNumObjPerEnt())) * (*(F->getNumValsPerObj()));
      int nField = nFVals * nFaces;
      dBuffData.resize(nField);
      std::copy(dRecvBuffer[itMap->first].begin() + foffset, dRecvBuffer[itMap->first].begin() + foffset + nField, dBuffData.begin());
      dAppender.setValues(&dBuffData);
      dAppender.modify(F->getValues());
      foffset += nField;
    }
    //populating node data
    iAppender.setValues(&(itMap->second[Node]));
    iAppender.modify(&nodeIDs);
    int nImpNodes = nNodes * dimSpace;
    dBuffData.resize(nImpNodes);
    std::copy(dRecvBuffer[itMap->first].begin(), dRecvBuffer[itMap->first].begin() + nImpNodes, dBuffData.begin());
    dAppender.setValues(&dBuffData);
    myMesh->modifyNodes(&dAppender);
    for(int k = 0; k < fieldMap[Node].size(); k++){
      Field * F = fieldMap[Node][k];
      nFVals = (*(F->getNumObjPerEnt())) * (*(F->getNumValsPerObj()));
      int nField = nFVals * nNodes;
      dBuffData.resize(nField);
      std::copy(dRecvBuffer[itMap->first].begin() + foffset, dRecvBuffer[itMap->first].begin() + foffset + nField, dBuffData.begin());
      dAppender.setValues(&dBuffData);
      dAppender.modify(F->getValues());
      foffset += nField;
    }
    myMesh->update();
  }
  emptyImportExportMaps();
  //compute the shared faces list
  computeSharedFaces();
  //update the mesh/fields with shared information
  if(rank == 0){
    std::cout << "... finished updating partition" << std::endl;
  }
};//update


void Partitioner::multiplyIndexes(int dim, const std::vector<int> * indexes, std::vector<int> * multiIndexes){
  multiIndexes->resize(dim*(indexes->size()));
  for(int i = 0; i < indexes->size(); i++){
    for(int j = 0; j < dim; j++){
      multiIndexes->at(i*dim + j) = indexes->at(i)*dim + j;
    }
  }
};

}
